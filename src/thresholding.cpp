#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export("rcpp_mat_sort")]]
mat rcpp_mat_sort(mat M){
  int r = M.n_rows, c=M.n_cols;
  mat f(r,c);
  
  for(int i=0;i<r; ++i){
    rowvec v = M.row(i);
    std::sort(v.begin(),v.end(),std::greater<double>());
    f.row(i) = v;
  }
  
  return f;
}

// [[Rcpp::export("rcpp_weight_thresholding")]]
mat rcpp_weight_thresholding(mat M, double w=0, std::string method="trim"){
  double r = M.n_rows, c = M.n_cols;

  if(method=="trim"){
    arma::uvec ids = find(M < w); 
    M.elem(ids).fill(0);
  }else if(method=="bin"){
    arma::uvec idslo = find(M < w); 
    arma::uvec idsup = find(M >= w); 
    M.elem(idslo).fill(0);
    M.elem(idsup).fill(1);    
  }else if(method=="resh"){
    for(int i=0; i <r; ++i){
      double eps = 0.01;
      rowvec v = M.row(i);
      double maxx = max(v), minx = w-eps;

      arma::uvec idslo = find(v < w);
      arma::uvec idsup = find(v >= w); 
      vec vresh = (v.elem(idsup)-minx)/(maxx-minx);
      v.elem(idsup) = vresh; 
      v.elem(idslo).fill(0);
      M.row(i) = v;
    }
  }    
  return M;
}

// [[Rcpp::export("rcpp_density_thresholding")]]
mat rcpp_density_thresholding(mat M, double w=0, std::string method="trim"){
  double r = M.n_rows, c = M.n_cols;
  mat Mnew(r,c);
  vec wvec(c);

  // Obtain weights for density-based approach
  mat Msort = rcpp_mat_sort(M);
  int d = std::round(w*c);
  if(d > 0){
    wvec = Msort.col(d-1);
  }else{
    wvec = Msort.col(0);
  }
  
  // Apply weights
  if(method=="trim"){
    for(int i=0; i<r; ++i){
      rowvec rowi = M.row(i);
      arma::uvec ids = find(rowi < wvec(i));
      rowi.elem(ids).fill(0);
      Mnew.row(i) = rowi;
    }
  }else if(method=="bin"){
    for(int i=0; i<r; ++i){
      rowvec rowi = M.row(i);
      arma::uvec idslo = find(rowi < wvec(i));
      arma::uvec idsup = find(rowi >= wvec(i));
      rowi.elem(idslo).fill(0);
      rowi.elem(idsup).fill(1);
      Mnew.row(i) = rowi;
    }
  }else if(method=="resh"){
    for(int i=0; i <r; ++i){
      double eps = 0.01;
      rowvec rowi = M.row(i);
      double maxx = max(rowi), minx = wvec(i)-eps;

      arma::uvec idslo = find(rowi < wvec(i));
      arma::uvec idsup = find(rowi >= wvec(i));
      vec vresh = (rowi.elem(idsup)-minx)/(maxx-minx);
      rowi.elem(idsup) = vresh;
      rowi.elem(idslo).fill(0);
      Mnew.row(i) = rowi;
    }
  }
  return Mnew;
}



// [[Rcpp::export("rcpp_vec_to_mat")]]
arma::mat rcpp_vec_to_mat(rowvec v, int p){
  
  mat M(p,p);
  // make empty matrices
  mat Z(p,p,fill::zeros);
  mat X(p,p,fill::zeros);
  
  // fill matrices with integers
  vec idx = linspace<mat>(1,p,p);
  X.each_col() += idx;
  Z.each_row() += trans(idx);
  
  // assign upper triangular elements
  // the >= allows inclusion of diagonal elements
  M.elem(find(Z>X)) = v;
  M  = M.t();
  M.elem(find(Z>X)) = v;
//  M.elem(find(Z<X)) = v;
  
  return M;
}

//[[Rcpp::export("rcpp_cc_func")]]
double rcpp_cc_func(rowvec v, int p, bool weighted=false){
  
  mat A = rcpp_vec_to_mat(v, p);

  int r = A.n_rows;
  vec n(r), plainsum(r), squaresum(r);

  if(!weighted){
    arma::uvec ids = find(A != 0); 
    A.elem(ids).fill(1);    
  }

  for(int i=0; i <r; ++i){
    rowvec v = A.row(i);
    rowvec C = v * A;
    n(i) = dot(C,v.t());

    rowvec v2 = pow(v,2);
    plainsum(i) = sum(v);
    squaresum(i) = sum(v2);
  }
  vec te = (pow(plainsum,2) - squaresum);
  vec cci = n / te;
  for(int j=0; j<cci.n_elem; ++j){
    if(!is_finite(cci(j))){cci(j)=0;}
  }
  return mean(cci);
}

//[[Rcpp::export("rcpp_wrapper_thresholding")]]
mat rcpp_wrapper_thresholding(mat M, double w, int p){
  int r = M.n_rows;

  vec tseq = linspace(0,100,2)/100;
  std::string tmeth = "trim";
  
  mat f(r, 2);
  mat wt_mat = rcpp_weight_thresholding(M=M, w=w, tmeth);
  mat dt_mat = rcpp_density_thresholding(M=M, w=w, tmeth);
    
  for(int i=0; i<r; ++i){
    rowvec v(2);
    rowvec rowi_wt = wt_mat.row(i);
    rowvec rowi_dt = dt_mat.row(i);
    
    v[0] = rcpp_cc_func(rowi_wt, p);
    v[1] = rcpp_cc_func(rowi_dt, p);
    
    f.row(i) = v;
    }

  return f;
}

/*** R
#set.seed(222)
x <- matrix(abs(rnorm(25, mean = 0, sd = 0.5)),10,10)
x[lower.tri(x)] = t(x)[lower.tri(x)]
diag(x)<-0
x[x>1] <-1
print(x)

rcpp_mat_sort(x)
rcpp_density_thresholding(x, w=0.5, method="trim")
rcpp_density_thresholding(x, w=0.5, method="bin")
rcpp_density_thresholding(x, w=0.5, method="resh")


# rcpp_cc_func(y)
# z=rcpp_thresholding(x, w=0.5, method="bin")
# mean(WGCNA::clusterCoef(z))

*/
