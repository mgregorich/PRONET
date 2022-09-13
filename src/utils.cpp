#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export("cpp_mat_sort")]]
mat cpp_mat_sort(mat M){
  int r = M.n_rows, c=M.n_cols;
  mat f(r,c);
  
  for(int i=0;i<r; ++i){
    rowvec v = M.row(i);
    std::sort(v.begin(),v.end(),std::greater<double>());
    f.row(i) = v;
  }
  
  return f;
}

// [[Rcpp::export("cpp_weight_thresholding")]]
mat cpp_weight_thresholding(mat M, double w=0, std::string method="trim"){
  int r = M.n_rows;

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

// [[Rcpp::export("cpp_density_thresholding")]]
mat cpp_density_thresholding(mat M, double w=0, std::string method="trim"){
  int r = M.n_rows, c = M.n_cols;
  mat Mnew(r,c);
  vec wvec(c);

  // Obtain weights for density-based approach
  mat Msort = cpp_mat_sort(M);
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


// [[Rcpp::export("cpp_vec_to_mat")]]
arma::mat cpp_vec_to_mat(rowvec x, int p){
  
  arma::mat out(p, p, arma::fill::zeros);
  arma::uvec lw_idx = arma::trimatl_ind( arma::size(out) , -1);
  out.elem(lw_idx) = x;
  out  = out.t();
  out.elem(lw_idx) = x;
  
  return out;
}

//[[Rcpp::export("cpp_cc_func")]]
double cpp_cc_func(rowvec v, int p, bool weighted=false){
  
  mat A = cpp_vec_to_mat(v, p);

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
  int lcci = cci.n_elem;
  for(int j=0; j<lcci; ++j){
    if(!is_finite(cci(j))){cci(j)=0;}
  }
  double cc = mean(cci);
  return cc;
}

//[[Rcpp::export("cpp_wrapper_thresholding")]]
Rcpp::List cpp_wrapper_thresholding(mat M, int p, double step_size){
  int r = M.n_rows;

  vec tseq = linspace(0,1,1/step_size+1);
  std::string tmeth = "trim";

  Rcpp::List out(tseq.n_elem);
  for(int j=0; j<tseq.n_elem;++j){
    mat cc(r, 2);
    double w = tseq(j);
    
    mat wt_mat = cpp_weight_thresholding(M=M, w, tmeth);
    mat dt_mat = cpp_density_thresholding(M=M, w, tmeth);

    for(int i=0; i<r; ++i){
      rowvec v(2);
      rowvec rowi_wt = wt_mat.row(i);
      rowvec rowi_dt = dt_mat.row(i);

      v[0] = cpp_cc_func(rowi_wt, p);
      v[1] = cpp_cc_func(rowi_dt, p);

      cc.row(i) = v;
    }
    // mat joined_mats = join_rows(tseq, cc);
    out(j) = cc.as_col();
  }
  return out;
}

//[[Rcpp::export("cpp_transform_to_beta")]]
Rcpp::NumericVector cpp_transform_to_beta(Rcpp::NumericVector eta, vec beta_pars, vec eta_pars){
  Rcpp::NumericVector p = Rcpp::pnorm(eta, eta_pars[0], eta_pars[1]);
  Rcpp::NumericVector q = Rcpp::qbeta(p, beta_pars[0], beta_pars[1]);
  return q;
}

/*** R
# set.seed(222)
# x <- matrix(abs(rnorm(25, mean = 0, sd = 0.5)),10,10)
# x[lower.tri(x)] = t(x)[lower.tri(x)]
# diag(x)<-0
# x[x>1] <-1
# print(x)
# v <- c(1:6)
# p <- 4
# cpp_vec_to_mat_new(v, p)
# 
# cpp_mat_sort(x)
# cpp_weight_thresholding(x, w=0.5, method="trim")
# cpp_density_thresholding(x, w=0.5, method="bin")
# cpp_density_thresholding(x, w=0.5, method="resh")
# 
# cpp_cc_func(as.numeric(GE.thres[1,]), p=p)
# z=cpp_thresholding(x, w=0.5, method="bin")
# mean(WGCNA::clusterCoef(z))
# test <- cpp_wrapper_thresholding(M=as.matrix(data.network), p=p)

# cpp_transform_to_beta(etai, beta_pars=c(2,5), eta_pars = c(5.5, 3.5))

*/
