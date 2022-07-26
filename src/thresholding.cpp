//sum.cpp
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins("cpp11")]]
#include <cmath>
#include <iostream>
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <iterator>
using namespace arma;
using namespace Rcpp;
using namespace Eigen;




// [[Rcpp::export]]
MatrixXd rcpp_DensityThresholding(MatrixXd result, double d=0, std::string method="trim"){
  
  int r= result.rows();
  int c = result.cols();
  double w;
  
  result.transposeInPlace();
  VectorXd v = (Map<VectorXd>(result.data(), r*c));
  std::sort(v.data(),v.data()+v.size()); 
  v = v.transpose();
  
  int s = v.size();
  
  if(s*d==0){
    result = (result.array() < w).select(0, result.unaryExpr([](double x) {
      x=0;
      return x;
    }));
  }else{
    int frac = round(s-s*d);
    w = v[frac];
    
    if(method=="trim"){
      result = (result.array() < w).select(0, result);
    }else if(method=="bin"){
      result = (result.array() < w).select(0, result.unaryExpr([](double x) {
        x=1;
        return x;
      }));
    }else if(method=="resh"){
      double maxx = result.maxCoeff();
      double minx = w;
      
      result = (result.array() < w).select(0, result.unaryExpr([&maxx, &minx](double x) {
        double y;
        y=(x - minx) / (maxx - minx);
        return y;
      }));
    }
  }

  // Return the vector to R
  return result;
}

// [[Rcpp::export]]
MatrixXd rcpp_WeightThresholding(MatrixXd result, double w=0, std::string method="trim"){
  
  // Creating a vector object
  if(method=="trim"){
    result = (result.array() < w).select(0, result);
  }else if(method=="bin"){
    result = (result.array() < w).select(0, result.unaryExpr([](double d) {
      d=1;
      return d;
    }));
  }else if(method=="resh"){
    double maxx = result.maxCoeff();
    double minx = w;

    result = (result.array() < w).select(0, result.unaryExpr([&maxx, &minx](double x) {
      double y;
      y=(x - minx) / (maxx - minx);
      return y;
    }));
  }
  
  // Return the vector to R
  return result;
  }



/*** R
# set.seed(222)
# x <- matrix(abs(rnorm(25, mean = 0, sd = 0.5)),5,5)
# print(x)
# 
# rcpp_WeightThresholding(x, w=0.2, method = "trim")
# rcpp_WeightThresholding(x, w=0.2, method = "resh")
# rcpp_WeightThresholding(x, w=0.2, method = "bin")
# 
# rcpp_DensityThresholding(x, d=0.5, method="bin")



*/
