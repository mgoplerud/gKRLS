#include "RcppEigen.h"
#include "calc_kernel.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
double kern_gauss(
    const Eigen::VectorXd X_one,
    const Eigen::VectorXd X_two,
    const double bandwidth,
    const bool raw
){
  double dist = (X_one - X_two).squaredNorm();
  double out;
  if (raw){
    out = dist;
  }else{
    out = std::exp(-dist/bandwidth);
  }
  return out;
}