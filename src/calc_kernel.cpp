#include "RcppEigen.h"
#include "calc_kernel.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
double kern_gauss(
    const Eigen::MatrixXd X_one,
    const Eigen::MatrixXd X_two,
    const double bandwidth
){
  double dist = (X_one - X_two).squaredNorm();
  double out = exp(-dist/bandwidth);
  return out;
}