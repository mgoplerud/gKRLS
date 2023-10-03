#include "RcppEigen.h"
#include "calc_kernel.h"
#include <cmath>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

//' Create the sketched kernel
//' @param X_test Test data
//' @param X_train Train data
//' @param S Sketch matrix
//' @param bandwidth Kernel bandwidth
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd create_sketched_kernel(
    const Eigen::Map<Eigen::MatrixXd> X_test,
    const Eigen::Map<Eigen::MatrixXd> X_train,
    const Eigen::MatrixXd S,
    const double bandwidth,
    const bool raw = false
){
  
  int N_test = X_test.rows();
  int N_train = X_train.rows();
  int M = S.rows();
  
  if (N_train != S.cols()){
    Rcpp::stop("ncol(S) must equal nrow(X_train)");
  }
  
  Eigen::MatrixXd KSt(N_test, M);
  
  // Loop over every row of the test data.
  for (int i = 0; i < N_test; i++){
    // Get the N_train \times 1 K(x_i, x_j) for x_j \in {1, \cdots, N_train\}
    Eigen::VectorXd kern_i(N_train);
    Eigen::VectorXd X_i = X_test.row(i);
    
    for (int j = 0; j < N_train; j++){
      Eigen::VectorXd X_j = X_train.row(j);
      kern_i(j) = kern_gauss(X_i, X_j, bandwidth, raw);      
    }
    
    // Get the sketched kernel S 
    KSt.row(i) = S * kern_i;
  }
  return KSt;
}