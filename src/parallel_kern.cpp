// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppParallel.h>

using namespace std;
using namespace Rcpp;
using namespace RcppParallel;
using namespace RcppEigen;

// https://github.com/RcppCore/RcppParallelTests/blob/master/testthat/cpp/distance.cpp
template <typename InputIterator1, typename InputIterator2>
inline double kernel_iter(InputIterator1 begin1, InputIterator1 end1, 
                       InputIterator2 begin2, const double bandwidth) {
  
  // value to return
  double rval = 0;
  
  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  
  // for each input item
  while (it1 != end1) {
    
    // take the value and increment the iterator
    double d1 = *it1++;
    double d2 = *it2++;
    rval += std::pow(d1 - d2, 2);
  }
  rval = std::exp(-rval / bandwidth);
  return rval;
}

// Helpful guidance from https://gist.github.com/JWiley/d9cba55603471f75d438
// on a worked example of RcppParallel

struct KernBuild : public Worker {

  // input
  const RMatrix<double> X_test;
  const RMatrix<double> X_train;
  const Eigen::MatrixXd& tS;
  const double& bandwidth;
  const unsigned int& M;
  const unsigned int& N_train;
  // output
  RMatrix<double> KS;
  
  // Key Function
  KernBuild(const NumericMatrix& X_test, const NumericMatrix& X_train, 
       const Eigen::MatrixXd& tS, const double& bandwidth,
       Rcpp::NumericMatrix KS)
    : X_test(X_test), X_train(X_train), 
      tS(tS), bandwidth(bandwidth), 
      M(tS.rows()), N_train(X_train.nrow()),
      KS(KS) {}

  void operator()(std::size_t begin, std::size_t end) {
    
    for(unsigned int i = begin; i < end; ++i){
      // Get X_i
      RMatrix<double>::Row X_i = X_test.row(i);
      Eigen::VectorXd k_i(N_train);
      for (unsigned int j = 0; j < N_train; j++){
        // Get X_j
        RMatrix<double>::Row X_j = X_train.row(j);
        double k_ij = kernel_iter(X_i.begin(), X_i.end(), X_j.begin(), bandwidth);
        k_i(j) = k_ij;
      }
      Eigen::VectorXd sketch_i = tS * k_i;
      for (unsigned int m = 0; m < M; m++){
        KS(i,m) = sketch_i(m);
      }
    }
  }
};

// [[Rcpp::export]]
Rcpp::NumericMatrix build_kern_parallel(
    const NumericMatrix X_train, 
    const NumericMatrix X_test,
    const Eigen::MatrixXd tS,
    double bandwidth,
    size_t grain_size = 100,
    int threads = -1) {
    
    Rcpp::NumericMatrix KS(X_test.nrow(), tS.rows());
  
    // pass input and output
    
    KernBuild kbuild(X_train, X_test, tS, bandwidth, KS);
    
#if RCPP_PARALLEL_USE_TBB
    Rcout << "Using " << threads << " threads to estimate sketched kernel." << std::endl;
#else
    Rcout << "Not using TBB; threads not respected" << std::endl;
    // parallelFor to do it
#endif
    
    parallelFor(0, X_test.nrow(), kbuild, grain_size, threads);
    
    return KS;
  }
  