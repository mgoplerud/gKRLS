#include "RcppEigen.h"
#include <algorithm>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// K-means++ initialization using vectorized Eigen distance computation
Eigen::MatrixXd kmeans_pp_init(
    const Eigen::MatrixXd& X,
    int k
) {
  int N = X.rows();
  int P = X.cols();
  Eigen::MatrixXd centers(k, P);

  // First center: random point
  int first = (int)(R::runif(0, 1) * N);
  if (first >= N) first = N - 1;
  centers.row(0) = X.row(first);

  // Squared distances to nearest center so far
  Eigen::VectorXd min_dist = (X.rowwise() - centers.row(0)).rowwise().squaredNorm();

  for (int c = 1; c < k; c++) {
    // Sample proportional to squared distance (D² weighting)
    double total = min_dist.sum();
    if (total <= 0) {
      // Degenerate case: pick randomly
      int idx = (int)(R::runif(0, 1) * N);
      if (idx >= N) idx = N - 1;
      centers.row(c) = X.row(idx);
    } else {
      double r = R::runif(0, 1) * total;
      double cumsum = 0;
      int chosen = N - 1;
      for (int i = 0; i < N; i++) {
        cumsum += min_dist(i);
        if (cumsum >= r) {
          chosen = i;
          break;
        }
      }
      centers.row(c) = X.row(chosen);
    }

    // Update min distances with new center
    Eigen::VectorXd new_dist = (X.rowwise() - centers.row(c)).rowwise().squaredNorm();
    min_dist = min_dist.cwiseMin(new_dist);
  }

  return centers;
}


//' Mini-batch k-means clustering for Nystrom landmark selection
//'
//' Uses k-means++ initialization and vectorized batch distance computation
//' via Eigen for fast landmark point selection.
//'
//' @param X Data matrix (N x P)
//' @param k Number of clusters (landmark points)
//' @param batch_size Size of each mini-batch
//' @param max_iter Maximum number of iterations
//' @param tol Convergence tolerance on centroid shift
//' @return Matrix of k centroids (k x P)
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd minibatch_kmeans_cpp(
    const Eigen::Map<Eigen::MatrixXd> X,
    int k,
    int batch_size,
    int max_iter,
    double tol
) {
  int N = X.rows();
  int P = X.cols();

  if (k <= 0) {
    Rcpp::stop("k must be positive");
  }
  if (k >= N) {
    return Eigen::MatrixXd(X);
  }

  int bs = std::min(batch_size, N);

  // K-means++ initialization
  Eigen::MatrixXd centers = kmeans_pp_init(X, k);
  Eigen::VectorXi counts = Eigen::VectorXi::Ones(k);

  // Pre-compute center squared norms for fast distance computation
  Eigen::VectorXd c_norms = centers.rowwise().squaredNorm();

  Eigen::VectorXi batch_idx(bs);
  Eigen::MatrixXd X_batch(bs, P);

  for (int iter = 0; iter < max_iter; iter++) {
    Eigen::MatrixXd old_centers = centers;

    // Sample mini-batch
    for (int i = 0; i < bs; i++) {
      batch_idx(i) = (int)(R::runif(0, 1) * N);
      if (batch_idx(i) >= N) batch_idx(i) = N - 1;
      X_batch.row(i) = X.row(batch_idx(i));
    }

    // Vectorized squared distance: D_ij = ||x_i||^2 + ||c_j||^2 - 2 x_i . c_j
    // This is the hot path; the matrix multiply is BLAS-optimized
    Eigen::VectorXd x_norms = X_batch.rowwise().squaredNorm();
    Eigen::MatrixXd D = -2.0 * X_batch * centers.transpose();
    D.colwise() += x_norms;
    D.rowwise() += c_norms.transpose();

    // Assign each batch point to nearest center
    Eigen::VectorXi assignments(bs);
    for (int i = 0; i < bs; i++) {
      Eigen::Index min_idx;
      D.row(i).minCoeff(&min_idx);
      assignments(i) = (int)min_idx;
    }

    // Accumulate per-cluster sums and counts, then do batch update
    Eigen::MatrixXd cluster_sums = Eigen::MatrixXd::Zero(k, P);
    Eigen::VectorXi n_assigned = Eigen::VectorXi::Zero(k);

    for (int i = 0; i < bs; i++) {
      int c = assignments(i);
      cluster_sums.row(c) += X_batch.row(i);
      n_assigned(c)++;
    }

    for (int c = 0; c < k; c++) {
      if (n_assigned(c) > 0) {
        counts(c) += n_assigned(c);
        double lr = (double)n_assigned(c) / counts(c);
        Eigen::RowVectorXd batch_mean = cluster_sums.row(c) / n_assigned(c);
        centers.row(c) += lr * (batch_mean - centers.row(c));
      }
    }

    // Update cached center norms
    c_norms = centers.rowwise().squaredNorm();

    // Check convergence: total squared shift of centroids
    double shift = (centers - old_centers).squaredNorm();
    if (shift < tol * tol * k) break;
  }

  return centers;
}
