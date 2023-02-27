#' Generalized Kernel Regularized Least Squares
#'
#' This page documents how to estimate \code{gKRLS} using \code{mgcv}'s
#' functions, e.g. \code{bam} or \code{gam}. \code{gKRLS} can be specified as
#' shown in the accompanying examples. Post-estimation functions to calculate
#' marginal effects are documented elsewhere, e.g. \link{calculate_effects}. 
#' 
#' @details 
#' 
#' The \code{gKRLS} function should not be called directly and is a control
#' argument to the smoother in \code{mgcv}, i.e. \code{s(..., bs = "gKRLS", xt =
#' gKRLS(...)}. Its arguments are described above Multiple kernels can be
#' included alongside other smooth arguments specified via \code{s(...)}.
#'
#' \bold{Note:} Variables must be separated with commas inside of \code{s(...)}.
#'
#' @encoding UTF-8
#' @param demean_kernel A logical value that indicates whether columns of the
#'   (sketched) kernel should be demeaned before estimation. The default is
#'   \code{FALSE}.
#' @param sketch_method A string that specifies which kernel sketching method
#'   should be used. Options include \code{"subsampling"} (sub-sampling),
#'   \code{"gaussian"}, \code{"bernoulli"}, or \code{"none"} (no sketching).
#'   Default is \code{"subsampling"}. See Drineas et al. (2005) and Yang et al.
#'   (2017) for details.
#' @param standardize A string that specifies how the data is standardized
#'   before distance between observations is calculated. The default is
#'   \code{"Mahalanobis"}. Other options are \code{"scaled"} (ensure all
#'   non-constant columns are mean zero and variance one) or \code{"none"} (no
#'   standardization).
#' @param sketch_multiplier By default, sketching size increases with \code{c *
#'   ceiling(nrow(X)^(1/3))} where \code{c} is the "multiplier". Default of 5;
#'   if results seems unstable, Chang and Goplerud (2023) find that 15 works
#'   well. See \code{sketch_size_raw} to directly set the sketching size.
#' @param sketch_size_raw Set the exact sketching size (independent of N).
#'   Exactly one of this or \code{sketch_multiplier} must be \code{NULL}.
#' @param sketch_prob For Bernoulli sketching, this sets the probability of
#'   \code{1}. See Yang et al. (2017) for details on this method.
#' @param rescale_penalty A logical value for whether the penalty should be
#'   rescaled for numerical stability. See documentation for
#'   \code{mgcv::smooth.spec} on the meaning of this term. The default is
#'   \code{TRUE}.
#' @param bandwidth The bandwidth for the kernel \eqn{D} where each element of
#'   the kernel is defined by \eqn{\exp(-||x_i - x_j||^2_2/D)}.
#' @param remove_instability A logical value that indicates whether numerical
#'   zeros (set via \code{truncate.eigen.tol}) should be removed when building
#'   the penalty matrix. The default is \code{TRUE}.
#' @param truncate.eigen.tol Remove columns of the penalty, i.e. \code{S^T K S},
#'   whose eigenvalue is below \code{truncate.eigen.tol}. This ensures a
#'   numerically positive-definite penalty. These columns are also removed from
#'   the sketched kernel. Default is \code{sqrt(.Machine$double.eps)}. Setting
#'   to 0 retains all numerically non-negative eigenvalues. This helps with the
#'   numerical stability of the algorithm. Removal can be disabled using
#'   \code{remove_instability = FALSE}.
#' @references 
#' 
#' Chang, Qing and Max Goplerud. 2023. "Generalized Kernel Regularized Least
#' Squares". \url{https://arxiv.org/abs/2209.14355}.
#' 
#' Drineas, Petros and Mahoney, Michael W and Nello Cristianini. 2005. "On the
#' Nyström Method for Approximating a Gram Matrix For Improved Kernel-Based
#' Learning". \emph{Journal of Machine Learning Research} 6(12):2153-2175.
#' 
#' Yang, Yun and Pilanci, Mert and Martin J. Wainwright. 2017. "Randomized
#' Sketches for Kernels: Fast and Optimal Nonparametric Regression".
#' \emph{Annals of Statistics} 45(3):991-1023.
#' 
#' @useDynLib gKRLS
#' @import Matrix
#' @export
#'
#' @examples
#' n <- 100
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' state <- sample(letters[1:5], n, replace = TRUE)
#' y <- 0.3 * x1 + 0.4 * x2 + 0.5 * x3 + rnorm(n)
#' data <- data.frame(y, x1, x2, x3, state)
#' data$state <- factor(data$state)
#' # A gKRLS model without fixed effects
#' gkrls_est <- mgcv::gam(y ~ s(x1, x2, x3, bs = "gKRLS"), data = data)
#' summary(gkrls_est)
#' # A gKRLS model with fixed effects
#' gkrls_fx <- mgcv::gam(y ~ state + s(x1, x2, x3, bs = "gKRLS"), data = data)
#' # Change default standardization to Mahalanobis, sketch method to Gaussian,
#' # and alter sketching multiplier
#' gkrls_mah <- mgcv::gam(y ~ s(x1, x2, x3,
#'   bs = "gKRLS",
#'   xt = gKRLS(
#'     standardize = "Mahalanobis",
#'     sketch_method = "gaussian",
#'     sketch_multiplier = 2
#'   )
#' ),
#' data = data
#' )
#' # A model with multiple kernels
#' gkrls_multiple <- mgcv::gam(y ~ s(x1, x2, bs = 'gKRLS') + s(x1, x3, bs = 'gKRLS'), data = data)
#' # calculate marginal effect
#' calculate_effects(gkrls_est, variables = "x1", continuous_type = "derivative")
gKRLS <- function(truncate.eigen.tol = sqrt(.Machine$double.eps),
                  demean_kernel = FALSE,
                  sketch_method = "subsampling",
                  standardize = "Mahalanobis",
                  sketch_multiplier = 5,
                  sketch_size_raw = NULL,
                  sketch_prob = NULL, 
                  rescale_penalty = TRUE,
                  bandwidth = NULL, 
                  remove_instability = TRUE) {
  if (length(sketch_method) == 1){
    sketch_method <- match.arg(sketch_method, c("subsampling", "gaussian", "bernoulli", "none"))
  }else{
    sketch_vector <- sketch_method
    sketch_method <- 'custom'
  }
  standardize <- match.arg(standardize, c("Mahalanobis", "scaled", "none"))
  if (!(rescale_penalty %in% c(TRUE, FALSE))){
    stop('rescale_penalty must be TRUE or FALSE.')
  }
  if (!is.null(sketch_size_raw) & !is.null(sketch_multiplier)) {
    stop('Only one of "sketch_size_raw" or "sketch_multiplier" may be provided.')
  } else if (is.null(sketch_size_raw) & is.null(sketch_multiplier)) {
    stop('One of "sketch_size_raw" or "sketch_multiplier" must be provided.')
  }

  return(mget(ls()))
}