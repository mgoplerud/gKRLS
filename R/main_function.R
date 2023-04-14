#' Generalized Kernel Regularized Least Squares
#'
#' This page documents how to use \code{gKRLS} as part of a model estimated with
#' \code{mgcv}. Post-estimation functions to calculate marginal effects are
#' documented elsewhere, e.g. \link{calculate_effects}.
#' 
#' @details 
#' 
#' \bold{Overview:} The \code{gKRLS} function should not be called directly. Its
#' options, described above, control how \code{gKRLS} is estimated. It should be
#' passed to \code{mgcv} as follows: \code{s(x1, x2, x3, bs = "gKRLS", xt =
#' gKRLS(...))}. Multiple kernels can be specified and have different
#' \code{gKRLS} arguments. It can also be used alongside the existing options
#' for \code{s()} in \code{mgcv}.
#' 
#' \bold{Default Settings:} By default, \code{bs = "gKRLS"} uses Mahalanobis
#' distance between the observations, random sketching using subsampling
#' sketching (i.e., where the kernel is constructed using a random sample of the
#' observations; Yang et al. 2017) and a sketching dimension of \code{5 *
#' ceiling(N^(1/3))} where \code{N} is the number of observations. Chang and
#' Goplerud (2023) provide an exploration of alternative options.
#' 
#' \bold{Notes:} Please note that variables must be separated with commas inside
#' of \code{s(...)} and that character variables should usually be passed as
#' factors to work smoothly with \code{mgcv}. When using this function with
#' \code{bam}, the sketching dimension uses \code{chunk.size} in place of
#' \code{N} and thus either \code{chunk.size} or \code{sketch_size_raw} must be used to cause
#' the sketching dimension to increase with \code{N}.
#' 
#' @encoding UTF-8
#' @param demean_kernel A logical value that indicates whether columns of the
#'   (sketched) kernel should be demeaned before estimation. The default is
#'   \code{FALSE}.
#' @param sketch_method A string that specifies which kernel sketching method
#'   should be used (default of \code{"subsampling"}). Options include
#'   \code{"subsampling"}, \code{"gaussian"}, \code{"bernoulli"}, or
#'   \code{"none"} (no sketching). Drineas et al. (2005) and Yang et al. (2017)
#'   provide more details on these options.
#'
#'   To force \code{"subsampling"} to select a specific set of observations, you
#'   can provide a vector of row positions to \code{sketch_method}. This
#'   manually sets the size of the sketching multiplier, implicitly overriding
#'   other options in \code{gKRLS}. The examples provide an illustration.
#' @param standardize A string that specifies how the data is standardized
#'   before calculating the distance between observations. The default is
#'   \code{"Mahalanobis"} (i.e., demeaned and transformed to have an identity
#'   covariance matrix). Other options are \code{"scaled"} (all columns are
#'   scaled to have mean zero and variance of one) or \code{"none"} (no
#'   standardization).
#' @param sketch_multiplier A number that sets the size of the sketching
#'   dimension: \code{sketch_multiplier * ceiling(N^(1/3))} where \code{N} is
#'   the number of observations. The default is 5; Chang and Goplerud (2023)
#'   find that increasing this to 15 may improve stability for certain complex
#'   kernels. \code{sketch_size_raw} can directly set the size of the sketching
#'   dimension.
#' @param sketch_size_raw A number to set the exact size of the sketching
#'   dimension. The default, \code{NULL}, means that this argument is not used
#'   and the size depends on the number of observations; see
#'   \code{sketch_multiplier}. Exactly one of \code{sketch_size_raw} or
#'   \code{sketch_multiplier} must be \code{NULL}.
#' @param sketch_prob A probability for an element of the sketching matrix to
#'   equal \code{1} when using Bernoulli sketching. Yang et al. (2017) provide
#'   more details.
#' @param rescale_penalty A logical value for whether the penalty should be
#'   rescaled for numerical stability. See documentation for
#'   \code{mgcv::smooth.spec} on the meaning of this term. The default is
#'   \code{TRUE}.
#' @param bandwidth A bandwidth \eqn{P} for the kernel where each element of
#'   the kernel \eqn{(i,j)} is defined by \eqn{\exp(-||x_i - x_j||^2_2/P)}.
#' @param remove_instability A logical value that indicates whether numerical
#'   zeros (set via \code{truncate.eigen.tol}) should be removed when building
#'   the penalty matrix. The default is \code{TRUE}.
#' @param truncate.eigen.tol A threshold to remove columns of the penalty
#'   \eqn{S K S^T} whose eigenvalues are small (below
#'   \code{truncate.eigen.tol}). These columns are removed from the sketched
#'   kernel and avoids instability due to numerically very small eigenvalues. The
#'   default is \code{sqrt(.Machine$double.eps)}. This adjustment can be
#'   disabled by setting \code{remove_instability = FALSE}.
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
#' @returns \code{gKRLS} returns a named list with the elements in "Arguments".
#' 
#' @useDynLib gKRLS
#' @import Matrix
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' state <- sample(letters[1:5], n, replace = TRUE)
#' y <- 0.3 * x1 + 0.4 * x2 + 0.5 * x3 + rnorm(n)
#' data <- data.frame(y, x1, x2, x3, state)
#' data$state <- factor(data$state)
#' # A gKRLS model without fixed effects
#' fit_gKRLS <- mgcv::gam(y ~ s(x1, x2, x3, bs = "gKRLS"), data = data)
#' summary(fit_gKRLS)
#' # A gKRLS model with fixed effects outside of the kernel
#' fit_gKRLS_FE <- mgcv::gam(y ~ state + s(x1, x2, x3, bs = "gKRLS"), data = data)
#' 
#' # HC3 is not available for mgcv; this uses the effective degrees of freedom
#' # instead of the number of columns; see ?estfun.gam for details
#' robust <- sandwich::vcovHC(fit_gKRLS, type = 'HC1')
#' cluster <- sandwich::vcovCL(fit_gKRLS, cluster = data$state)
#' 
#' # Change default standardization to "scaled", sketch method to Gaussian,
#' # and alter sketching multiplier
#' fit_gKRLS_alt <- mgcv::gam(y ~ s(x1, x2, x3,
#'   bs = "gKRLS",
#'   xt = gKRLS(
#'     standardize = "scaled",
#'     sketch_method = "gaussian",
#'     sketch_multiplier = 2
#'   )
#' ),
#' data = data
#' )
#' # A model with multiple kernels
#' fit_gKRLS_2 <- mgcv::gam(y ~ s(x1, x2, bs = 'gKRLS') + s(x1, x3, bs = 'gKRLS'), data = data)
#' # A model with a custom set of ids for sketching
#' id <- sample(1:n, 5)
#' fit_gKRLS_custom <- mgcv::gam(y ~ s(x1, bs = 'gKRLS', xt = gKRLS(sketch_method = id)), data = data)
#' # Note that the ids of the sampled observations can be extracted 
#' # from the fitted mgcv object
#' stopifnot(identical(id, fit_gKRLS_custom$smooth[[1]]$subsampling_id))
#' # calculate marginal effect (see ?calculate_effects for more examples)
#' calculate_effects(fit_gKRLS, variables = "x1")
gKRLS <- function(sketch_method = "subsampling",
                  standardize = "Mahalanobis",
                  bandwidth = NULL, 
                  sketch_multiplier = 5,
                  sketch_size_raw = NULL,
                  sketch_prob = NULL, 
                  rescale_penalty = TRUE,
                  truncate.eigen.tol = sqrt(.Machine$double.eps),
                  demean_kernel = FALSE,
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