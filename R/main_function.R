#' gKRLS
#'
#' Fit a generalized KRLS model using mgcv
#'
#' The function `gKRLS` should be given as a control argument to any relevant
#' function for fitting a generalized additive model from `mgcv` (e.g. `bam`,
#' `gam`, `gamm4`). Its arguments are described below with a simple example.
#'
#' A kernel can be specified by `s(..., bs = "gKRLS")`.
#'
#' NOTE: Variables must be separated by commas inside of the `s(..)`.
#'
#' @param demean_kernel A logical variable ``True'' or ``False'' indicates whether
#' the kernel should be demeaned. The default is False. If True is given, column mean
#' will be used to calculate the demean kernel matrix.
#' @param sketch_method String vector that specifies which kernel sketch method
#'   should be used. Options include ``gaussian'': gaussian kernel approximation,
#'    ``nystrom'':  Nystrom approximation, ``bernoulli'', and no sketch. The default is
#'   Nystrom approximation.
#' @param standardize A string vector that specifies which standardization
#'   method should be used. Must be one of ``scaled'', ``Mahalanobis'', or ``none'', which
#'   menas no standardization. The default is ``scaled''.
#' @param sketch_multiplier By default, sketching size increases with c ``ceiling(nrow(X)^(1/3)''
#'   where c is the "multiplier". By default, set to 5; if results seem
#'   unstable, try increasing to around 15.
#' @param sketch_size_raw If desired, set the exact sketching size (independent
#'   of N). Exactly one of this or sketch_multiplier must be NULL.
#' @param sketch_prob For bernoulli sketching, what is probability of "1"?
#' @param rescale_penalty Rescale penalty for numerical stability; see documentation for
#'   \code{mgcv::smooth.spec} on the meaning of this term. Default of "TRUE".
#' @param remove_instability A logical variable indicates whether numerical
#'   zeros (set via truncate.eigen.tol) should be removed when building the
#'   penalty matrix. The default is ``true''.
#' @param truncate.eigen.tol It determines how much eigenvalues should be truncated.
#' If truncate.eigen.tol set to 1e-6, this means we only keep eigenvalues greater or
#' equal to 1e-6 in the penalty. The default is sqrt(.Machine$double.eps).
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
#'
#' # calculate marginal effect
#' calculate_effects(gkrls_est, variables = "x1", continuous_type = "derivative")
gKRLS <- function(truncate.eigen.tol = sqrt(.Machine$double.eps),
                  demean_kernel = FALSE,
                  sketch_method = "nystrom",
                  standardize = "Mahalanobis",
                  sketch_multiplier = 5,
                  sketch_size_raw = NULL,
                  sketch_prob = NULL, 
                  rescale_penalty = TRUE,
                  remove_instability = TRUE) {
  sketch_method <- match.arg(sketch_method, c("nystrom", "gaussian", "bernoulli", "none"))
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
