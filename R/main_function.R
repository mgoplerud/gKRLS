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
#' @param no.rescale  
#' @param standardize A string vector that specifies which standardization
#'   method should be used. Must be one of ``scaled'', ``Mahalanobis'', or ``none'', which
#'   menas no standardization. The default is ``scaled''.
#' @param sketch_size Specify the dimension of the randomized sketches of the
#'   kernel matrix. The default size is ``ceiling(nrow(X)^(1/3)) * 5,'' where X
#'   is kernel matrix. If desires to use a larger matrix, change 5 to a larger number, 
#'   e.g. 15.
#' @param remove_instability A logical variable indicates whether 0 should be removed 
#' from the eigenvector when building the kernel matrix. The default is ``True''
#' @param truncate.eigen.tol It determines how much eigenvalues should be truncated. 
#' If truncate.eigen.tol set to 1e-6, this means we only keep eigenvalue greater or 
#' equal to 1e-6. The default is sqrt(.Machine$double.eps), where .Machine$double.eps 
#' is the smallest positive floating-point number x such that 1 + x != 1.
#' @useDynLib gKRLS
#' @import Matrix
#' @export
#' 
#' @examples
#' n <- 5000
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' state <- sample(letters, n, replace = T)
#' y = 0.3*x1 + 0.4*x2 +0.5*x3 + rnorm(n)
#' data <- data.frame(y, x1, x2, x3, state)
#' 
#' # A gKRLS model without fixed effects
#' gkrls_est <- gam(y ~ s(x1,x2,x3, bs="gKRLS"), data = data)
#' summary(gkrls_est)
#' # A gKRLS model with fixed effects
#' gkrls_fx <- gam(y ~ factor(state) + s(x1,x2,x3, bs="gKRLS"), data = data)
#' # Change default standardization to Mahalanobis, sketch method to Nystrom, and larger kernel matrix
#' gkrls_mah <- gam(y ~ s(x1,x2,x3, bs="gKRLS", xt = gKRLS(standardize = 'Mahalanobis',
#'                                                          sketch_method = "nystrom",
#'                                                          sketch_size = function(N) {ceiling(N^(1/3)) * 15 })),
#'                                                           data = data)
#' 
#'  # calculate marginal effect
#'  calculate_effects(gkrls_est, variables = "x1", continuous_type = 'derivative')
#' 
gKRLS <- function(truncate.eigen.tol = sqrt(.Machine$double.eps),
                              demean_kernel = FALSE,
                              sketch_method = 'nystrom',
                              no.rescale = FALSE, standardize = 'mahalanobis',
                              sketch_size = function(N){ceiling(N^(1/3)) * 5},
                              remove_instability = TRUE){
  return(mget(ls()))
}

