#' Marginal Effects by Derivative
#'
#' @description Find the (average) marginal effect and standard error by taking the
#' analytical derivative of the conditional expectation function. This method is
#' designed to provide comparability with the original \code{KRLS} package, but
#' is rather limited in the scenarios where it can be applied. Specifically, it can
#' only compute effects with a single kernel and where variables are not both in
#' the kernel and outside the kernel.
#'
#' Because of these restrictions, users should rely on \code{calculate_effects}
#' and this function may be removed in future updates. If a derivative is really
#' desired, this can be approximated as described in the
#' \link{calculate_effects} documentation.
#'
#' @param object A fitted gKRLS model.
#' @param newdata A new data frame that used to calculate the marginal effect.
#' @param keep A character string indicating which variables to calculate marginal
#' effects. Default is \code{NULL}, which calculates marginal effects for all.
#'
#' @return The \code{legacy_marginal_effect} a list that contains the following elements:
#' \itemize{
#' \item{"ME_pointwise": } The marginal effects for each observation.
#' \item{"ME_pointwise_var": } The variance for each pointwise maringal effect
#' in "ME_pointwise".
#' \item{"AME_pointwise": } The average marginal effect, i.e. the column
#' averages of "ME_pointwise".
#' \item{"AME_pointwise_var": } The variance of each average marginal effect.
#' }
#' 
#' @examples
#' n <- 500
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' state <- sample(letters, n, replace = TRUE)
#' y <- 0.3 * x1 + 0.4 * x2 + 0.5 * x3 + rnorm(n)
#' data <- data.frame(y, x1, x2, x3, state)
#'
#' # A gKRLS model
#' gkrls_est <- mgcv::gam(y ~ s(x1, x2, x3, bs = "gKRLS"), data = data)
#' # calculate marginal effect using derivative
#' legacy_marginal_effect(gkrls_est, newdata = data)
#' @importFrom stats plogis dnorm pnorm predict coef
#' @export
legacy_marginal_effect <- function(object, newdata, keep = NULL) {
  
  standardize <- object$internal$standardize

  N_eff <- length(object$y) - sum(object$edf)
  N <- length(object$y)

  full_lp <- predict(object, newdata = newdata, type = "lpmatrix")
  model_offset <- attr(full_lp, "model.offset")

  class_smooth <- lapply(object$smooth, FUN = function(i) {
    class(i)
  })
  stopifnot(all(lengths(class_smooth) == 2))
  class_smooth <- sapply(class_smooth, FUN = function(i) {
    i[1]
  })

  if (sum(class_smooth == "gKRLS.smooth") != 1) {
    stop("marginal_effects can only be run with one kernel...for now...")
  }
  kern_smooth <- object$smooth[[which(class_smooth == "gKRLS.smooth")]]

  # Get the fixed effects
  fe_id <- setdiff(
    seq_len(ncol(full_lp)),
    unlist(lapply(object$smooth, FUN = function(i) {
      i$first.par:i$last.par
    }))
  )
  if (length(fe_id) > 0) {
    FE_matrix_test <- full_lp[, fe_id, drop = F]
  } else {
    FE_matrix_test <- matrix(nrow = nrow(full_lp), ncol = 0)
  }
  colnames(FE_matrix_test) <- colnames(full_lp)[fe_id]

  kernel_id <- unlist(kern_smooth[c("first.para", "last.para")])
  kernel_id <- kernel_id[1]:kernel_id[2]
  # All *other* smoothes, e.g. random effects
  Z <- full_lp[, -c(fe_id, kernel_id), drop = F]

  standardize <- kern_smooth$xt$standardize

  std_whiten <- kern_smooth$std_train$whiten
  std_mean <- kern_smooth$std_train$mean
  W_Matrix <- kern_smooth$std_train$W_Matrix

  std_X_train <- kern_smooth$X_train
  # WX_train <- sweep(std_X_train %*% MASS::ginv(as.matrix(kern_smooth$std_train$whiten)), 2,
  #     kern_smooth$std_train$mean, "+") %*% W_Matrix
  WX_train <- sweep(
    std_X_train %*% t(kern_smooth$std_train$whiten), 2,
    kern_smooth$std_train$mean %*% W_Matrix, "+"
  )

  kern_smooth$xt$return_raw <- TRUE
  kern_smooth$std_train <- NULL

  raw_X_test <- PredictMat(kern_smooth, data = newdata)
  WX_test <- raw_X_test %*% W_Matrix

  family <- object$family

  fe_names <- colnames(FE_matrix_test)
  names_mfx <- c(fe_names, kern_smooth$term)
  
  if (!identical(kern_smooth$term, colnames(kern_smooth$X_train))) {
    if (standardize != "Mahalanobis") {
      message("Names of kern and std_X_train seem misaligned...")
      message('Names from mgcv')
      message(paste0(kern_smooth$term, collapse=', '))
      message('Names from std_x_train')
      message(paste0(colnames(std_X_train), collapse=', '))
      stop('Error in legacy_predict. This usually occurs when factors are provided to mgcv.\n"calculate_effects" should work as expected.')
    }else{
      if (any(!sapply(kern_smooth$term_levels, is.null))){
        stop('Error in legacy_predict. This usually occurs when factors are provided to mgcv.\n"calculate_effects" should work as expected.')
      }
    }
  }


  std_X_test <- sweep(raw_X_test, 2, std_mean, FUN = "-")
  std_X_test <- as.matrix(std_X_test %*% std_whiten)

  tS <- t(kern_smooth$sketch_matrix)
  bandwidth <- kern_smooth$bandwidth

  all_mean <- matrix(coef(object))
  if (length(fe_id) > 0) {
    fe_mean <- matrix(coef(object)[fe_id])
    re_mean <- all_mean[-fe_id, , drop = F]
  } else {
    re_mean <- all_mean
    fe_mean <- numeric(0)
  }
  vcov_ridge <- vcov(object)
  fmt_sd_y <- 1

  N_train <- nrow(std_X_train)
  N_test <- nrow(std_X_test)

  if (length(model_offset) == 1) {
    offset <- rep(0, N_test)
  } else {
    if (length(model_offset) != N_test) {
      stop("Model offest is incorrect length.")
    }
    offset <- model_offset
  }

  if (ncol(Z) == 0) {
    any_Z <- FALSE
    Z <- sparseMatrix(i = 1, j = 1, x = 0)
  } else {
    any_Z <- TRUE
    offset_mean <- matrix(coef(object)[-c(fe_id, kernel_id)])
    re_mean <- matrix(coef(object)[kernel_id])
    offset <- as.vector(offset + Z %*% offset_mean)
    Z <- drop0(Z)
  }
  Sc <- t(tS) %*% re_mean

  fd_flag <- kern_smooth$fd_flag

  if (family$family == "binomial" & family$link == "logit") {
    family <- "logit"
  } else if (family$family == "binomial" & family$link == "probit") {
    family <- "probit"
  } else if (family$family == "gaussian") {
    family <- "gaussian"
  } else if (family$family == "poisson") {
    family <- "poisson"
  } else {
    stop('Invalid family for legacy; use "calculate_effects"!')
  }


  ME_pointwise <- ME_pointwise_var <- matrix(
    NA, nrow(std_X_test),
    ncol(FE_matrix_test) + ncol(W_Matrix)
  )
  AME_grad <- matrix(0, nrow = nrow(all_mean), ncol = ncol(FE_matrix_test) + ncol(W_Matrix))
  SIZE_FE <- ncol(FE_matrix_test)
  SIZE_KERNEL <- ncol(W_Matrix)

  # Flag the columns that should be analyzed using first differences
  fd_flag <- names(fd_flag)
  fd_matrix <- matrix(data = 0, ncol = 2, nrow = SIZE_KERNEL + SIZE_FE)
  rownames(fd_matrix) <- names_mfx
  # The values to set the first difference to
  fd_matrix[fd_flag, 1] <- 0
  fd_matrix[fd_flag, 2] <- 1
  fd_flag <- (rownames(fd_matrix) %in% fd_flag)

  std_fd_matrix <- sweep(t(fd_matrix), 2, c(rep(0, SIZE_FE), std_mean), FUN = "-")
  std_fd_matrix <- std_fd_matrix %*%
    bdiag(Diagonal(x = rep(0, SIZE_FE)), std_whiten)
  std_fd_matrix <- t(std_fd_matrix)
  std_fd_matrix <- as.matrix(std_fd_matrix)

  if (standardize == "Mahalanobis") {
    std_fd_matrix[, ] <- 0
  } else {
    rownames(std_fd_matrix) <- names_mfx
  }

  type_mfx <- c(rep("deriv_FE", ncol(FE_matrix_test)), rep("deriv_Kern", nrow(std_whiten)))
  type_mfx[which(fd_flag)] <- "FD"
  mahal <- standardize == "Mahalanobis"
  std_whiten <- as.matrix(std_whiten)
  vcov_ridge <- as.matrix(vcov_ridge)
  W_Matrix <- as.matrix(W_Matrix)
  WX_test <- as.matrix(WX_test)
  WX_train <- as.matrix(WX_train)

  if (!is.null(keep)) {
    # Variables to calculate AME for
    kern_to_fit <- intersect(kern_smooth$term, keep)
    fe_to_fit <- intersect(fe_names, keep)
    missing_provided <- setdiff(setdiff(keep, kern_to_fit), fe_to_fit)
    all_to_fit <- c(fe_to_fit, kern_to_fit)
    if (length(missing_provided) != 0) {
      warning('Some variables in "keep" not located in data.')
    }
    fit_position <- which(rownames(fd_matrix) %in% all_to_fit)
    mfx_counter <- seq_len(length(fit_position))
  } else {
    fit_position <- seq_len(nrow(fd_matrix))
    mfx_counter <- seq_len(nrow(fd_matrix))
  }

  cpp_fit_position <- as.integer(fit_position - 1)
  cpp_mfx_counter <- as.integer(mfx_counter - 1)

  out_ME <- cpp_gkrls_me(
    std_X_train = std_X_train, std_X_test = std_X_test,
    bandwidth = bandwidth, family = family,
    tZ = t(Z), offset = offset, any_Z = any_Z,
    mahal = mahal, sd_y = fmt_sd_y, tS = tS, fe_mean = fe_mean,
    re_mean = as.vector(re_mean),
    SIZE_PARAMETER = nrow(all_mean), vcov_ridge = vcov_ridge,
    FE_matrix_test = FE_matrix_test,
    W_Matrix = W_Matrix, WX_test = WX_test, WX_train = WX_train,
    raw_X_test = raw_X_test,
    std_mean = std_mean, std_whiten = std_whiten,
    fd_matrix = fd_matrix, std_fd_matrix = std_fd_matrix,
    type_mfx = type_mfx, fit_position = cpp_fit_position,
    mfx_counter = cpp_mfx_counter
  )

  names(out_ME$AME_pointwise) <- colnames(out_ME$ME_pointwise) <- colnames(out_ME$ME_pointwise_var) <- names_mfx[fit_position]
  colnames(out_ME$AME_grad) <- names(out_ME$AME_pointwise_var) <- names_mfx[fit_position]

  out_ME$type <- type_mfx[fit_position[mfx_counter]]
  out_ME$mfx_type <- "legacy"
  out_ME$N_eff <- N_eff
  out_ME$N <- N

  class(out_ME) <- c("gKRLS_mfx")

  return(out_ME)
}
