
#' @importFrom Matrix sparseMatrix
create_data_gKRLS <- function(term_levels, term_class, data_term, terms, allow.missing.levels = FALSE){

  out <- mapply(term_levels, term_class, data_term, terms, SIMPLIFY = FALSE, 
    FUN=function(l_i, c_i, d_i, n_i){
    if (is.null(l_i)){
      if (!(c_i %in% c('numeric', 'integer'))){
        warning(paste0(n_i, ' is neither numeric nor integer but is not a factor. Check that it is parsed correctly.'))
      }
      d_i <- as.numeric(d_i)
      if (any(is.na(d_i))){
        stop(paste0('Coercing ', n_i, ' to numeric induced NAs.'))
      }
      return(d_i)
    }else{
      if (!any(c_i %in% 'factor')){
        warning(paste0(n_i, ' has levels but does not have a factor class. Check that it is parsed correctly or convert to factor.'))
      }
      d_i_j <- match(d_i, l_i)
      d_i_i <- seq_len(length(d_i))
      if (allow.missing.levels){
        new_levels_ij <- which(is.na(d_i_j))
        if (length(new_levels_ij) > 0){
          d_i_j <- d_i_j[-new_levels_ij]
          d_i_i <- d_i_i[-new_levels_ij]
          message(paste0('New levels found in factor ', n_i,'. See warnings from mgcv.\nGiving a zero value in the kernel.'))
        }
      }else{
        if (any(is.na(d_i_j))){
          stop(paste0('Missing values found in ', n_i,' when looking at levels.'))
        }
      }
      transf_d_i <- as.matrix(sparseMatrix(i = d_i_i, 
                                           j = d_i_j, 
                                           x = 1, 
                                           dims = c(length(d_i), length(l_i))
      ))
      colnames(transf_d_i) <- paste0(n_i, l_i)
      return(transf_d_i)
    }
  })
  out <- do.call('cbind', out)
  return(out)
  
}

#' Constructor for gKRLS smooth
#' @import mgcv
#' @importFrom Rcpp sourceCpp
#' @keywords internal
#' 
#' See \link{gKRLS} for details.
#' @param object a smooth object; see documentation for other methods in
#'   \code{mgcv}.
#' @param data a data.frame; see documentation for other methods in
#'   \code{mgcv}.
#' @param knots not used
#' @export
smooth.construct.gKRLS.smooth.spec <- function(object, data, knots) {
  
  if (is.null(object$xt)) {
    object$xt <- gKRLS()
  }
  if (!("return_raw" %in% names(object$xt))) {
    object$xt$return_raw <- FALSE
  }

  if (!is.null(object$pc)) {
    stop('non-NULL "pc" not set up yet.')
  }
  if (!is.null(object$id)) {
    stop('custom "id" not set up yet for "kern".')
  }
  if (!is.na(object$p.order)) {
    stop('m not used in "kern".')
  }

  if (object$bs.dim != "-1") {
    stop('k should not be modified directly. Set sketch size via "xt"')
  }

  term_levels <- lapply(data[object$term], levels)
  term_class <- lapply(data[object$term], class)
  
  if (length(object$term) > 1) {
    length_data <- sapply(data, length)
    if (length(unique(length_data)) != 1) {
      stop("Error: All input variables should have same length. The syntax is s(a,b,c,d)")
    }
  }
  
  X <- create_data_gKRLS(term_levels = term_levels, 
    term_class = term_class, data_term = data[object$term],
    terms = object$term)
    

  if (NCOL(X) == 1 & is.null(ncol(X))) {
    X <- matrix(X)
  }

  fd_flag <- which(apply(X, MARGIN = 2, FUN = function(i) {
    all(i %in% c(0, 1))
  }))

  std_X <- standardize_design(kernel_X = X, standardize = object$xt$standardize)
  X <- std_X$std_kernel_X
  std_train <- std_X$std_train
  std_X_zero <- std_X$zero_columns
  rm(std_X)
  gc()

  bandwidth <- object$xt$bandwith
  if (is.null(bandwidth)) {
    bandwidth <- ncol(X)
  }

  # Create the Kernel
  N <- nrow(X)

  sketch_size_raw <- object$xt$sketch_size_raw
  sketch_multiplier <- object$xt$sketch_multiplier
  # Either multiply N^(1/3) by multiplier *or* use raw value
  if (!is.null(sketch_multiplier)) {
    sketch_size <- floor(ceiling(N^(1 / 3)) * sketch_multiplier)
  } else {
    sketch_size <- sketch_size_raw
  }


  X_train <- X

  if (is.na(sketch_size)) {
    sketch_matrix <- diag(N)
  } else {
    if (sketch_size >= N) {
      warning("Sketch size exceeds size of data")
    }
    if (sketch_size < 0) { # Do N - s
      sketch_size <- N + sketch_size
    }

    if (object$xt$sketch_method == "nystrom") {
      if (sketch_size > N) {
        stop("Nystrom requires sketch_size < N.")
      }
      nystrom_id <- sample(1:N, sketch_size)
      X_train <- X[nystrom_id, , drop = F]
      sketch_matrix <- diag(length(nystrom_id))
    } else {
      sketch_matrix <- create_sketch_matrix(N, sketch_size, object$xt$sketch_prob, object$xt$sketch_method)
    }
  }
  
  KS <- create_sketched_kernel(
    X_test = X,
    X_train = X_train, tS = t(sketch_matrix), bandwidth = bandwidth
  )
  
  if (!object$fixed) {
    if (!is.na(sketch_size) & object$xt$sketch_method == "nystrom") {
      Penalty <- t(KS[nystrom_id, ]) %*% sketch_matrix
    } else {
      Penalty <- t(KS) %*% sketch_matrix
    }

    # Penalty <- Penalty/sum(Penalty^2)

    if (object$xt$remove_instability) {
      old_size_S <- ncol(KS)

      truncate_eigen_penalty <- object$xt$truncate.eigen.tol

      eigen_sketch <- eigen(Penalty, symmetric = TRUE)
      eigen_sketch$values <- ifelse(eigen_sketch$values < truncate_eigen_penalty, 0, eigen_sketch$values)
      nonzero <- which(eigen_sketch$values != 0)

      if (length(nonzero) == 0) {
        stop('After truncation, all eigenvectors are removed. Decrease "truncate.eigen.tol" or set "remove_instability" = FALSE to proceed.')
      }
      KS <- as.matrix(KS %*% (eigen_sketch$vectors[, nonzero] %*%
        Diagonal(x = 1 / sqrt(eigen_sketch$values[nonzero]))))
      sketch_matrix <- sketch_matrix %*% as.matrix(eigen_sketch$vectors[, nonzero] %*%
        Diagonal(x = 1 / sqrt(eigen_sketch$values[nonzero])))
      Penalty <- diag(ncol(KS))
      if (length(nonzero) != length(eigen_sketch$values)) {
        message(paste("Sketch dimension decreased from", old_size_S, "to", ncol(KS), "because of singular penalty matrix."))
      }
    }

    object$S <- list()
    object$S[[1]] <- Penalty
  } else {
    object$S <- NULL
    warning("fx=TRUE may result in severe overfitting. Use with caution.")
  }

  if (object$xt$demean_kernel) {
    KS_mean <- colMeans(KS)
    KS <- apply(KS, MARGIN = 2, FUN = function(i) {
      i - mean(i)
    })
  } else {
    KS_mean <- rep(0, ncol(KS))
  }

  # Required elements for kernel prediction
  object$KS_mean <- KS_mean
  object$X <- KS
  object$copy <- KS
  object$X_train <- X_train
  object$bandwidth <- bandwidth
  object$sketch_matrix <- sketch_matrix
  object$std_train <- std_train
  object$fd_flag <- fd_flag
  object$term_levels <- term_levels
  object$term_class <- term_class
  
  # Required elements for "gam"
  object$rank <- ncol(KS)
  object$null.space.dim <- 0
  object$df <- ncol(KS)
  # If rescale_penalty is NOT true, then set "no.rescale" to TRUE, i.e. to not rescale.
  # otherwise, leave null.
  if (!object$xt$rescale_penalty) {
    object$no.rescale <- TRUE
  }
  object$te.ok <- 0
  object$plot.me <- FALSE
  object$C <- matrix(nrow = 0, ncol = ncol(KS))

  class(object) <- "gKRLS.smooth" # Give object a class
  object
}

#' Predict Methods for gKRLS smooth
#' @keywords internal
#' @param object a smooth object; see documentation for other methods in
#'   \code{mgcv}.
#' @param data a data.frame; see documentation for other methods in
#'   \code{mgcv}.
#' @export
Predict.matrix.gKRLS.smooth <- function(object, data) {
  
  X_test <- create_data_gKRLS(
    term_levels = object$term_levels, 
    term_class = object$term_class, 
    data_term = data[object$term], 
    terms = object$term, 
    allow.missing.levels = TRUE)

  if (!is.null(object$std_train)) {
    X_test <- sweep(X_test, 2, object$std_train$mean, FUN = "-")
    X_test <- X_test %*%
      object$std_train$whiten
  }

  X_test <- as.matrix(X_test)

  if (object$xt$return_raw) {
    return(X_test)
  }

  KS_test <- create_sketched_kernel(
    X_test = X_test, X_train = object$X_train,
    tS = t(object$sketch_matrix),
    bandwidth = object$bandwidth
  )
  
  KS_test <- sweep(KS_test, 2, object$KS_mean, FUN = "-")

  return(KS_test)
}
