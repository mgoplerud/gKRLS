
#' @importFrom Matrix sparseMatrix
create_data_gKRLS <- function(term_levels, term_class, data_term, terms, allow.missing.levels = FALSE){

  out <- mapply(term_levels, term_class, data_term, terms, SIMPLIFY = FALSE, 
    FUN=function(l_i, c_i, d_i, n_i){
    if (is.null(l_i)){
      if (!(any(c_i %in% c('numeric', 'integer', 'double', 'float')))){
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
#' @importFrom stats runif
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

  bandwidth <- object$xt$bandwidth
  if (is.null(bandwidth)) {
    bandwidth <- ncol(X)
  }
  
  # Create the Kernel
  N <- nrow(X)

  if (object$xt$sketch_method != 'custom'){
    sketch_size_raw <- object$xt$sketch_size_raw
    sketch_multiplier <- object$xt$sketch_multiplier
    # Either multiply N^(1/3) by multiplier *or* use raw value
    if (!is.null(sketch_multiplier)) {
      sketch_size <- floor(ceiling(N^(1 / 3)) * sketch_multiplier)
    } else {
      sketch_size <- sketch_size_raw
    }
  }else{
    sketch_size <- sketch_size_raw <- NA
  }


  X_train <- X

  subsampling_id <- NULL
  if (is.na(sketch_size) & object$xt$sketch_method != 'custom') {
    sketch_matrix <- diag(N)
  } else {
    
    if (object$xt$sketch_method != 'custom'){
      if (sketch_size >= N) {
        warning("Sketch size exceeds size of data")
      }
      if (sketch_size < 0) { # Do N - s
        sketch_size <- N + sketch_size
      }
    }
    # Remove leverage code for now...    
    # else if (object$xt$sketch_method == 'subsampling_leverage'){
    #   
    #   if (sketch_size > N) {
    #     stop("Subsampling requires sketch_size < N.")
    #   }
    #   leverage_dim <- min(c(nrow(X), 10 * sketch_size))
    #   message(paste0('Computing leverage scores with rank ', leverage_dim))
    #   leverage_scores <- function(X, bandwidth, k){stop('SET UP LEVERAGE SCORES')}
    #   lscores <- leverage_scores(X = X, bandwidth = bandwidth, k = leverage_dim)
    #   subsampling_id <- which(sketch_size * lscores >= runif(nrow(X)))
    #   message(paste0(length(subsampling_id), ' sampled using leverage scores: ', sketch_size, ' was requested.'))
    #   X_train <- X[subsampling_id, , drop = F]
    #   sketch_matrix <- diag(length(subsampling_id)) * sqrt(N/sketch_size)
    #   
    # }
    if (object$xt$sketch_method == 'custom'){
      subsampling_id <- object$xt$sketch_vector
      X_train <- X[subsampling_id, , drop = F]
      sketch_size <- length(subsampling_id)
      sketch_matrix <- diag(length(subsampling_id)) * sqrt(N/sketch_size)
    }else if (object$xt$sketch_method == "subsampling") {
      if (sketch_size > N) {
        stop("Subsampling requires sketch_size < N.")
      }
      subsampling_id <- sample(1:N, sketch_size)
      X_train <- X[subsampling_id, , drop = F]
      sketch_matrix <- diag(length(subsampling_id)) * sqrt(N/sketch_size)
      
    } else {
      sketch_matrix <- create_sketch_matrix(N, sketch_size, object$xt$sketch_prob, object$xt$sketch_method)
    }
    
  }

  KSt <- create_sketched_kernel(
    X_test = X,
    X_train = X_train, S = sketch_matrix, bandwidth = bandwidth
  )
  
  ev_orig <- NULL
  P_orig <- NULL
  if (!object$fixed) {
    # S K S^T is the Penalty Term
    if (!is.na(sketch_size) & object$xt$sketch_method %in% c("custom", "subsampling", "subsampling_leverage") ) {
      Penalty <- t(KSt[subsampling_id, ]) %*% t(sketch_matrix)
    } else {
      Penalty <- t(KSt) %*% t(sketch_matrix)
    }
    P_orig <- Penalty
    # Penalty <- Penalty/sum(Penalty^2)

    if (object$xt$remove_instability) {
      old_size_S <- ncol(KSt)

      truncate_eigen_penalty <- object$xt$truncate.eigen.tol

      eigen_sketch <- eigen(Penalty, symmetric = TRUE)
      ev_orig <- eigen_sketch$values
      eigen_sketch$values <- ifelse(eigen_sketch$values < truncate_eigen_penalty, 0, eigen_sketch$values)
      nonzero <- which(eigen_sketch$values != 0)
    
      if (length(nonzero) == 0) {
        stop('After truncation, all eigenvectors are removed. Decrease "truncate.eigen.tol" or set "remove_instability" = FALSE to proceed.')
      }
      KSt <- as.matrix(KSt %*% (eigen_sketch$vectors[, nonzero] %*%
        Diagonal(x = 1 / sqrt(eigen_sketch$values[nonzero]))))
      # P = S K S^T = Q Lambda Q^T
      # S = t(t(S) %*% Q %*% Lambda^{-1/2})
      # S = Lambda^{-1/2} %*% Q^T %*% S
      # as note that KS^T = KS^T Q Lambda^{-1/2} after transformation.
      sketch_matrix <- t( t(sketch_matrix) %*% as.matrix(eigen_sketch$vectors[, nonzero] %*%
        Diagonal(x = 1 / sqrt(eigen_sketch$values[nonzero]))) )
      Penalty <- diag(ncol(KSt))
      if (length(nonzero) != length(eigen_sketch$values)) {
        message(paste("Sketch dimension decreased from", old_size_S, "to", ncol(KSt), "because of singular penalty matrix."))
      }
    }

    object$S <- list()
    object$S[[1]] <- Penalty
  } else {
    object$S <- NULL
    warning("fx=TRUE may result in severe overfitting. Use with caution.")
  }

  if (object$xt$demean_kernel) {
    KSt_mean <- colMeans(KSt)
    KSt <- apply(KSt, MARGIN = 2, FUN = function(i) {
      i - mean(i)
    })
  } else {
    KSt_mean <- rep(0, ncol(KSt))
  }
  # Required elements for kernel prediction
  object$ev_orig <- P_orig
  object$KSt_mean <- KSt_mean
  object$X <- KSt
  object$copy <- KSt
  object$X_train <- X_train
  object$bandwidth <- bandwidth
  object$sketch_matrix <- sketch_matrix
  object$std_train <- std_train
  object$fd_flag <- fd_flag
  object$term_levels <- term_levels
  object$term_class <- term_class
  object$subsampling_id <- subsampling_id
  # Required elements for "gam"
  object$rank <- ncol(KSt)
  object$null.space.dim <- 0
  object$df <- ncol(KSt)
  # If rescale_penalty is NOT true, then set "no.rescale" to TRUE, i.e. to not rescale.
  # otherwise, leave null.
  if (!object$xt$rescale_penalty) {
    object$no.rescale <- TRUE
  }
  object$te.ok <- 0
  object$plot.me <- FALSE
  object$C <- matrix(nrow = 0, ncol = ncol(KSt))

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

  KSt_test <- create_sketched_kernel(
    X_test = X_test, X_train = object$X_train,
    S = object$sketch_matrix,
    bandwidth = object$bandwidth
  )
  
  if (object$xt$demean_kernel) {
    KSt_test <- sweep(KSt_test, 2, object$KSt_mean, FUN = "-")
  }
  
  return(KSt_test)
}

#' @importFrom stats poly
#' @export
smooth.construct.unregpoly.smooth.spec <- function(object, data, knots) {
  if (length(knots) != 0){stop('"knots" not used for unregularized-polynomial.')}
  if (is.null(object$xt)) {
    object$xt <- NULL
  }
  X <- do.call('poly', c(list(x = data[[object$term]]), object$xt)) 
  object$orig_poly <- X
  object$X <- as.matrix(X)
  object$rank <- NA
  object$null.space.dim <- 0
  object$df <- ncol(X)
  object$no.rescale <- TRUE
  object$te.ok <- 0
  object$plot.me <- FALSE
  object$C <- matrix(nrow = 0, ncol = ncol(X))
  object$S <- list()
  object$fixed <- TRUE
  class(object) <- 'unregpoly.smooth'
  return(object)
}

#' @export
Predict.matrix.unregpoly.smooth <- function(object, data) {
  X_new <- predict(object = object$orig_poly, newdata = data[[object$term]])
  return(X_new)
}

# leverage_scores <- function(X, bandwidth, k){
#   
#   kern_prod <- function(x, args){
#     create_sketched_kernel(X_test = args$X, X_train = args$X, S = x, bandwidth = args$bandwidth)
#   }
#   
#   eig_func <- RSpectra::eigs_sym(
#     kern_prod, k = k, 
#     which = 'LM', n = nrow(X),
#     args = list(X = X, bandwidth = bandwidth))
#   
#   leverage <- rowMeans(eig_func$vectors^2)
#   return(leverage)
# }
