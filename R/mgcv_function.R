#' @import mgcv
#' @import RcppEigen
#' @importFrom Rcpp sourceCpp
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

  if (length(object$term) > 1) {
    length_data <- sapply(data, length)
    if (length(unique(length_data)) != 1) {
      stop("Error: All input variables should have same length. The syntax is s(a,b,c,d)")
    }
    X <- do.call("cbind", data[object$term])
  } else {
    X <- matrix(data[[object$term]])
  }

  if (NCOL(X) == 1 & is.null(ncol(X))) {
    X <- matrix(X)
  }

  fd_flag <- which(apply(X, MARGIN = 2, FUN = function(i) {
    all(i %in% c(0, 1))
  }))

  std_X <- standardize_design(kernel_X = X, standardize = object$xt$standardize)
  X <- std_X$std_kernel_X
  std_train <- std_X$std_train
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
  if (object$xt$force_base) {
    KS <- exp(-base_kernel(X, X_train) / bandwidth) %*% sketch_matrix
  } else {
    KS <- create_sketched_kernel(
      X_test = X,
      X_train = X_train, tS = t(sketch_matrix), bandwidth = bandwidth
    )
  }

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

  # Required elements for "gam"
  object$rank <- ncol(KS)
  object$null.space.dim <- 0
  object$df <- ncol(KS)
  if (object$xt$no.rescale) {
    object$no.rescale <- TRUE
  }
  object$te.ok <- 0
  object$plot.me <- FALSE
  object$C <- matrix(nrow = 0, ncol = ncol(KS))

  class(object) <- "gKRLS.smooth" # Give object a class
  object
}

#' @export
Predict.matrix.gKRLS.smooth <- function(object, data) {
  if (length(object$term) > 1) {
    X_test <- do.call("cbind", data[object$term])
  } else {
    X_test <- matrix(data[[object$term]])
  }

  if (!is.null(object$std_train)) {
    X_test <- sweep(X_test, 2, object$std_train$mean, FUN = "-")
    X_test <- X_test %*%
      object$std_train$whiten
  }

  X_test <- as.matrix(X_test)

  if (object$xt$return_raw) {
    return(X_test)
  }

  if (object$xt$force_base) {
    KS_test <- exp(-base_kernel(X_test, object$X_train) / object$bandwidth) %*% object$sketch_matrix
  } else {
    KS_test <- create_sketched_kernel(
      X_test = X_test, X_train = object$X_train,
      tS = t(object$sketch_matrix),
      bandwidth = object$bandwidth
    )
  }
  KS_test <- sweep(KS_test, 2, object$KS_mean, FUN = "-")

  return(KS_test)
}
