#' @importFrom stats rnorm rbinom
create_sketch_matrix <- function(N, sketch_size, sketch_prob = NULL, sketch_method) {

  if (sketch_method == "gaussian") {
    
    S <- t(matrix(rnorm(N * sketch_size), nrow = N))
    S <- S * sqrt(1 / sqrt(sketch_size))
    
  } else if (sketch_method == "bernoulli") {
    
    if (is.null(sketch_prob)) {
      stop('sketch method "bernoulli" requires a probability.')
    }

    S <- matrix(rbinom(N * sketch_size, 1, prob = sketch_prob), nrow = N)
    S <- S * 1 / sqrt(sketch_prob * (1 - sketch_prob))
    S <- S * 1 / sqrt(sketch_size)
    S <- t(S)
    
  } else if (sketch_method == "none") {
    
    S <- as.matrix(Diagonal(n = N))
    
  }

  S <- as.matrix(S)
  
  return(S)
}

#' @importFrom stats cov var
standardize_design <- function(kernel_X, standardize) {
  names_kx <- colnames(kernel_X)

  zero_columns <- NULL
  
  if (standardize == "Mahalanobis") {
    std_mean_X <- colMeans(kernel_X)
    # var(X) = A -> var(X A) = A^T var(X) A
    # A Q L Q^T A -> A = Q sqrt(L)^{-1}

    eigen_kernel_X <- eigen(cov(kernel_X))
    eigen_kernel_X$values <- ifelse(
      eigen_kernel_X$values < sqrt(.Machine$double.eps), 0,
      eigen_kernel_X$values
    )

    std_whiten_X <- with(
      eigen_kernel_X,
      vectors %*% Diagonal(x = ifelse(values == 0, 0, 1 / sqrt(values)))
    )

    W_Matrix <- with(
      eigen_kernel_X,
      vectors %*% Diagonal(x = ifelse(values == 0, 0, 1 / values)) %*% t(vectors)
    )

    zero_columns <- apply(std_whiten_X, MARGIN = 2, FUN = function(i) {
      all(i == 0)
    })
    zero_columns <- which(zero_columns == TRUE)
    if (length(zero_columns) > 0) {
      message("Original X matrix is not full rank. Using an effective rank of ", ncol(kernel_X) - length(zero_columns), " for bandwidth throughout.")
      std_whiten_X <- std_whiten_X[, -zero_columns, drop = F]
    }
    colnames(std_whiten_X) <- names_kx <- paste0("m_", 1:ncol(std_whiten_X))
  } else if (standardize == "scaled") {
    std_mean_X <- colMeans(kernel_X)
    vr_X <- apply(kernel_X, MARGIN = 2, var)
    std_whiten_X <- Diagonal(x = 1 / sqrt(ifelse(vr_X == 0, 1, vr_X)))
    W_Matrix <- Diagonal(x = ifelse(vr_X == 0, 0, 1 / vr_X))
  } else if (standardize == "none") {
    std_mean_X <- rep(0, ncol(kernel_X))
    std_whiten_X <- Diagonal(n = ncol(kernel_X))
    W_Matrix <- Diagonal(n = ncol(kernel_X))
  } else {
    stop("Invalid standardization method.")
  }
  # demean X
  std_kernel_X <- sweep(kernel_X, 2, std_mean_X, FUN = "-")
  # standardize
  std_kernel_X <- std_kernel_X %*% std_whiten_X
  std_kernel_X <- as.matrix(std_kernel_X)

  std_train <- list(
    mean = std_mean_X,
    whiten = std_whiten_X,
    W_Matrix = W_Matrix
  )

  colnames(std_kernel_X) <- names_kx

  return(list(std_kernel_X = std_kernel_X, 
              std_train = std_train,
              zero_columns = zero_columns))
}

# Calculate kernel in base R
base_kernel <- function(X, Y) {
  edist <- function(x, y) {
    sum((x - y)^2)
  }
  out <- outer(1:nrow(X), 1:nrow(Y), FUN = Vectorize(function(x, y) {
    edist(X[x, ], Y[y, ])
  }))
  return(out)
}

# Calibrate b modifying the procedure in Hartman, Hazlett, and Sterbenz (2021)
# for sketched kernels; see "gKRLS_addendum.pdf" on the GitHub repo for more
# information.
#' @importFrom stats optimize
calibrate_bandwidth <- function(X, S = NULL, id_S = NULL, tol = sqrt(.Machine$double.eps)){
  
  # subsampling sketch as only "id_S" is provided
  if (!is.null(id_S)){ 
    design_X <- X[id_S,]
    raw_KSt <- create_sketched_kernel(
      X_test = X, X_train = design_X,
      S = diag(length(id_S)),
      bandwidth = 1, raw = T)
    type <- 'subsampling'
  }else if (is.null(id_S) & !is.null(S)){
    type <- 'direct'
  }else{stop('both id_S and S cannot be NULL')}
  
  f <- function(ln_b, raw_KSt, type){
    b <- exp(ln_b)
    if (type == 'subsampling'){
      KSt <- exp(-raw_KSt/b)
      reduced_K <- KSt[id_S,]
    }else{
      KSt <- create_sketched_kernel(
        X_test = X, X_train = X,
        S = S, bandwidth = b)
      reduced_K <- S %*% KSt
    }
    eigen_K <- eigen(reduced_K)
    if (any(eigen_K$values < tol)){
      nonzero_ev <- which(eigen_K$values >= tol)
      eigen_K$vectors <- eigen_K$vectors[,nonzero_ev]
      eigen_K$values <- eigen_K$values[nonzero_ev]
    }
    A <- eigen_K$vectors %*% Diagonal(x=1/sqrt(eigen_K$values))
    KSt_A <- KSt %*% A
    AtA <- t(KSt_A) %*% KSt_A
    
    N <- nrow(KSt_A)^2 
    sum_K <- sum( colSums( (KSt_A) )^2 )
    sum_Ksq <-  sum((KSt_A %*% AtA) * KSt_A)
    # sum_Ksq <- sum(
    #   apply(KSt_A, MARGIN = 1, FUN=function(i){sum( (KSt_A %*% i)^2 )})
    # )
    var_K <- sum_Ksq/N - (sum_K/N)^2
    var_K <- var_K * N/(N-1)
    return(var_K)
  }
  
  opt <- optimize(f = f, raw_KSt = raw_KSt, type = type,
                  interval = c(-10, 10), maximum = TRUE)
  opt$maximum <- exp(opt$maximum)
  dist_from_boundary <- min(abs(opt$maximum - exp(c(-10, 10))))
  if (dist_from_boundary < 1e-4){warning('kernel bandwidth optimized near boundary; try rescaling covariates')}
  return(opt$maximum)
}