create_sketch_matrix <- function(N, sketch_size, sketch_prob = NULL, sketch_method){
  if (sketch_method == 'gaussian'){
    
    S <- matrix(rnorm(N * sketch_size), nrow = N)
    S <- S * 1/sqrt(sketch_size)
    
  }else if (sketch_method == 'bernoulli'){
    if (is.null(sketch_prob)){stop('sketch method "bernoulli" requires a probability.')}
    
    S <- matrix(rbinom(N * sketch_size, 1, prob = sketch_prob), nrow = N)
    S <- S * 1/sqrt(sketch_prob * (1 - sketch_prob))
    S <- S * 1/sqrt(sketch_size)
    
  }else if (sketch_method == 'none'){
    S <- as.matrix(Diagonal(n = N))
  }
  
  S <- as.matrix(S)  
}

standardize_design <- function(kernel_X, standardize){

  names_kx <- colnames(kernel_X)
  
  if (standardize == 'Mahalanobis'){
    
    std_mean_X <- colMeans(kernel_X)
    # var(X) = A -> var(X A) = A^T var(X) A 
    # A Q L Q^T A -> A = Q sqrt(L)^{-1}
    
    eigen_kernel_X <- eigen(cov(kernel_X))
    eigen_kernel_X$values <- ifelse(
      eigen_kernel_X$values < sqrt(.Machine$double.eps), 0, 
      eigen_kernel_X$values
    )
    
    std_whiten_X <- with(eigen_kernel_X, 
                         vectors %*% Diagonal(x = ifelse(values == 0, 0, 1/sqrt(values))))
    
    W_Matrix <- with(eigen_kernel_X, 
                     vectors %*% Diagonal(x = ifelse(values == 0, 0, 1/values)) %*% t(vectors)
    )
    
    zero_columns <- apply(std_whiten_X, MARGIN = 2, FUN=function(i){all(i == 0)})
    zero_columns <- which(zero_columns == TRUE)  
    if (length(zero_columns) > 0){
      print(paste('Original X matrix is not full rank. Using an effective rank of', ncol(kernel_X) - length(zero_columns), 'throughout.'))
      std_whiten_X <- std_whiten_X[,-zero_columns,drop=F]
    }
    colnames(std_whiten_X) <- names_kx <- paste0('m_', 1:ncol(std_whiten_X))
  }else if (standardize == 'scaled'){
    
    std_mean_X <- colMeans(kernel_X)
    vr_X <- apply(kernel_X, MARGIN = 2, var)
    std_whiten_X <- Diagonal(x = 1/sqrt(ifelse(vr_X == 0, 1, vr_X)))
    W_Matrix <- Diagonal(x = ifelse(vr_X == 0, 0, 1/vr_X))
    
  }else if (standardize == 'none'){
    
    std_mean_X <- rep(0, ncol(kernel_X))
    std_whiten_X <- Diagonal(n = ncol(kernel_X))
    W_Matrix <- Diagonal(n = ncol(kernel_X))
    
  }else{stop("Invalid standardization method.")}
  # demean X
  std_kernel_X <- sweep(kernel_X, 2, std_mean_X, FUN = "-")
  # standardize
  std_kernel_X <- std_kernel_X %*% std_whiten_X
  std_kernel_X <- as.matrix(std_kernel_X)  
  
  std_train <- list(
    mean = std_mean_X,
    whiten = std_whiten_X,
    W_Matrix = W_Matrix)
  
  colnames(std_kernel_X) <- names_kx
  
  return(list(std_kernel_X = std_kernel_X, std_train = std_train))
}