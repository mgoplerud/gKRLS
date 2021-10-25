
#' Estimate marginal effects
#' @export
marginal_effect_base <- function(object, newdata, newkernel_X, method = 'cpp'){
  
  standardize <- object$internal$standardize
  
  f <- function(x,family){
    if (family == 'logit'){
      return(plogis(x))
    }else if (family == 'gaussian'){
      return(x)
    }else if (family == 'probit'){
      return(pnorm(x))
    }else if (family == 'poisson'){
      return(exp(x))
    }else{stop('Invalid family.')}
  }
  
  f_prime <- function(x,family){
    if (family == 'logit'){
      p <- plogis(x)
      return(p * (1-p))
    }else if (family == 'gaussian'){
      return(rep(1, length(x)))
    }else if (family == 'probit'){
      return(dnorm(x))
    }else if (family == 'poisson'){
      return(exp(x))
    }else{stop('Invalid family.')}
  }
  
  f_double_prime <- function(x,family){
    if (family == 'logit'){
      p <- plogis(x)
      return(p * (1-p) * (1-2 * p))
    }else if (family == 'gaussian'){
      return(rep(0, length(x)))
    }else if (family == 'probit'){
      return(dnorm(x) * -x)
    }else if (family == 'poisson'){
      return(exp(x))
    }else{stop('Invalid family.')}
  }
  
  
  prepped_data <- prepare_predict_data(object = object, newdata = newdata, 
                                       newkernel_X = newkernel_X,
                                       allow_missing_levels = allow_missing_levels)
  
  Z <- prepped_data$Z
  newdata_FE <- prepped_data$newdata_FE
  total_obs <- prepped_data$total_obs
  obs_in_both <- prepped_data$obs_in_both
  newdataKS <- prepped_data$newdataKS
  std_X_test <- prepped_data$std_newkernel_X
  std_X_train <- prepped_data$std_kernel_X
  pos_re <- prepped_data$n_pos
  rm(prepped_data); gc()
  
  family <- object$internal$family
  if (family$family == 'binomial' & family$link == 'logit'){
    family <- 'logit'
  }else if (family$family == 'binomial' & family$link == 'probit'){
    family <- 'probit'
  }else if (family$family == 'gaussian'){
    family <- 'gaussian'
  }else if (family$family == 'poisson'){
    family <- 'poisson'
  }else{stop('Invalid family!')}
  
  names_mfx <- c(colnames(newdata_FE), colnames(object$internal$kernel_X_train))
  
  W_Matrix <- object$internal$std_train$W_Matrix
  WX_train <- object$internal$kernel_X_train %*% W_Matrix
  WX_test <- newkernel_X %*% W_Matrix
  
  tS <- t(object$internal$sketch)
  bandwidth <- object$internal$bandwidth
  
  fe_mean <- object$fe$mean
  re_mean <- object$re$mean
  
  if (length(fe_mean) > 0){
    all_mean <- rbind(fe_mean, re_mean)
  }else{
    all_mean <- matrix(re_mean)
  }
  vcov_ridge <- object$vcov_ridge
  FE_matrix_test <- newdata_FE
  
  fmt_sd_y <- object$internal$sd_y
  
  pos_kern <- pos_re[which(names(pos_re) == '1 | kernel_RE')]
  stopifnot(names(pos_re)[length(pos_re)] == '1 | kernel_RE')
  
  #####################
  ### Begin Estimation
  # Uses the following arguments
  # std_X_test, std_X_train, 
  # WX_train, WX_test,
  # tS, bandwidth, FE_matrix_test,
  # fe_mean, re_mean, all_mean,
  # vcov_ridge
  #####################

  N_train <- nrow(std_X_train)
  N_test <- nrow(std_X_test)
  
  offset <- rep(0, N_test)

  orig_re <- re_mean
  offset_re <- orig_re[seq_len(length(orig_re) - pos_kern),,drop=F]
  # Are there any "pure" REs in the model? If not, simplify
  
  if (nrow(offset_re) == 0){
    re_mean <- orig_re
    any_Z <- FALSE
    Z <- sparseMatrix(i = 1, j = 1, x = 0)
  }else{
    any_Z <- TRUE
    re_mean <- orig_re[-seq_len(length(orig_re) - pos_kern), , drop = F]
    offset <- as.vector(Z[,seq_len(length(orig_re) - pos_kern), , drop = F] %*% offset_re)
    Z <- Z[,seq_len(length(orig_re) - pos_kern), , drop = F]
  }

  Sc <- t(tS) %*% re_mean
  ME_pointwise <- ME_pointwise_var <- matrix(NA, nrow(std_X_test), 
                                             ncol(FE_matrix_test) + ncol(W_Matrix))
  AME_grad <- matrix(0, nrow = nrow(all_mean), ncol = ncol(FE_matrix_test) + ncol(W_Matrix))
  SIZE_FE <- ncol(FE_matrix_test)
  SIZE_KERNEL <- ncol(W_Matrix)
  
  # Flag the columns that should be analyzed using first differences
  fd_flag <- which(apply(object$internal$kernel_X_train, MARGIN = 2, FUN=function(i){
    all(i %in% c(0,1))
  }))
  fd_flag <- names(fd_flag)
  fd_matrix <- matrix(data = 0, ncol = 2, nrow = SIZE_KERNEL + SIZE_FE)
  rownames(fd_matrix) <- names_mfx
  # The values to set the first difference to
  fd_matrix[fd_flag, 1] <- 0
  fd_matrix[fd_flag, 2] <- 1
  fd_flag <- (rownames(fd_matrix) %in% fd_flag)
  
  
  std_whiten <- object$internal$std_train$whiten
  std_mean <- object$internal$std_train$mean
  
  std_fd_matrix <- sweep(t(fd_matrix), 2, c(rep(0, SIZE_FE), std_mean), FUN = "-")
  std_fd_matrix <- std_fd_matrix %*% 
    bdiag(Diagonal(x = rep(0, SIZE_FE)), std_whiten)
  std_fd_matrix <- t(std_fd_matrix)
  std_fd_matrix <- as.matrix(std_fd_matrix)
  
  if (standardize == "Mahalanobis"){
    std_fd_matrix[,] <- 0
  }else{
    rownames(std_fd_matrix) <- names_mfx
  }
  
  type_mfx <- c(rep('deriv_FE', ncol(FE_matrix_test)), rep('deriv_Kern', ncol(std_X_test)))
  type_mfx[which(fd_flag)] <- 'FD'
  mahal <- standardize == 'Mahalanobis'
  raw_X_test <- newkernel_X
  std_whiten <- as.matrix(std_whiten)
  vcov_ridge <- as.matrix(vcov_ridge)
  W_Matrix <- as.matrix(W_Matrix)
  WX_test <- as.matrix(WX_test)
  WX_train <- as.matrix(WX_train)
  
  if (method == 'cpp'){
    out_ME <- cpp_gkrls_me(std_X_train = std_X_train, std_X_test = std_X_test,
                 bandwidth = bandwidth, family = family,
                 tZ = t(Z), offset = offset, any_Z = any_Z,
                 mahal = mahal, sd_y = fmt_sd_y, tS = tS, fe_mean = fe_mean, re_mean = as.vector(re_mean), 
                 SIZE_PARAMETER = nrow(all_mean), vcov_ridge = vcov_ridge, FE_matrix_test = FE_matrix_test, 
                 W_Matrix = W_Matrix, WX_test = WX_test, WX_train = WX_train, raw_X_test = raw_X_test, 
                 std_mean = std_mean, std_whiten = std_whiten, 
                 fd_matrix = fd_matrix, std_fd_matrix = std_fd_matrix, 
                 type_mfx = type_mfx)
  }else{
    
    for (i in seq_len(N_test)){
      
      std_X_i <- std_X_test[i,]
      WX_i <- WX_test[i,]
      if (any_Z){
        z_i <- Z[i,]
      }else{
        z_i <- NULL
      }
      offset_i <- offset[i]
      X_FE_i <- FE_matrix_test[i,]
      k_i <- sapply(seq_len(N_train), FUN=function(j){
        exp(-sum( (std_X_i - std_X_train[j,])^2)/bandwidth)
      })   
      tilde_k_i <- tS %*% k_i
      fe_i <- sum(fe_mean * X_FE_i)
      linpred_i <- fe_i + sum(tilde_k_i * re_mean) + offset_i
      f_prime_i <- f_prime(linpred_i, family)
      f_double_prime_i <- f_double_prime(linpred_i, family)
      
      # Get Marginal Effects for FE
      me_fe_i <- f_prime_i * fe_mean
      # Get Marginal Effects for Kernel Columns
      D_i <- lapply(seq_len(SIZE_KERNEL), FUN=function(p){WX_i[p] - WX_train[,p]})
      
      me_kern_i <- sapply(D_i, FUN=function(dip){
        as.numeric(f_prime_i * -2/bandwidth * (t(k_i) %*% Diagonal(x = dip) %*% Sc))
      })
      
      # Store all pointwise ME
      ME_pointwise[i,] <- c(me_fe_i, me_kern_i)
      
      # Get Gradient w.r.t. each MFX for the FE
      for (p in setdiff(seq_len(SIZE_FE), which(fd_flag))){
        
        grad_ME_FE_p_beta <- sapply(seq_len(SIZE_FE), FUN=function(p.prime){
          f_double_prime_i * fe_mean[p] * X_FE_i[p.prime] + 
            (p == p.prime) * f_prime_i
        })
        grad_ME_FE_p_c <- f_double_prime_i * fe_mean[p] * -2/bandwidth *
          (tS %*% k_i)
        if (any_Z){
          grad_ME_FE_p_alpha <- f_double_prime_i * fe_mean[p] * z_i
        }else{
          grad_ME_FE_p_alpha <- NULL
        }
        
        grad_ME_FE_p_c <- as.vector(grad_ME_FE_p_c)
        if (i == 1){stopifnot(names(pos_re)[length(pos_re)] == '1 | kernel_RE')}
        grad_ME_FE_p <- c(grad_ME_FE_p_beta, grad_ME_FE_p_alpha, grad_ME_FE_p_c)
        
        ME_pointwise_var[i,p] <- as.vector(t(grad_ME_FE_p) %*% vcov_ridge %*% grad_ME_FE_p)
        AME_grad[,p] <- AME_grad[,p] + grad_ME_FE_p
      }
      
      # Get the Gradient w.r.t. each MFX for the Kernel
      for (pkern in setdiff(seq_len(SIZE_KERNEL), -SIZE_FE + which(fd_flag))){
        
        meat_ip <- as.numeric(t(k_i) %*% Diagonal(x = D_i[[pkern]]) %*% Sc)
        grad_ME_K_p_c <- f_double_prime_i * (tS %*% k_i) * 
          (-2/bandwidth * meat_ip) +
          f_prime_i * -2/bandwidth * tS %*% Diagonal(x = D_i[[pkern]]) %*% k_i
        grad_ME_K_p_c <- as.vector(grad_ME_K_p_c)
        grad_ME_K_p_beta <- f_double_prime_i * (-2/bandwidth * meat_ip) * X_FE_i
        
        if (any_Z){
          grad_ME_K_p_alpha <- f_double_prime_i * -2/bandwidth * meat_ip * z_i
        }else{
          grad_ME_K_p_alpha <- NULL
        }
        
        if (i == 1){stopifnot(names(pos_re)[length(pos_re)] == '1 | kernel_RE')}
        
        grad_ME_K_p <- c(grad_ME_K_p_beta, grad_ME_K_p_alpha, grad_ME_K_p_c)
        
        ME_pointwise_var[i, SIZE_FE + pkern] <- as.vector(t(grad_ME_K_p) %*% vcov_ridge %*% grad_ME_K_p)
        AME_grad[ , SIZE_FE + pkern] <- AME_grad[ , SIZE_FE + pkern] + grad_ME_K_p
        
      }
      
      # Get the Gradient w.r.t. each MFX for the RE
      
      for (p_fd in which(fd_flag)){

        if (p_fd <= SIZE_FE){
          # If the flagged FD variable is in the fixed effect
          fe_i_no_p_fd <- fe_i - X_FE_i[p_fd] * fe_mean[p_fd]
          fe_i_0 <- fe_i_no_p_fd + fd_matrix[p_fd, 1] * fe_mean[p_fd]
          fe_i_1 <- fe_i_no_p_fd + fd_matrix[p_fd, 2] * fe_mean[p_fd]
          k_i_0 <- k_i
          k_i_1 <- k_i
        }else{
          # If the flagged FD variable is in the kernel  
          if (standardize != 'Mahalanobis'){
            
            std_X_train_p_fd <- std_X_train[,p_fd - SIZE_FE]
            std_X_train_i_p_df <- std_X_i[p_fd - SIZE_FE]
            
            # Remove the contribution of the FD column
            k_i_no_p_fd <-  k_i * exp( (std_X_train_i_p_df - std_X_train_p_fd)^2 / bandwidth)
            # Create the counterfactual kernels
            k_i_0 <- k_i_no_p_fd * exp( -(std_fd_matrix[p_fd,1] - std_X_train_p_fd)^2 / bandwidth)
            k_i_1 <- k_i_no_p_fd * exp( -(std_fd_matrix[p_fd,2] - std_X_train_p_fd)^2 / bandwidth)
          }else{
            
            raw_i <- newkernel_X[i,]
            if (fd_matrix[p_fd, 1] != raw_i[p_fd - SIZE_FE]){
              raw_i_0 <- raw_i
              raw_i_0[p_fd - SIZE_FE] <- fd_matrix[p_fd, 1]
              std_X_i_0 <- t(std_whiten) %*% (raw_i_0 - std_mean)
              k_i_0 <- sapply(seq_len(N_train), FUN=function(j){
                exp(-sum( (std_X_i_0 - std_X_train[j,])^2)/bandwidth)
              })   
            }else{
              k_i_0 <- k_i
            }
            
            if (fd_matrix[p_fd, 2] != raw_i[p_fd - SIZE_FE]){
              raw_i_1 <- raw_i
              raw_i_1[p_fd - SIZE_FE] <- fd_matrix[p_fd, 2]
              std_X_i_1 <- t(std_whiten) %*% (raw_i_1 - std_mean)
              k_i_1 <- sapply(seq_len(N_train), FUN=function(j){
                exp(-sum( (std_X_i_1 - std_X_train[j,])^2)/bandwidth)
              })   
            }else{
              k_i_1 <- k_i
            }
          }
          fe_i_0 <- fe_i
          fe_i_1 <- fe_i
        }
        linpred_i_1 <- fe_i_1 + as.numeric(t(k_i_1) %*% t(tS) %*% re_mean) + offset_i
        linpred_i_0 <- fe_i_0 + as.numeric(t(k_i_0) %*% t(tS) %*% re_mean) + offset_i
        
        ME_pointwise[i, p_fd] <- f(linpred_i_1, family) - f(linpred_i_0, family)
        # Get the gradient
        grad_ME_i_fd_beta <- (f_prime(linpred_i_1, family) - f_prime(linpred_i_0, family)) * X_FE_i
        grad_ME_i_fd_c <- tS %*% (f_prime(linpred_i_1, family) * k_i_1 - f_prime(linpred_i_0, family) * k_i_0)
        if (any_Z){
          grad_ME_i_fd_alpha <- (f_prime(linpred_i_1, family) - f_prime(linpred_i_0, family)) * z_i
        }else{
          grad_ME_i_fd_alpha <- NULL
        }
        if (i == 1){stopifnot(names(pos_re)[length(pos_re)] == '1 | kernel_RE')}
        
        
        grad_ME_FD <- c(grad_ME_i_fd_beta, grad_ME_i_fd_alpha, grad_ME_i_fd_c)
        ME_pointwise_var[i,p_fd] <- as.vector(t(grad_ME_FD) %*% vcov_ridge %*% grad_ME_FD)
        AME_grad[,p_fd] <- AME_grad[,p_fd] + grad_ME_FD
        
      }
      
    }
    
    AME_grad <- 1/N_test * AME_grad
    
    AME_pointwise_var <- apply(AME_grad, MARGIN = 2, 
                               FUN=function(i){as.numeric(t(i) %*% vcov_ridge %*% i)})
    AME_pointwise <- colMeans(ME_pointwise)
    
    ME_pointwise_var <- ME_pointwise_var * fmt_sd_y^2
    ME_pointwise <- ME_pointwise * fmt_sd_y
    
    AME_pointwise_var <- AME_pointwise_var * fmt_sd_y^2
    AME_pointwise <- AME_pointwise * fmt_sd_y
    
    AME_grad <- AME_grad * fmt_sd_y
    
    out_ME <- list(
      AME_grad = AME_grad,
      AME_pointwise = AME_pointwise,
      AME_pointwise_var = AME_pointwise_var,
      ME_pointwise = ME_pointwise,
      ME_pointwise_var = ME_pointwise_var
    )
    
  }
  
  names(out_ME$AME_pointwise) <- colnames(out_ME$ME_pointwise) <- colnames(out_ME$ME_pointwise_var) <- names_mfx
  colnames(out_ME$AME_grad) <- names(out_ME$AME_pointwise_var) <- names_mfx
  
  class(out_ME) <- c('gKRLS_ME')
  return(out_ME)
}