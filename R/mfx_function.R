
#' Marginal Effects by Derivative
#' 
#' Find the (average) marginal effect and standard error by taking the
#' analytical derivative of the conditional expectation function. This method is
#' designed to provide comparability with `KRLS` but is limited in the scenarios
#' where it can be applied. Specifically, it can only compute effects with a
#' single kernel and where variables are not both in the kernel and outside the
#' kernel.  
#' 
#' Because of these restrictions, users should rely on `calculate_effects` and
#' this function may be removed in future updates. If a derivative is really
#' desired, this can be approximated as described in the `calculate_effects`
#' documentation.
#' @export
legacy_marginal_effect <- function(object, newdata, keep = NULL,
  method = 'cpp'){
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

  full_lp <- predict(object, newdata = newdata, type = 'lpmatrix')
  model_offset <- attr(full_lp, 'model.offset')
  
  class_smooth <- lapply(object$smooth, FUN=function(i){class(i)})
  stopifnot(all(lengths(class_smooth) == 2))
  class_smooth <- sapply(class_smooth, FUN=function(i){i[1]})
  
  if (sum(class_smooth == 'gKRLS.smooth') != 1){
    stop('marginal_effects can only be run with one kernel...for now...')
  }
  kern_smooth <- object$smooth[[which(class_smooth == 'gKRLS.smooth')]]
  
  # Get the fixed effects
  fe_id <- setdiff(seq_len(ncol(full_lp)),
        unlist(lapply(object$smooth, FUN=function(i){i$first.par:i$last.par}))
  )
  if (length(fe_id) > 0){
    FE_matrix_test <- full_lp[, fe_id, drop = F]
  }else{
    FE_matrix_test <- matrix(nrow = nrow(full_lp), ncol = 0)
  }
  colnames(FE_matrix_test) <- colnames(full_lp)[fe_id]
  
  kernel_id <- unlist(kern_smooth[c('first.para', 'last.para')])
  kernel_id <- kernel_id[1]:kernel_id[2]
  # All *other* smoothes, e.g. random effects
  Z <- full_lp[, -c(fe_id, kernel_id), drop = F]
  
  # Get the raw design and UNDO standardization
  std_X_train <- kern_smooth$X_train
  raw_X_train <- sweep(std_X_train %*% solve(kern_smooth$std_train$whiten), 
    2, kern_smooth$std_train$mean, FUN = "+")

  standardize <- kern_smooth$xt$standardize  
  
  std_whiten <- kern_smooth$std_train$whiten
  std_mean <- kern_smooth$std_train$mean
  W_Matrix <- kern_smooth$std_train$W_Matrix
  
  kern_smooth$xt$return_raw <- TRUE
  kern_smooth$std_train <- NULL
  
  raw_X_test <- PredictMat(kern_smooth, data = newdata)
  family <- object$family
  
  kern_names <- colnames(std_X_train)
  fe_names <- colnames(FE_matrix_test)
  names_mfx <- c(fe_names, kern_names)
  if (!identical(kern_smooth$term, colnames(std_X_train))){
    message('Names of kern and std_X_train seem misaligned...')
    print(kern_smooth$term)
    print(colnames(std_X_train))
    stop()
  }
  
  WX_train <- raw_X_train %*% W_Matrix
  WX_test <- raw_X_test %*% W_Matrix
  
  std_X_test <- sweep(raw_X_test, 2, std_mean, FUN = "-")
  std_X_test <- as.matrix(std_X_test %*% std_whiten)
  
  tS <- t(kern_smooth$sketch_matrix)
  bandwidth <- kern_smooth$bandwidth
  
  all_mean <- matrix(coef(object))
  if (length(fe_id) > 0){
    fe_mean <- matrix(coef(object)[fe_id])
    re_mean <- all_mean[-fe_id,,drop=F]
  }else{
    re_mean <- all_mean
    fe_mean <- numeric(0)
  }
  vcov_ridge <- vcov(object)
  fmt_sd_y <- 1
  
  N_train <- nrow(std_X_train)
  N_test <- nrow(std_X_test)

  if (length(model_offset) == 1){
    offset <- rep(0, N_test)
  }else{
    if (length(model_offset) != N_test){stop('Model offest is incorrect length.')}
    offset <- model_offset
  }
  
  if (ncol(Z) == 0){
    any_Z <- FALSE
    Z <- sparseMatrix(i = 1, j = 1, x = 0)
  }else{
    any_Z <- TRUE
    offset_mean <- matrix(coef(object)[-c(fe_id, kernel_id)])
    re_mean <- matrix(coef(object)[kernel_id])
    offset <- as.vector(offset + Z %*% offset_mean)
    Z <- drop0(Z)
  }
  Sc <- t(tS) %*% re_mean
  
  fd_flag <- kern_smooth$fd_flag
  
  if (family$family == 'binomial' & family$link == 'logit'){
    family <- 'logit'
  }else if (family$family == 'binomial' & family$link == 'probit'){
    family <- 'probit'
  }else if (family$family == 'gaussian'){
    family <- 'gaussian'
  }else if (family$family == 'poisson'){
    family <- 'poisson'
  }else{stop('Invalid family!')}
  
  
  ME_pointwise <- ME_pointwise_var <- matrix(NA, nrow(std_X_test), 
      ncol(FE_matrix_test) + ncol(W_Matrix))
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
  
  if (standardize == "Mahalanobis"){
    std_fd_matrix[,] <- 0
  }else{
    rownames(std_fd_matrix) <- names_mfx
  }
  
  type_mfx <- c(rep('deriv_FE', ncol(FE_matrix_test)), rep('deriv_Kern', ncol(std_X_test)))
  type_mfx[which(fd_flag)] <- 'FD'
  mahal <- standardize == 'Mahalanobis'
  std_whiten <- as.matrix(std_whiten)
  vcov_ridge <- as.matrix(vcov_ridge)
  W_Matrix <- as.matrix(W_Matrix)
  WX_test <- as.matrix(WX_test)
  WX_train <- as.matrix(WX_train)
  
  if (!is.null(keep)){
    # Variables to calculate AME for
    kern_to_fit <- intersect(kern_names, keep)
    fe_to_fit <- intersect(fe_names, keep)
    missing_provided <- setdiff(setdiff(keep, kern_to_fit), fe_to_fit)
    all_to_fit <- c(fe_to_fit, kern_to_fit)
    if (length(missing_provided) != 0){
      warning('Some variables in "keep" not located in data.')
    }
    fit_position <- which(rownames(fd_matrix) %in% all_to_fit)
    mfx_counter <- seq_len(length(fit_position))
  }else{
    fit_position <- seq_len(nrow(fd_matrix))
    mfx_counter <- seq_len(nrow(fd_matrix))
  }

  if (method == 'cpp'){
    
    cpp_fit_position <- as.integer(fit_position - 1)
    cpp_mfx_counter <- as.integer(mfx_counter - 1)
    
    out_ME <- cpp_gkrls_me(std_X_train = std_X_train, std_X_test = std_X_test,
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
                mfx_counter = cpp_mfx_counter)
    }else{

    AME_grad <- AME_grad[,seq_len(length(mfx_counter)), drop = F]
    ME_pointwise <- ME_pointwise[,seq_len(length(mfx_counter)), drop = F]
    ME_pointwise_var <- ME_pointwise_var[,seq_len(length(mfx_counter)), drop = F]
    
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

      for (raw_pos in mfx_counter){
        p <- fit_position[raw_pos]
        if (type_mfx[p] == 'deriv_FE'){
          
          # Get Marginal Effects for FE
          me_fe_ip <- f_prime_i * fe_mean[p]
          
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
          
          ME_pointwise[i, raw_pos] <- me_fe_ip
          ME_pointwise_var[i,raw_pos] <- as.vector(t(grad_ME_FE_p) %*% vcov_ridge %*% grad_ME_FE_p)
          AME_grad[,raw_pos] <- AME_grad[,raw_pos] + grad_ME_FE_p
          
        }else if (type_mfx[p] == 'deriv_Kern'){
          pkern <- p - SIZE_FE
          
          D_ip <- WX_i[pkern] - WX_train[,pkern]
  
          meat_ip <- as.numeric(t(k_i) %*% Diagonal(x = D_ip) %*% Sc)
          me_kern_ip <- meat_ip * f_prime_i * -2/bandwidth
          
          grad_ME_K_p_c <- f_double_prime_i * (tS %*% k_i) * 
            (-2/bandwidth * meat_ip) +
            f_prime_i * -2/bandwidth * tS %*% Diagonal(x = D_ip) %*% k_i
          grad_ME_K_p_c <- as.vector(grad_ME_K_p_c)
          grad_ME_K_p_beta <- f_double_prime_i * (-2/bandwidth * meat_ip) * X_FE_i
          
          if (any_Z){
            grad_ME_K_p_alpha <- f_double_prime_i * -2/bandwidth * meat_ip * z_i
          }else{
            grad_ME_K_p_alpha <- NULL
          }
          
          if (i == 1){stopifnot(names(pos_re)[length(pos_re)] == '1 | kernel_RE')}
          
          grad_ME_K_p <- c(grad_ME_K_p_beta, grad_ME_K_p_alpha, grad_ME_K_p_c)
          
          ME_pointwise[i, raw_pos] <- me_kern_ip
          ME_pointwise_var[i, raw_pos] <- as.vector(t(grad_ME_K_p) %*% vcov_ridge %*% grad_ME_K_p)
          AME_grad[ , raw_pos] <- AME_grad[ , raw_pos] + grad_ME_K_p
          
        }else if (type_mfx[p] == 'FD'){
          
          p_fd <- p - SIZE_FE
          
          if (p <= SIZE_FE){
            # If the flagged FD variable is in the fixed effect
            fe_i_no_p_fd <- fe_i - X_FE_i[p] * fe_mean[p]
            fe_i_0 <- fe_i_no_p_fd + fd_matrix[p, 1] * fe_mean[p]
            fe_i_1 <- fe_i_no_p_fd + fd_matrix[p, 2] * fe_mean[p]
            k_i_0 <- k_i
            k_i_1 <- k_i
          }else{
            # If the flagged FD variable is in the kernel  
            if (standardize != 'Mahalanobis'){
              
              std_X_train_p_fd <- std_X_train[, p_fd]
              std_X_train_i_p_df <- std_X_i[p_fd]
              
              # Remove the contribution of the FD column
              k_i_no_p_fd <-  k_i * exp( (std_X_train_i_p_df - std_X_train_p_fd)^2 / bandwidth)
              # Create the counterfactual kernels
              k_i_0 <- k_i_no_p_fd * exp( -(std_fd_matrix[p,1] - std_X_train_p_fd)^2 / bandwidth)
              k_i_1 <- k_i_no_p_fd * exp( -(std_fd_matrix[p,2] - std_X_train_p_fd)^2 / bandwidth)
            }else{
              
              raw_i <- newkernel_X[i,]
              
              if (fd_matrix[p, 1] != raw_i[p_fd]){
                raw_i_0 <- raw_i
                raw_i_0[p_fd] <- fd_matrix[p, 1]
                std_X_i_0 <- t(std_whiten) %*% (raw_i_0 - std_mean)
                k_i_0 <- sapply(seq_len(N_train), FUN=function(j){
                  exp(-sum( (std_X_i_0 - std_X_train[j,])^2)/bandwidth)
                })   
              }else{
                k_i_0 <- k_i
              }
              
              if (fd_matrix[p, 2] != raw_i[p_fd]){
                raw_i_1 <- raw_i
                raw_i_1[p_fd] <- fd_matrix[p, 2]
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
        
        ME_pointwise[i, raw_pos] <- f(linpred_i_1, family) - f(linpred_i_0, family)
        ME_pointwise_var[i, raw_pos] <- as.vector(t(grad_ME_FD) %*% vcov_ridge %*% grad_ME_FD)
        AME_grad[ , raw_pos] <- AME_grad[ , raw_pos] + grad_ME_FD
        
        }else{stop('type_mfx INVALID')}
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
  
  names(out_ME$AME_pointwise) <- colnames(out_ME$ME_pointwise) <- colnames(out_ME$ME_pointwise_var) <- names_mfx[fit_position]
  colnames(out_ME$AME_grad) <- names(out_ME$AME_pointwise_var) <- names_mfx[fit_position]
  
  out_ME$type <- type_mfx[fit_position[mfx_counter]]
  out_ME$mfx_type <- 'legacy'
  
  class(out_ME) <- c('gKRLS_mfx')
  
  return(out_ME)
}



