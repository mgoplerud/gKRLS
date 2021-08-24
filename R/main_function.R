
#' Fit a generalized KRLS model using Laplace approximation
#' @useDynLib gKRLS
#' @import Matrix
#' @export
gKRLS <- function(y, 
      formula, data,
      kernel_X, family, 
      sketch_size = ceiling(nrow(kernel_X)^(1/3)) * 5, 
      sketch_method = 'gaussian', sketch_prob = NULL, bandwidth = NULL,
      standardize = 'scaled',
      control = list(init_var =  1, intercept = TRUE,
                     truncate_eigen = sqrt(.Machine$double.eps))){
  
  if (!inherits(family, 'family')){
    stop('Provide "family" as in glm, glmer, etc.')
  }
  
  # Parse data and get response
  data <- model.frame(subbars(formula), data)
  nobs_complete <- nrow(data)
  y <- model.response(data)
  rownames(y) <- NULL

  if (family$family == 'gaussian'){
    sd_y <- sd(y)
    mean_y <- mean(y)
    fmt_y <- (y - mean(y))/sd(y)
    data[[1]] <- fmt_y
    rm(fmt_y)
  }else{
    sd_y <- 1
    mean_y <- 0
    
  }
  
  if (standardize == 'Mahalanobis'){
    
    std_mean_X <- colMeans(kernel_X)
    std_var_X <- rep(1, ncol(kernel_X))
    # var(X) = A -> var(X A) = A^T var(X) A 
    # A Q L Q^T A -> A = Q sqrt(L)^{-1}
    std_whiten_X <- with(eigen(cov(kernel_X)), 
          vectors %*% Diagonal(x = ifelse(values < sqrt(.Machine$double.eps), 0, 1/sqrt(values))))
    zero_columns <- apply(std_whiten_X, MARGIN = 2, FUN=function(i){all(i == 0)})
    zero_columns <- which(zero_columns == TRUE)  
    if (length(zero_columns) > 0){
      print(paste('Original X matrix is not full rank. Using an effective rank of', ncol(X) - length(zero_columns), 'throughout.'))
      std_whiten_X <- std_whiten_X[,-zero_columns,drop=F]
    }
  }else if (standardize == 'scaled'){
    
    std_mean_X <- colMeans(kernel_X)
    std_var_X <- apply(kernel_X, MARGIN = 2, var)
    std_whiten_X <- Diagonal(n = ncol(kernel_X)) 
    
    # Set "0" variance to one.
    std_var_X <- ifelse(std_var_X == 0, 1, std_var_X)
    
  }else if (standardize == 'none'){
    
    std_mean_X <- rep(0, ncol(kernel_X))
    std_var_X <- rep(1, ncol(kernel_X))
    std_whiten_X <- Diagonal(n = ncol(kernel_X))
    
  }
  # demean X
  kernel_X <- sweep(kernel_X, 2, std_mean_X, FUN = "-")
  # standardize
  kernel_X <- kernel_X %*% Diagonal(x = 1/sqrt(std_var_X)) %*% std_whiten_X
  kernel_X <- as.matrix(kernel_X)  
  
  # Get some metadata
  N <- nrow(kernel_X)
  P <- ncol(kernel_X)
  
  if (is.null(bandwidth)){
    bandwidth <- 2 * P
  }
  
  
  # Create the sketching matrix S
  if (sketch_method == 'gaussian'){
    
    S <- matrix(rnorm(N * sketch_size), nrow = N)
    S <- S * 1/sqrt(sketch_size)
    
  }else if (sketch_method == 'nystrom'){
    
    # Sample rows from the design
    S <- t(as.matrix(sparseMatrix(i = 1:sketch_size,
          j = sample(1:N, sketch_size, replace = T),
          x = 1, dims = c(sketch_size, N))))
    
  }else if (sketch_method == 'bernoulli'){
    if (is.null(sketch_prob)){stop('sketch method "bernoulli" requires a probability.')}
    
    S <- matrix(rbinom(N * sketch_size, 1, prob = sketch_prob), nrow = N)
    S <- S * 1/sqrt(sketch_prob * (1 - sketch_prob))
    S <- S * 1/sqrt(sketch_size)
    
  }else if (sketch_method == 'none'){
    S <- as.matrix(Diagonal(n = N))
  }
  
  S <- as.matrix(S)
  
  # Create the sketched kernel for the *training* data
  KS <- create_sketched_kernel(X_test = kernel_X, X_train = kernel_X, 
      tS = t(S), bandwidth = bandwidth)
  # Get the eigen-decomposition of the penalty 
  eigen_sketch <- eigen(t(KS) %*% S)
  # Set all negative eigenvalues to zero. This should only occur
  # because of numerical errors as S^T K S should be positive semi-definite...
  eigen_sketch$values <- ifelse(eigen_sketch$values < control$truncate_eigen, 0, eigen_sketch$values)
  # Note the positions of the nonzero eigenvalues. Zero eigenvalues can be 
  # ignored from the analysis.
  nonzero <- which(eigen_sketch$values != 0)
  
  # Create and scale the design given to the mixed effect model.
  
  projected_data <- KS %*% eigen_sketch$vectors[,nonzero]
  projected_data <- as.matrix(projected_data %*% 
    Diagonal(x = 1/sqrt(eigen_sketch$values[nonzero])))

  # Fit the KRLS model using lme4
  fit_krls <- internal_fit_krls(
    kernel_data = projected_data, 
    formula = formula,
    data = data,
    init_var = control$init_var,
    intercept = control$intercept, family = family)
  
  # Adjust fitted values
  fit_krls$fitted <- fit_krls$fitted * sd_y + mean_y
  # Adjust RE to be dimension of sketched data.
  fit_krls$re$mean <- eigen_sketch$vectors %*% 
    c(fit_krls$re$mean * sqrt(1/eigen_sketch$values[nonzero]),
      rep(0, nrow(eigen_sketch$vectors) - length(fit_krls$re$mean)))
  # Adjust RE variance to match sketched data
  diag_ev <- Diagonal(x = sqrt(1/eigen_sketch$values[nonzero]))
  
  n_FE_p <- length(fit_krls$fe$mean)
  ev <- eigen_sketch$vectors
  ev <- bdiag(Diagonal(x = rep(1, n_FE_p)), ev)
  aug_diag_ev <- bdiag(Diagonal(x = rep(1, n_FE_p)), diag_ev)
  
  meat_ridge <- bdiag(aug_diag_ev %*% fit_krls$vcov_ridge %*% aug_diag_ev, 
                      Diagonal(x = rep(0, nrow(eigen_sketch$vectors) + 
                                         -length(nonzero))))

  fit_krls$vcov_ridge <- ev %*% meat_ridge %*% t(ev)
  fit_krls$re$var <- eigen_sketch$vectors %*% 
    bdiag(diag_ev %*% fit_krls$re$var %*% diag_ev, 
        Diagonal(x = rep(0, nrow(eigen_sketch$vectors) + 
          -length(nonzero)))) %*% t(eigen_sketch$vectors)

  fit_krls$control <- control
  fit_krls$internal <- list(sd_y = sd_y, mean_y = mean_y, 
                            bandwidth = bandwidth,
                            kernel_X_train = kernel_X,
                            eigen_sketch = eigen_sketch,
                            eigen_nonzero = nonzero,
                            sketch = S,
                            std_train = list(
                               mean = std_mean_X,
                               var = std_var_X, 
                               whiten = std_whiten_X),
                            family = family)
  
  
  class(fit_krls) <- 'gKRLS'
  
  return(fit_krls)
}

#' @importFrom lme4 fixef ranef
#' @export
fixef.gKRLS <- function(object){
  object$fe$mean
}

#' @export
ranef.gKRLS <- function(object, type = 'mean'){
  if (type == 'mean'){
    object$re$mean
  }else if (type == 'variance'){
    object$re$var
  }else{stop("type must be mean or variance")}
}

# Load fixef, ranef from lme4
#' @export
lme4::fixef
#' @export
lme4::ranef

#' @importFrom lme4 glFormula lFormula mkLmerDevfun mkGlmerDevfun mkMerMod
#'   VarCorr optimizeGlmer updateGlmerDevfun optimizeLmer
internal_fit_krls <- function(kernel_data, formula, data, y, init_var,
                              intercept, family){
  
  dim_kernel <- ncol(kernel_data)
  N <- nrow(data)
  if (N != nrow(kernel_data)){stop('FE and Kernel have different #s of observations.')}
  # Create a "fake" RE with two levels to trick lmer pre-process to run.
  kernel_RE <- c(1, rep(2, N-1))
  formula <- update.formula(formula, '. ~ . + (1 | kernel_RE)')
  data$kernel_RE <- kernel_RE
  # Generate the data  
  lmod <- glFormula(formula = formula, data = data, family = family)
  
  # Process and wrap the data from lmer
  # This reasonably closely follows: 
  # https://bbolker.github.io/mixedmodels-misc/notes/varmats.html
  mod_flist <- lmod$reTrms$flist
  mod_flist$kernel_RE[-1] <- 1
  mod_flist$kernel_RE <- factor(mod_flist$kernel_RE) 
  
  # Set initial variance to give to optimizer
  fake_reTrms <- list(
    Zt = drop0(t(kernel_data)),
    theta = init_var,
    Lambdat = sparseMatrix(i = 1:dim_kernel, j = 1:dim_kernel, x = init_var),
    Lind = rep(1, dim_kernel),
    flist = mod_flist,
    lower = lmod$reTrms$lower,
    Gp = as.integer(c(0, dim_kernel)),
    cnms = lmod$reTrms$cnms
  )
  # Prepare the lme4 object
  man_lmod <- list(
    fr = lmod$fr,
    X = lmod$X,
    family = lmod$family,
    formula = NA,
    reTrms = fake_reTrms,
    wmsgs = lmod$wmsgs
  )
  
  if (family$family == 'gaussian'){
    
    devfun <- do.call(mkLmerDevfun, man_lmod)
    opt_gKRLS <- optimizeLmer(devfun)

  }else{
    # GHrule <- lme4:::GHrule
    devfun <- do.call(mkGlmerDevfun, man_lmod)
    opt_gKRLS <- optimizeGlmer(devfun)
    devfun <- updateGlmerDevfun(devfun, man_lmod$reTrms)
    opt_gKRLS <- optimizeGlmer(devfun, stage=2)
  }
  
  # Turn into a lme4 style object
  fmt_gKRLS <- mkMerMod(environment(devfun), opt_gKRLS, 
              man_lmod$reTrms, fr = man_lmod$fr)
  
  fmla <- formula(fmt_gKRLS)
  
  # Custom internal function to extract the REs
  extract_re_KRLS <- function(object){
    ans <- object@pp$b(1)
    ans.var <- lme4:::condVar(fmt_gKRLS, scaled = TRUE)
    return(list(mean = ans, var = ans.var))
  }
  # Extract the REs
  re_all <- extract_re_KRLS(fmt_gKRLS)
  
  fe_all <- list(mean = fixef(fmt_gKRLS))
  if (length(fe_all$mean) > 0){
    fe_all$vcov <- vcov(fmt_gKRLS)
  }
  
  vcov_ridge <- get_vcov_ridge(fmt_gKRLS, family)
  
  out <- list(sigma = sigma(fmt_gKRLS),
    par_gKRLS = opt_gKRLS$par,
    fitted = fitted(fmt_gKRLS),
    fmt_varcorr = VarCorr(fmt_gKRLS),
    fe = fe_all,
    re = re_all,
    vcov_ridge = vcov_ridge,
    formula = fmla
  )

  return(out)
}

#' Internal Standard Errors
#' 
get_vcov_ridge <- function(object, family){
  
  # Get various data components
  
  X <- getME(object, 'X')
  Z <- getME(object, 'Z')
  XZ <- drop0(cbind(X,Z))
  y <- getME(object, 'y')
  mu <- getME(object, 'mu')
  re_mean <- getME(object, 'b')
  fe_mean <- getME(object, 'beta')  
  sigma <- sigma(object)
  Lambda <- getME(object, 'Lambda')
  # Add the flat prior on the FE
  Ridge <- bdiag(Diagonal(x = rep(0, ncol(X))), solve(Lambda))
  Ridge <- crossprod(Ridge)
  # XB + ZA - linear predictor of all terms
  eta <- as.vector(X %*% fe_mean) + as.vector(Z %*% re_mean)
  
  if (family$family == 'binomial'){
    if (!(family$link %in% c('logit', 'probit'))){
      stop('Standard errors not set up for binomial that is not logit/probit.')
    }
  }else if (family$family == 'poisson'){
    if (family$link != 'log'){
      stop('Standard errors not set up for non-canonical Poisson.')
    }
  }else if (family$family == 'gaussian'){
    if (family$link != 'identity'){
      stop('Standard errors not set up for non-identity link Gaussian.')
    }
  }
  
  if (family$family == 'binomial' & family$link == 'logit'){
    p <- plogis(eta) 
    weight <- p * (1-p)
    stopifnot(all.equal(mu, plogis(eta)))
    stopifnot(sigma == 1)
  }else if (family$family == 'gaussian' & family$link == 'identity'){
    weight <- rep(1, nrow(XZ))
  }else if (family$family == 'binomial' & family$link == 'probit'){
    weight <- (2 * y - 1) * dnorm(eta)/pnorm(eta * (2 * y - 1))
    weight <- weight * (eta + weight)
  }else if (family$family == 'poisson' & family$link == 'log'){
    weight <- exp(eta)
  }else{
    stop('Not recogized family and link.')
  }
  
  vcov_ridge <- sigma^2 * solve(crossprod(Diagonal(x = sqrt(weight)) %*% XZ) + Ridge)
  return(vcov_ridge)
}

