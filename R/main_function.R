#' gKRLS
#' 
#' Fit a generalized KRLS model using Laplace approximation
#' @param y A numeric vector that contains the values of dependent variable or response variable. Missing Values not allowed.
#' @param X A matrix of numeric observations of independent variables. Factors and missing values not allowed. Also, no intercept is required.
#' @param family Provide a family as in the Generalized Linear Models (GLM). For example: gaussian or poisson.
#' @param sketch_size Specify the dimension of the randomized sketches of the kernel matrix. The default size is ``ceiling(nrow(X)^(1/3)) * 5,'' where X is matrix of the independent variable.
#' @param sketch_method String vector that specifies which kernel sketch method should be used. Options include ``gaussian'': gaussian kernel approximation, 
#' ``nystrom'': , bernoulli, and no sketch. The default is gaussian.
#' @param bandwidth A numerical number indicating the bandwidth for kernel. The default is  ``2*ncol(X).'' 
#' @param standardize A string vector that specifies which standardization method should be used. Must be one of scaled, Mahalanobis, or none, which menas no standardization. The default is ``Mahalanobis''.
#' @param control A list containing control parameters. The init_var specify the initial value for covariance parameter in lme4 models. See details about theta in lme4 package.
#' The default is 1. The truncate_eigen determines how much eigenvalues should be truncated. If truncate_eigen set to 1e-6, this means we only keep eigenvalue greater or equal to 1e-6.
#' The default is sqrt(.Machine$double.eps), where .Machine$double.eps is the smallest positive floating-point number x such that 1 + x != 1.
#' @param intercept Specify whether add intercept during lme4 estimation. 
#' @useDynLib gKRLS
#' @import Matrix
#' @export
gKRLS <- function(y, X, family, sketch_size = ceiling(nrow(X)^(1/3)) * 5, 
      sketch_method = 'gaussian', sketch_prob = NULL, bandwidth = NULL,
      standardize = 'Mahalanobis',
      control = list(init_var =  1, 
                     truncate_eigen = sqrt(.Machine$double.eps)),
      intercept = FALSE){
  
  if (!inherits(family, 'family')){
    stop('Provide "family" as in glm, glmer, etc.')
  }
  if (intercept){stop('Does not accept intercept, yet..')}

  
  if (family$family == 'gaussian'){
    sd_y <- sd(y)
    mean_y <- mean(y)
    fmt_y <- (y - mean(y))/sd(y)
  }else{
    sd_y <- 1
    mean_y <- 0
    fmt_y <- y
  }
  
  if (standardize == 'Mahalanobis'){
    
    std_mean_X <- colMeans(X)
    std_var_X<- rep(1, ncol(X))
    # var(X) = A -> var(X A) = A^T var(X) A 
    # A Q L Q^T A -> A = Q sqrt(L)^{-1}
    std_whiten_X <- with(eigen(cov(X)), 
          vectors %*% Diagonal(x = ifelse(values < sqrt(.Machine$double.eps), 0, 1/sqrt(values))))
    zero_columns <- apply(std_whiten_X, MARGIN = 2, FUN=function(i){all(i == 0)})
    zero_columns <- which(zero_columns == TRUE)  
    if (length(zero_columns) > 0){
      print(paste('Original X matrix is not full rank. Using an effective rank of', ncol(X) - length(zero_columns), 'throughout.'))
      std_whiten_X <- std_whiten_X[,-zero_columns,drop=F]
    }
  }else if (standardize == 'scaled'){
    
    std_mean_X <- colMeans(X)
    std_var_X <- apply(X, MARGIN = 2, var)
    std_whiten_X <- Diagonal(n = ncol(X)) 
    
    # Set "0" variance to one.
    std_var_X <- ifelse(std_var_X == 0, 1, std_var_X)
    
  }else if (standardize == 'none'){
    
    std_mean_X <- rep(0, ncol(X))
    std_var_X <- rep(1, ncol(X))
    std_whiten_X <- Diagonal(n = ncol(X))
    
  }
  # demean X
  X <- sweep(X, 2, std_mean_X, FUN = "-")
  # standardize
  X <- X %*% Diagonal(x = 1/sqrt(std_var_X)) %*% std_whiten_X
  X <- as.matrix(X)  
  
  # Get some metadata
  N <- nrow(X)
  P <- ncol(X)
  
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
  KS <- create_sketched_kernel(X_test = X, X_train = X, 
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
  fit_krls <- internal_fit_krls(data = projected_data, 
    y = fmt_y, N = N, init_var = control$init_var,
    intercept = intercept, family = family)
  
  # Adjust fitted values
  fit_krls$fitted <- fit_krls$fitted * sd_y + mean_y
  # Adjust RE to be dimension of sketched data.
  fit_krls$re$mean <- eigen_sketch$vectors %*% 
    c(fit_krls$re$mean * sqrt(1/eigen_sketch$values[nonzero]),
      rep(0, nrow(eigen_sketch$vectors) - length(fit_krls$re$mean)))
  # Adjust RE variance to match sketched data
  diag_ev <- Diagonal(x = sqrt(1/eigen_sketch$values[nonzero]))
  
  fit_krls$re$var <- eigen_sketch$vectors %*% 
    bdiag(diag_ev %*% fit_krls$re$var %*% diag_ev, 
        Diagonal(x = rep(0, nrow(eigen_sketch$vectors) + 
          -length(nonzero)))) %*% t(eigen_sketch$vectors)

  fit_krls$control <- control
  fit_krls$internal <- list(sd_y = sd_y, mean_y = mean_y, 
                            bandwidth = bandwidth,
                            X_train = X,
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
internal_fit_krls <- function(data, y, N, init_var,
                              intercept, family){
  
  dim_kernel <- ncol(data)
  # Create a "fake" RE with two levels to trick lmer to running.
  kernel_RE <- c(1, rep(2, N-1))
  
  if (intercept){
    lmod <- glFormula(y ~ 1 + (1 | kernel_RE), family = family)
  }else{
    lmod <- glFormula(y ~ 0 + (1 | kernel_RE), family = family)
  }
  # Process and wrap the data from lmer
  # This reasonably closely follows: 
  # https://bbolker.github.io/mixedmodels-misc/notes/varmats.html
  mod_flist <- lmod$reTrms$flist
  mod_flist$kernel_RE[-1] <- 1
  mod_flist$kernel_RE <- factor(mod_flist$kernel_RE) 
  
  # Set initial variance to give to optimizer
  fake_reTrms <- list(
    Zt = drop0(t(data)),
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

    # # Custom implementation...
    # devfun <- do.call(mkLmerDevfun, man_lmod)
    # # Set an objective function that calls this on the *log scale*
    # # for somewhat more stable fitting
    # eval_gKRLS <- function(log_theta){
    #   devfun(exp(log_theta))
    # }
    # # Call optim using 'L-BFGS-B' with boundary constraints.
    # opt_gKRLS <- optim(par=log(init_var), eval_gKRLS, method = 'L-BFGS-B',
    #                    lower = -5, upper = 5)
    # # Wrap the output together
    # opt_gKRLS <- with(opt_gKRLS,
    #                   list(par=exp(par), 
    #                        fval=value, conv=convergence, message=message))
  }else{
    GHrule <- lme4:::GHrule
    devfun <- do.call(mkGlmerDevfun, man_lmod)
    opt_gKRLS <- optimizeGlmer(devfun)
    devfun <- updateGlmerDevfun(devfun, man_lmod$reTrms)
    opt_gKRLS <- optimizeGlmer(devfun, stage=2)

  }
  # Turn into a lme4 style object
  fmt_gKRLS <- mkMerMod(environment(devfun), opt_gKRLS, 
              man_lmod$reTrms, fr = man_lmod$fr)
  # Custom internal function to extract the REs
  extract_re_KRLS <- function(object){
    ans <- object@pp$b(1)
    ans.var <- lme4:::condVar(fmt_gKRLS, scaled = TRUE)
    return(list(mean = ans, var = ans.var))
  }
  # Extract the REs
  re_all <- extract_re_KRLS(fmt_gKRLS)
  
  out <- list(sigma = sigma(fmt_gKRLS),
              var = opt_gKRLS$par,
              fitted = fitted(fmt_gKRLS),
              fmt_varcorr = VarCorr(fmt_gKRLS),
              fe = list(mean = fixef(fmt_gKRLS), vcov = vcov(fmt_gKRLS)),
              re = re_all
  )
  
  return(out)
}