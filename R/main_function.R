#' gKRLS
#' 
#' Fit a generalized KRLS model using Laplace approximation
#' @param formula A formula containing the outcome and any fixed effects. Use `y
#'   ~ 1` for a model with only the kernel and an intercept.
#' @param data data.frame containing the outcome and any fixed effects
#' @param kernel_X matrix containing the variables used to build the kernel 
#' @param X A matrix of numeric observations of independent variables. Factors
#'   and missing values not allowed. Also, no intercept is required.
#' @param family Provide a family as in the Generalized Linear Models (GLM). For
#'   example: gaussian or poisson.
#' @param sketch_size Specify the dimension of the randomized sketches of the
#'   kernel matrix. The default size is ``ceiling(nrow(X)^(1/3)) * 5,'' where X
#'   is matrix of the independent variable.
#' @param sketch_method String vector that specifies which kernel sketch method
#'   should be used. Options include ``gaussian'': gaussian kernel
#'   approximation, ``nystrom'': , bernoulli, and no sketch. The default is
#'   gaussian.
#' @param bandwidth A numerical number indicating the bandwidth for kernel. The
#'   default is  ``2*ncol(X).''
#' @param standardize A string vector that specifies which standardization
#'   method should be used. Must be one of scaled, Mahalanobis, or none, which
#'   menas no standardization. The default is ``Mahalanobis''.
#' @param control A list containing control parameters. The init_var specify the
#'   initial value for covariance parameter in lme4 models. See details about
#'   theta in lme4 package. The default is 1. The truncate_eigen determines how
#'   much eigenvalues should be truncated. If truncate_eigen set to 1e-6, this
#'   means we only keep eigenvalue greater or equal to 1e-6. The default is
#'   sqrt(.Machine$double.eps), where .Machine$double.eps is the smallest
#'   positive floating-point number x such that 1 + x != 1.
#' @useDynLib gKRLS
#' @importFrom RcppParallel RcppParallelLibs
#' @import Matrix
#' @importFrom lme4 subbars 
#' @export
gKRLS <- function(
      formula, data,
      kernel_X, family, 
      sketch_size = ceiling(nrow(kernel_X)^(1/3)) * 5, 
      sketch_method = 'gaussian', sketch_prob = NULL, bandwidth = NULL,
      standardize = 'scaled', prior_stabilize = TRUE,
      verbose = FALSE,
      control = list(init_var =  1,
                     truncate_eigen = sqrt(.Machine$double.eps))){
  
  if (!inherits(family, 'family')){
    stop('Provide "family" as in glm, glmer, etc.')
  }
  
  # Parse data and get response
  data <- model.frame(subbars(formula), data)
  nobs_complete <- nrow(data)
  response <- model.response(data)
  
  # From lm to deal with factors
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), 
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$formula <- nobars(formula)
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  response_2 <- model.response(mf, "numeric")
  design_FE <- model.matrix(mt, mf, contrasts)
  xlevels <- .getXlevels(mt, mf)
  
  fe_design_options <- list(contrasts = attr(design_FE, 'contrasts'),
       xlevels = xlevels,
       terms = mt,
       na.action = attr(mf, 'na.action'))
  rm(xlevels, mt, mf, m); gc()
  
  if (!isTRUE(all.equal(response_2, response))){
    stop('response parsed incorrectly. Maybe NA in data?')
  }

  if (family$family == 'gaussian'){
    sd_y <- sd(response)
    mean_y <- mean(response)
    fmt_y <- (response - mean_y)/sd_y
    # data[[1]] <- fmt_y
    data[[1]] <- NA
    response <- fmt_y
    rm(fmt_y)
  }else{
    sd_y <- 1
    mean_y <- 0
    
  }
  
  if (verbose){
    print('Standardizing Data')
    print(Sys.time())
  }
  
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
  
  # Get some metadata
  N <- nrow(std_kernel_X)
  P <- ncol(std_kernel_X)
  
  if (is.null(bandwidth)){
    bandwidth <- P
  }
  
  if (verbose){
    print('Create Sketching Matrix')
    print(Sys.time())
  }
  
  full_kernel_original <- std_kernel_X
  # Create the sketching matrix S
  if (sketch_method == 'gaussian'){
    
    S <- matrix(rnorm(N * sketch_size), nrow = N)
    S <- S * 1/sqrt(sketch_size)
    
  }else if (sketch_method == 'nystrom'){
    
    nystrom_id <- sample(1:N, sketch_size, replace = T)

    std_kernel_X <- full_kernel_original[nystrom_id,, drop = F]
    S <- as.matrix(Diagonal(n = length(nystrom_id)))
    
  }else if (sketch_method == 'bernoulli'){
    if (is.null(sketch_prob)){stop('sketch method "bernoulli" requires a probability.')}
    
    S <- matrix(rbinom(N * sketch_size, 1, prob = sketch_prob), nrow = N)
    S <- S * 1/sqrt(sketch_prob * (1 - sketch_prob))
    S <- S * 1/sqrt(sketch_size)
    
  }else if (sketch_method == 'none'){
    S <- as.matrix(Diagonal(n = N))
  }
  
  S <- as.matrix(S)
  
  if (verbose){
    print('Create Sketched Kernel')
    print(Sys.time())
  }
  # Create the sketched kernel for the *training* data
  KS <- create_sketched_kernel(
      X_test = full_kernel_original, 
      X_train = std_kernel_X, 
      tS = t(S), bandwidth = bandwidth)
  
  if (verbose){
    print('Eigendecompose Penalty')
    print(Sys.time())
  }
  # Get the eigen-decomposition of the penalty 
  if (sketch_method == 'nystrom'){
    eigen_sketch <- eigen(t(KS[nystrom_id,]) %*% S)
  }else{
    eigen_sketch <- eigen(t(KS) %*% S)
  }
  # Set all negative eigenvalues to zero. This should only occur
  # because of numerical errors as S^T K S should be positive semi-definite...
  eigen_sketch$values <- ifelse(eigen_sketch$values < control$truncate_eigen, 0, eigen_sketch$values)
  # Note the positions of the nonzero eigenvalues. Zero eigenvalues can be 
  # ignored from the analysis.
  nonzero <- which(eigen_sketch$values != 0)
  
  # Create and scale the design given to the mixed effect model.
  if (verbose){
    print('Project Sketched Data')
    print(Sys.time())
  }
  projected_data <- KS %*% eigen_sketch$vectors[,nonzero]
  projected_data <- as.matrix(projected_data %*% 
    Diagonal(x = 1/sqrt(eigen_sketch$values[nonzero])))

  if (verbose){
    print('Fit lmer')
    print(Sys.time())
  }
  
  # Fit the KRLS model using lme4
  fit_krls <- internal_fit_krls(
    kernel_data = projected_data, 
    response = response,
    design_FE = design_FE,
    init_var = control$init_var,
    family = family,
    formula = formula,
    prior_stabilize = prior_stabilize)

  d_j <- fit_krls$d_j
  if (names(d_j)[length(d_j)] != 'kernel_RE'){
    stop('kernel_RE must come last')
  }
  n_normal_RE <- sum(d_j[-length(d_j)])
  n_kernel <- d_j['kernel_RE']
  # Adjust fitted values
  fit_krls$fitted <- fit_krls$fitted * sd_y + mean_y
  
  if (verbose){
    print('Process lmer')
    print(Sys.time())
  }
  
  # Adjust RE to be dimension of sketched data.
  fit_krls$re$mean <- bdiag(Diagonal(n = n_normal_RE), eigen_sketch$vectors) %*% 
    c(fit_krls$re$mean[seq_len(n_normal_RE)], 
      fit_krls$re$mean[n_normal_RE + seq_len(n_kernel)] * sqrt(1/eigen_sketch$values[nonzero]),
      rep(0, nrow(eigen_sketch$vectors) - n_kernel))
  # Adjust RE variance to match sketched data
  diag_ev <- Diagonal(x = sqrt(1/eigen_sketch$values[nonzero]))
  
  n_FE_p <- length(fit_krls$fe$mean)
  ev <- eigen_sketch$vectors
  ev <- bdiag(Diagonal(x = rep(1, n_FE_p + n_normal_RE)), ev)
  aug_diag_ev <- bdiag(Diagonal(x = rep(1, n_FE_p + n_normal_RE)), diag_ev)
  
  
  meat_ridge <- bdiag(aug_diag_ev %*% fit_krls$vcov_ridge$vcov_ridge %*% aug_diag_ev, 
                      Diagonal(x = rep(0, nrow(eigen_sketch$vectors) + 
                                         -length(nonzero))))

  effective_N <- nrow(kernel_X) - fit_krls$vcov_ridge$effective_df
  

  fit_krls$effective_df <- fit_krls$vcov_ridge$effective_df
  fit_krls$effective_N <- effective_N
  fit_krls$ll <- fit_krls$vcov_ridge$ll
  fit_krls$vcov_ridge <- ev %*% meat_ridge %*% t(ev)
  int_aug_ev <- bdiag(Diagonal(x = rep(1, n_normal_RE)), diag_ev)
  fit_krls$re$var <- bdiag(Diagonal(n = n_normal_RE), eigen_sketch$vectors) %*% 
    bdiag(int_aug_ev %*% fit_krls$re$var %*% int_aug_ev, 
        Diagonal(x = rep(0, nrow(eigen_sketch$vectors) + 
          -length(nonzero)))) %*% t(bdiag(Diagonal(n = n_normal_RE), eigen_sketch$vectors))

  if (sketch_method == 'nystrom'){
    kernel_X <- kernel_X[nystrom_id,]
  }
  fit_krls$control <- control
  fit_krls$formula <- formula
  fit_krls$internal <- list(sd_y = sd_y, mean_y = mean_y, 
                            bandwidth = bandwidth,
                            kernel_X_train = kernel_X,
                            eigen_sketch = eigen_sketch,
                            eigen_nonzero = nonzero,
                            sketch = S,
                            fe_options = fe_design_options,
                            standardize = standardize,
                            response = response,
                            std_train = list(
                               mean = std_mean_X,
                               whiten = std_whiten_X,
                               W_Matrix = W_Matrix),
                            family = family)
  
  if (verbose){
    print('Completed')
    print(Sys.time())
  }
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
#'   VarCorr optimizeGlmer updateGlmerDevfun optimizeLmer nloptwrap GHrule
internal_fit_krls <- function(kernel_data, design_FE, response, 
  init_var, family, formula, prior_stabilize){
  
  copy_init <- as.numeric(init_var + 1 - 1)
  dim_kernel <- ncol(kernel_data)
  N <- nrow(design_FE)
  if (N != nrow(kernel_data)){stop('FE and Kernel have different #s of observations.')}
  # Create a "fake" RE with two levels to trick lmer pre-process to run.
  kernel_RE <- c(1, rep(2, N-1))

  # formula <- as.formula(paste0(paste(deparse(formula), collapse=' '), ' + (1 | kernel_RE)'))
  # data$kernel_RE <- kernel_RE
  # # Generate the data  
  # lmod <- glFormula(formula = formula, data = data, family = family)
  
  extra_re <- paste(sapply(findbars(formula), FUN=function(i){paste0('(', deparse(i), ')')}), collapse = ' + ')
  if (extra_re == ""){extra_re <- NULL}
  
  if (ncol(design_FE) == 0){
    formula <- as.formula(paste(c('response ~ 0 + (1 | kernel_RE)', extra_re), collapse = ' + '), env = environment())
  }else{
    formula <- as.formula(paste(c('response ~ 0 + design_FE + (1 | kernel_RE)', extra_re), collapse = ' + '), env = environment())
  }
  lmod <- glFormula(formula = formula, data = NULL, family = family) 
  if (ncol(design_FE) != 0 & max(abs(lmod$X - design_FE)) != 0){
    stop('lmer processing of FE failed.')
  }
  colnames(lmod$X) <- gsub(colnames(lmod$X), pattern='^design_FE', replacement = '')
  if (ncol(lmod$X) == 1){
    if (colnames(lmod$X)[1] == ""){
      colnames(lmod$X) <- '(Intercept)'
    }
  }
  # Process and wrap the data from lmer
  # This reasonably closely follows: 
  # https://bbolker.github.io/mixedmodels-misc/notes/varmats.html
  mod_flist <- lmod$reTrms$flist
  mod_flist$kernel_RE[-1] <- 1
  mod_flist$kernel_RE <- factor(mod_flist$kernel_RE) 
  
  
  d_j <- diff(lmod$reTrms$Gp)
  names(d_j) <- names(lmod$reTrms$cnms)
  
  if (names(lmod$reTrms$cnms)[length(d_j)] != 'kernel_RE'){
    stop("kernel_RE isn't at the end")
  }
  d_j[which(names(lmod$reTrms$cnms) == 'kernel_RE')] <- dim_kernel
  Gp <- cumsum(c(0, d_j))
  
  Zt_aug <- lmod$reTrms$Ztlist
  Zt_aug[["1 | kernel_RE"]] <- drop0(t(kernel_data))
  Zt_aug <- do.call('rbind', Zt_aug)
  
  theta_aug <- lmod$reTrms$theta
  theta_aug[which(names(lmod$reTrms$cnms) == 'kernel_RE')] <- copy_init
  
  Lind_aug <- rep(1:length(d_j), d_j)
  Lambdat_aug <- sparseMatrix(i = 1:sum(d_j), j = 1:sum(d_j), x = rep(theta_aug, d_j))

  
  
  # Set initial variance to give to optimizer
  fake_reTrms <- list(
    Zt = Zt_aug,
    theta = theta_aug,
    Lambdat = Lambdat_aug,
    Lind = Lind_aug,
    flist = mod_flist,
    lower = lmod$reTrms$lower,
    Gp = as.integer(Gp),
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
    # Fixing variance
    # https://stackoverflow.com/questions/39718754/fixing-variance-values-in-lme4/39732683#39732683
    GHrule <- lme4::GHrule
    
    if (!prior_stabilize){
      devfun <- do.call(mkGlmerDevfun, man_lmod)
      opt_gKRLS <- optimizeGlmer(devfun)
      devfun <- updateGlmerDevfun(devfun, man_lmod$reTrms)
      opt_gKRLS <- optimizeGlmer(devfun, stage = 2)
    }else{
      devfun <- do.call(mkGlmerDevfun, man_lmod)
      theta.lwr <- environment(devfun)$lower  ## this changes after update!
      devfun <- updateGlmerDevfun(devfun, man_lmod$reTrms)
      rho <- environment(devfun)
      nbeta <- ncol(rho$pp$X)
      theta <- rho$pp$theta
      if (any(theta < 1e-6)){
        theta[which(theta < 1e-6)] <- theta_aug[which(theta < 1e-6)]
        placeholder <- devfun(c(theta, rep(0, nbeta)))
      }
      if (any(!is.finite(theta))){
        theta[which(!is.finite(theta))] <- theta_aug[which(!is.finite(theta))]
        placeholder <- devfun(c(theta, rep(0, nbeta)))
      }
      n_theta <- seq_len(length(theta_aug))
      penalized_dev <- function(par){
        par[n_theta] <- exp(par[n_theta])
        return(devfun(par) - sum(log(par[n_theta])))
      }
      
      opt_gKRLS <- tryCatch(lme4::nloptwrap(par=c(log(theta),rep(0,nbeta)),
        fn=penalized_dev, lower = rep(-Inf, length(theta_aug) + nbeta), upper = rep(Inf, length(theta_aug) + nbeta)), error = function(e){NULL})
      if (is.null(opt_gKRLS)){stop('Estimation failed')}
      opt_gKRLS$par[n_theta] <- exp(opt_gKRLS$par[n_theta])
    }
    
    
  }
  
  # Turn into a lme4 style object
  fmt_gKRLS <- mkMerMod(environment(devfun), opt_gKRLS, 
              man_lmod$reTrms, fr = man_lmod$fr)
  # fmla <- formula(fmt_gKRLS)

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
  
  levels_of_RE <- lapply(lmod$reTrms$Ztlist, rownames)
  
  out <- list(
    sigma = sigma(fmt_gKRLS),
    par_gKRLS = opt_gKRLS$par,
    fitted = fitted(fmt_gKRLS),
    fmt_varcorr = VarCorr(fmt_gKRLS),
    fe = fe_all,
    re = re_all,
    d_j = d_j,
    levels_of_RE = levels_of_RE,
    AIC = AIC(fmt_gKRLS),
    BIC = BIC(fmt_gKRLS),
    logLik = logLik(fmt_gKRLS),
    vcov_ridge = vcov_ridge
  )

  return(out)
}

#' Internal Standard Errors
#' @importFrom lme4 getME
get_vcov_ridge <- function(object, family){
  
  # Get various data components
  
  X <- getME(object, 'X')
  Z <- getME(object, 'Z')
  XZ <- drop0(cbind(X,Z))
  y <- getME(object, 'y')
  re_mean <- getME(object, 'b')
  fe_mean <- getME(object, 'beta')  
  sigma <- sigma(object)
  Lambda <- getME(object, 'Lambda')
  
  if (all(diag(Lambda) == 0)){
    if (!isDiagonal(Lambda)){stop("Unusual non-invertible Lambda")}
    Ridge <- bdiag(Diagonal(x = c(rep(0, ncol(X)), rep(Inf, ncol(Lambda)))))
  }else{
    Ridge <- bdiag(Diagonal(x = rep(0, ncol(X))), solve(Lambda))
  }

  # Add the flat prior on the FE
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
    stopifnot(sigma == 1)
    ll <- sum(ifelse(y == 1, log(p), log(1-p)))
  }else if (family$family == 'gaussian' & family$link == 'identity'){
    weight <- rep(1, nrow(XZ))
    ll <- sum(dnorm(y, mean = eta, sd = sigma, log = TRUE))
  }else if (family$family == 'binomial' & family$link == 'probit'){
    weight <- (2 * y - 1) * dnorm(eta)/pnorm(eta * (2 * y - 1))
    weight <- weight * (eta + weight)
    p <- pnorm(eta)
    ll <- sum(ifelse(y == 1, log(p), log(1-p)))
  }else if (family$family == 'poisson' & family$link == 'log'){
    weight <- exp(eta)
    ll <- sum(dpois(y, lambda = exp(eta), log = TRUE))
  }else{
    stop('Not recogized family and link.')
  }
  
  LL_Meat <- crossprod(Diagonal(x = sqrt(weight)) %*% XZ)
  vcov_ridge <- sigma^2 * solve(LL_Meat + Ridge)
  
  effective_df <- sum(Matrix::diag( vcov_ridge %*% (1/sigma^2 * LL_Meat) ))
  
  return(list(vcov_ridge = vcov_ridge, LL_Meat = LL_Meat,
              effective_df = effective_df, ll = ll))
  
}

