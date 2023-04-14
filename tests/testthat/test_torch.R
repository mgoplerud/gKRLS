# Test to only be run when manually called
# as it depends on torch that is not a requires/suggestion
# for gKRLS as it is only used in this test.
# Thus, it is excluded using .Rbuildignore and .Rinstignore
# but would be tested using devtools::test()

if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  env_test <- "CRAN"
  set.seed(140) # CRAN SEED
}else{
  # If on local machine
  env_test <- 'local'
}

if (interactive() & (env_test == 'local')){

  test_that(" Test AME versus automatic differentation using torch ", {
    
    require(torch)
    # From mgcv
    n <- 1000
    f1 <- function(x) sin(3*pi*x)*exp(-x)
    f2 <- function(x) x^3
    f3 <- function(x) .5*exp(-x^2)-.2
    f4 <- function(x) 1
    x1 <- runif(n);x2 <- runif(n)
    eta1 <- 2*(f1(x1) + f2(x2))-.5
    eta2 <- 2*(f3(x1) + f4(x2))-1
    p <- exp(cbind(0,eta1,eta2))
    p <- p/rowSums(p) ## prob. of each category 
    cp <- t(apply(p,1,cumsum)) ## cumulative prob.
    ## simulate multinomial response with these probabilities
    ## see also ?rmultinom
    y <- apply(cp,1,function(x) min(which(x>runif(1))))-1
    y[sample(1:1000, 100)] <- 3
    dat <- data.frame(x1, x2, y = y + 1)
    dat$new_y <- dat$y - 1
    rm(x1, x2, y)
    
    # From mgcv: Confirm predictions line up for
    # single linear predictor models
    
    for (f in list(multinom(K=3), multinom(K=2), ocat(R=4), binomial(link = 'cloglog'), 
                   binomial(link = 'logit'), binomial(link = 'probit'), 
                   poisson(), gaussian(), tw(), nb(), scat(), ziP())){
      # print(f)
      if (grepl(f$family, pattern='multinom')){
        if (f$nlp == 3){
          b <- gam(list(new_y ~ te(x1, x2),
                        ~ te(x1,x2) + s(x1,x2),
                        ~ s(x1, x2, bs = 'gKRLS')), 
                   data = dat,
                   family = multinom(K = 3))
        }else{
          new_dat <- subset(dat, new_y < 3)
          b <- gam(list(new_y ~ te(x1, x2), 
                        ~ s(x1, x2, bs = 'gKRLS')), 
                   data = new_dat,
                   family = multinom(K = 2))
        }
      }else if (!grepl(f$family, pattern='^Zero|^binomial$')){
        b <- suppressWarnings(gam(y~ te(x1,x2),family=f,data=dat))
      }else if (grepl(f$family, pattern='binomial')){
        b <- suppressWarnings(gam(as.numeric(y > 1) ~ te(x1,x2), family=f,data=dat))
      }else{
        b <- suppressWarnings(gam(new_y~ te(x1,x2),family=f,data=dat))
      }
      pred_gKRLS <- calculate_effects(model = b, continuous_type = 'predict')
      pred_mgcv <- predict(b, type = 'response')
      
      torch_X <- torch_tensor(model.matrix(b))
      torch_beta <- torch_tensor(matrix(coef(b)), requires_grad = TRUE)
      torch_lp <- torch_mm(torch_X, torch_beta)
      
      if (f$family == 'multinom'){
        lpi <- attr(model.matrix(b), 'lpi')
        torch_lp <- lapply(lpi, FUN=function(i){
          torch_mm(torch_X[,i], torch_beta[i,])
        })
        torch_lp <- c(list(torch_tensor(matrix(rep(0, nrow(model.matrix(b)))))), torch_lp)
        torch_eta <- torch_cat(torch_lp, dim = 2)
        torch_p <- nn_softmax(2)(torch_eta)
      }else if (f$family == 'zero inflated Poisson'){
        torch_theta <- torch_tensor(b$family$getTheta(TRUE))
        torch_eta <- torch_theta[1] + torch_theta[2] * torch_lp
        torch_p_nonzero <- 1 - torch_exp(-torch_exp(torch_eta))
        torch_lambda <- torch_exp(torch_lp)
        torch_mu_nonzero <- torch_lambda/(1 - torch_exp(-torch_lambda))
        torch_p <- torch_p_nonzero * torch_mu_nonzero
      }else if (f$family == 'binomial'){
        if (f$link == 'logit'){
          torch_p <- torch_sigmoid(torch_lp)
        }else if (f$link == 'probit'){
          torch_ncdf <- distr_normal(loc = 0, scale = 1)
          torch_p <- torch_ncdf$cdf(torch_lp)
        }else if (f$link == 'cloglog'){
          torch_p <- 1 - torch_exp(-torch_exp(torch_lp))
        }else{stop('Unsetup link')}
      }else if (f$family == 'Ordered Categorical'){
        torch_theta <- torch_tensor(c(-Inf, b$family$getTheta(TRUE), Inf))
        torch_cumul_p <- torch_sigmoid(torch_add(-torch_lp, torch_theta))
        torch_p <- torch_diff(torch_cumul_p)
      }else if (f$family %in% c('poisson', 'Tweedie', 'negative binomial')){
        torch_p <- torch_exp(torch_lp)
      }else if (f$family %in% c('gaussian', 'scaled t')){
        torch_p <- torch_lp
      }else{stop('.')}
      
      torch_pred <- torch_mean(torch_p, dim = 1)
      est_se_torch <- sapply(1:nrow(torch_pred), FUN=function(i){
        grad_beta <- as_array(autograd_grad(torch_pred[i], torch_beta, retain_graph = TRUE)[[1]])
        sqrt(as.numeric(t(grad_beta) %*% vcov(b) %*% grad_beta))
      })
      est_torch <- as_array(torch_pred)
      
      if (is.matrix(pred_mgcv)){
        avg_mgcv <- colMeans(pred_mgcv)
      }else{
        avg_mgcv <- mean(pred_mgcv)
      }
      expect_equal(avg_mgcv, pred_gKRLS$est)
      expect_equal(est_torch, pred_gKRLS$est, tol = 1e-6, scale = 1)
      expect_equal(est_se_torch, pred_gKRLS$se, tol = 1e-6, scale = 1)
      rm(pred_gKRLS, pred_mgcv)
      rm(list = ls(pattern='^torch_'))
    }
  })
  
}
