if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  env_test <- "CRAN"
  set.seed(125) # CRAN SEED
}else{
  # If on local machine
  env_test <- 'local'
}

test_that(" Test for prediction/SE for complex families ", {
  
  # From mgcv
  n <- 400
  dat <- gamSim(1,n=n, verbose = FALSE)
  dat$f <- dat$f - mean(dat$f)
  
  alpha <- c(-Inf,-1,0,5,Inf)
  R <- length(alpha)-1
  y <- dat$f
  u <- runif(n)
  u <- dat$f + log(u/(1-u)) 
  for (i in 1:R) {
    y[u > alpha[i]&u <= alpha[i+1]] <- i
  }
  dat$y <- y
  dat$new_y <- dat$y - 1
  
  # From mgcv: Confirm predictions line up for
  # single linear predictor models

  if (env_test == 'CRAN'){# Slightly fewer links for CRAN
    list_link <- list(ocat(R=R), poisson(link = 'identity'), nb(), scat(), ziP())
  }else{
    list_link <- list(ocat(R=R), poisson(link = 'identity'), gaussian(), tw(), nb(), scat(), ziP())
  }
  for (f in list_link){
    if (!grepl(f$family, pattern='^Zero')){
      b <- suppressWarnings(gam(y~ s(x0) + s(x1) + s(x2) + s(x3),family=f,data=dat))
    }else{
      b <- suppressWarnings(gam(new_y~ s(x0) + s(x1) + s(x2) + s(x3),family=f,data=dat))
    }
    pred_gKRLS <- calculate_effects(model = b, continuous_type = 'predict', individual = T)
    pred_mgcv <- predict(b, newdata = dat, se.fit = TRUE, type = 'response')    
    expect_equivalent(get_individual_effects(pred_gKRLS)$est, as.vector(pred_mgcv$fit))
    expect_equivalent(get_individual_effects(pred_gKRLS)$se, as.vector(pred_mgcv$se.fit))
    rm(pred_gKRLS, pred_mgcv)
  }
  
  b <- gam(list(new_y ~ s(x0), 
                ~ x0 + x1 + x2, 
                ~ s(x3, x2, bs = 'gKRLS')), 
           data = dat,
           family = multinom(K = 3))

  pred_gKRLS <- calculate_effects(model = b, continuous_type = 'predict', individual = T)
  pred_mgcv <- predict(b, newdata = dat, se.fit = TRUE, type = 'response')

  expect_equivalent(get_individual_effects(pred_gKRLS)$est, as.vector(pred_mgcv$fit))
  if (packageVersion('mgcv') > '1.8-42'){
    expect_equivalent(get_individual_effects(pred_gKRLS)$se, as.vector(pred_mgcv$se.fit), tol = 1e-6, scale = 1)
  }else{
    expect_true(!isTRUE(all.equal(get_individual_effects(pred_gKRLS)$se, as.vector(pred_mgcv$se.fit))))
  }
  
  # Test when K = 2 (i.e. 3 categories)
  new_dat <- subset(dat, new_y < 3)
  b <- gam(list(new_y ~ s(x0), 
                ~ s(x3, x2, bs = 'gKRLS')), 
           data = new_dat,
           family = multinom(K = 2))
  
  pred_gKRLS <- calculate_effects(model = b, continuous_type = 'predict', individual = T)
  pred_mgcv <- predict(b, newdata = new_dat, se.fit = TRUE, type = 'response')
  
  expect_equivalent(get_individual_effects(pred_gKRLS)$est, as.vector(pred_mgcv$fit), tol = 1e-6, scale = 1)
  expect_equivalent(get_individual_effects(pred_gKRLS)$se, as.vector(pred_mgcv$se.fit), tol = 1e-6, scale = 1)
  
  if (env_test != 'CRAN'){
    # Test that multinomial and logit agree when 2 categories (K=1) 
    new_dat_2 <- subset(dat, new_y < 2)
    b <- gam(list(new_y ~ x0 + x1 + x2 * x3), 
             data = new_dat_2,
             family = multinom(K = 1), method = 'REML')
    b2 <- gam(new_y ~ x0 + x1 + x2 * x3,
              data = new_dat_2,
              family = binomial(), method = 'REML')
    
    fit_multinom <- calculate_effects(b, variables = 'x3')
    fit_logit <- calculate_effects(b2, variables = 'x3')
    expect_equivalent(
      fit_multinom[2,-5],
      fit_logit,
      tol = 1e-4, 
    )
    expect_equivalent(vcov(b), vcov(b2), tol = 1e-4, scale = 1)
    expect_equivalent(coef(b), coef(b2), tol = 1e-4, scale = 1)
  }
})

test_that("calculate_interactions works with complex familes", {
  
  n <- 100
  dat <- gamSim(1,n=n, verbose = FALSE)
  dat$f <- dat$f - mean(dat$f)
  
  alpha <- c(-Inf,-1,0,5,Inf)
  R <- length(alpha)-1
  y <- dat$f
  u <- runif(n)
  u <- dat$f + log(u/(1-u)) 
  for (i in 1:R) {
    y[u > alpha[i]&u <= alpha[i+1]] <- i
  }
  dat$y <- y
  dat$new_y <- dat$y - 1
  dat$x0 <- as.numeric(dat$x0 > median(dat$x0))
  dat$x1 <- cut(dat$x1, 3)

  for (f in list(ocat(R=R), nb(), scat(), ziP())){
    # print(f)
    if (!grepl(f$family, pattern='^Zero')){
      b <- suppressWarnings(gam(y~ x0 * x1 + s(x2) + s(x3),family=f,data=dat))
    }else{
      b <- suppressWarnings(gam(new_y~ x0 * x1 + s(x2) + s(x3),family=f,data=dat))
    }
    # Check that flipping doesn't affect estimates
    pred_gKRLS <- calculate_interactions(model = b, variables = list(c('x0', 'x1')))
    pred_gKRLS_alt <- calculate_interactions(model = b, variables = list(c('x1', 'x0')))
    expect_equal(
      subset(pred_gKRLS[-c(2,3)], QOI %in% c('AIE', 'AMIE')), 
      subset(pred_gKRLS_alt[-c(2,3)], QOI %in% c('AIE', 'AMIE'))
    )
    # Check the effects are correct on CI or local
    if (env_test != 'CRAN'){
      
      safe_cm <- function(x){if (is.matrix(x)){colMeans(x)}else{mean(x)}}
      AME_x1 <- mapply(levels(dat$x1), FUN=function(i){
        copy <- dat
        copy$x1 <- i
        return(safe_cm(predict(b, newdata = copy, type = 'response')))
      })
      AME_x0 <- sapply(c(0,1), FUN=function(i){
        copy <- dat
        copy$x0  <- i
        return(safe_cm(predict(b, newdata = copy, type = 'response')))
      })
      ACE <- sapply(levels(dat$x1), FUN=function(i){
        copy <- dat
        copy$x0 <- 0
        copy$x1 <- levels(dat$x1)[1]
        pred_0 <- safe_cm(predict(b, newdata = copy, type = 'response'))
        copy$x0 <- 1
        copy$x1 <- i
        pred_1 <- safe_cm(predict(b, newdata = copy, type = 'response'))
        return(pred_1 - pred_0)
      })
      if (!is.matrix(AME_x1)){
        ACE <- t(matrix(ACE))
        AME_x1 <- t(matrix(AME_x1))
        AME_x1 <- AME_x1 - AME_x1[1]
        AME_x0 <- diff(AME_x0)
      }else{
        AME_x1 <- sweep(AME_x1, MARGIN = 1, STATS = AME_x1[,1], FUN = '-')
        AME_x0 <- as.vector(diff(t(AME_x0)))
      }
      ACE <- ACE[,-1, drop = F]
      AME_x1 <- AME_x1[,-1, drop = F]
      AMIE <- ACE - AME_x1 - kronecker(matrix(AME_x0), t(matrix(rep(1,2))))
      
      AIE_low <- sapply(levels(dat$x1), FUN=function(i){
        copy <- dat
        copy$x1 <- levels(dat$x1)[1]
        copy$x0 <- 0
        pred_0 <- safe_cm(predict(b, newdata = copy, type = 'response'))

        copy <- dat
        copy$x1 <- i
        copy$x0 <- 0
        pred_1 <- safe_cm(predict(b, newdata = copy, type = 'response'))
        return(pred_1 - pred_0)
      })
      
      AIE_high <- sapply(levels(dat$x1), FUN=function(i){
        copy <- dat
        copy$x1 <- levels(dat$x1)[1]
        copy$x0 <- 1
        pred_0 <- safe_cm(predict(b, newdata = copy, type = 'response'))
        
        copy <- dat
        copy$x1 <- i
        copy$x0 <- 1
        pred_1 <- safe_cm(predict(b, newdata = copy, type = 'response'))
        return(pred_1 - pred_0)
      })
      
      if (is.matrix(AIE_high)){
        AIE <- (AIE_high - AIE_low)[,-1, drop = F]
      }else{
        AIE <- (AIE_high - AIE_low)[-1]
      }
      # Confirm all equal
      expect_equal(
        subset(pred_gKRLS, QOI == 'AME')$est,
        c(as.vector(t(cbind(AME_x0, AME_x1[,1,drop=F]))), AME_x1[,2]) 
      )
      expect_equal(
        subset(pred_gKRLS, QOI == 'ACE')$est,
        as.vector(ACE) 
      )
      expect_equal(
        subset(pred_gKRLS, QOI == 'AMIE')$est,
        as.vector(AMIE) 
      )
      expect_equal(
        subset(pred_gKRLS, QOI == 'AIE')$est,
        as.vector(AIE) 
      )
      
    }
  }
  
})

test_that("calculate interactions works with lpi (e.g., gaulss)", {
  
  n <- 400

  dat <- gamSim(eg = 3,  n = n, dist = 'gamma', verbose = FALSE)
  
  if (env_test == 'CRAN'){
    link_list <- list(gaulss())
  }else{
    link_list <- list(gaulss(), gevlss())
  }
  
  for (f in link_list){
    # print(f)
    fmla <- list(y ~ s(x1, x2, bs = 'gKRLS'), ~ s(x1, x2, bs = 'gKRLS'))
    if (f$family == 'gevlss'){
      fmla <- c(fmla, ~x1 + x2)
    }
    fit_mgcv <- gam(fmla, data = dat, family = f)
    predict_mgcv <- calculate_effects(fit_mgcv, continuous_type = 'predict')
    # Check predictions align
    expect_equal(predict_mgcv$est, 
      colMeans(predict(fit_mgcv, data = dat, type = 'response'))
    )
    # Shorten for CRAN by only checking interactions on CI or local
    if (env_test != 'CRAN'){
      predict_inter <- calculate_interactions(fit_mgcv, continuous_type = 'minmax',
                                              variables = list(c('x1', 'x2')))
      
      AME_x1 <- sapply(range(dat$x1), FUN=function(i){
        copy <- dat
        copy$x1 <- i
        return(colMeans(predict(fit_mgcv, newdata = copy, type = 'response')))
      })
      AME_x2 <- sapply(range(dat$x2), FUN=function(i){
        copy <- dat
        copy$x2 <- i
        return(colMeans(predict(fit_mgcv, newdata = copy, type = 'response')))
      })
      ACE <- mapply(range(dat$x1), range(dat$x2), FUN=function(i,j){
        copy <- dat
        copy$x1 <- i
        copy$x2 <- j
        return(colMeans(predict(fit_mgcv, newdata = copy, type = 'response')))
      })
      
      AME_x1 <- as.vector(diff(t(AME_x1)))
      AME_x2 <- as.vector(diff(t(AME_x2)))
      ACE <- as.vector(diff(t(ACE)))
      AMIE <- ACE - AME_x1 - AME_x2
      
      AIE_low <- sapply(range(dat$x2), FUN=function(i){
        copy <- dat
        copy$x1 <- min(dat$x1)
        copy$x2 <- i
        return(colMeans(predict(fit_mgcv, newdata = copy, type = 'response')))
      })
      AIE_high <- sapply(range(dat$x2), FUN=function(i){
        copy <- dat
        copy$x1 <- max(dat$x1)
        copy$x2 <- i
        return(colMeans(predict(fit_mgcv, newdata = copy, type = 'response')))
      })
      
      AIE <- as.vector(diff(t(AIE_high))) - as.vector(diff(t(AIE_low)))
      expect_equal(predict_inter$est, c(as.vector(t(cbind(AME_x1, AME_x2))), ACE, AMIE, AIE))
    }
    
  }

  
})
  