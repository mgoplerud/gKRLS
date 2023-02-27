context("tests for general/extended families in mgcv")

test_that(" Test for prediction/SE for complex families ", {
  
  # From mgcv
  n <- 400
  dat <- gamSim(1,n=n)
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
  
  # From mgcv: Confirm predictions line up for
  # single linear predictor models

  for (f in list(ocat(R=R), poisson(link = 'identity'), tw(), nb(), scat(), ziP())){
    print(f)
    b <- suppressWarnings(gam(y~ s(x0) + s(x1) + s(x2) + s(x3),family=f,data=dat))
    
    pred_gKRLS <- calculate_effects(model = b, continuous_type = 'predict', individual = T)
    pred_mgcv <- predict(b, newdata = dat, se.fit = TRUE, type = 'response')    
    expect_equivalent(pred_gKRLS$individual$est, as.vector(pred_mgcv$fit))
    expect_equivalent(pred_gKRLS$individual$se, as.vector(pred_mgcv$se.fit))
    rm(pred_gKRLS, pred_mgcv)
  }
  
  # Still need to set up for multinomial
  # dat$new_y <- dat$y - 1
  # b <- gam(list(new_y ~ s(x0), ~ x0 + x1 + x2, ~ s(x3, x2, bs = 'gKRLS')), data = dat,
  #          family = multinom(K = 3))
  # 
  # pred_gKRLS <- calculate_effects(model = b, continuous_type = 'predict', individual = T)
  # pred_mgcv <- predict(b, newdata = dat, se.fit = TRUE, type = 'response')
  # 
  # expect_equivalent(pred_gKRLS$individual$est, as.vector(pred_mgcv$fit))
  # expect_equivalent(pred_gKRLS$individual$se, as.vector(pred_mgcv$se.fit))
})
