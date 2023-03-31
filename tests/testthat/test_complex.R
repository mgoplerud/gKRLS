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
  dat$new_y <- dat$y - 1
  
  # From mgcv: Confirm predictions line up for
  # single linear predictor models

  for (f in list(ocat(R=R), poisson(link = 'identity'), gaussian(), tw(), nb(), scat(), ziP())){
    if (!grepl(f$family, pattern='^Zero')){
      b <- suppressWarnings(gam(y~ s(x0) + s(x1) + s(x2) + s(x3),family=f,data=dat))
    }else{
      b <- suppressWarnings(gam(new_y~ s(x0) + s(x1) + s(x2) + s(x3),family=f,data=dat))
    }
    pred_gKRLS <- calculate_effects(model = b, continuous_type = 'predict', individual = T)
    pred_mgcv <- predict(b, newdata = dat, se.fit = TRUE, type = 'response')    
    expect_equivalent(pred_gKRLS$individual$est, as.vector(pred_mgcv$fit))
    expect_equivalent(pred_gKRLS$individual$se, as.vector(pred_mgcv$se.fit))
    rm(pred_gKRLS, pred_mgcv)
  }
  
  # Multinomial Test: Note that this will *not* work with 1.8-42 or below when K >= 3
  b <- gam(list(new_y ~ s(x0), 
                ~ x0 + x1 + x2, 
                ~ s(x3, x2, bs = 'gKRLS')), 
           data = dat,
           family = multinom(K = 3))

  pred_gKRLS <- calculate_effects(model = b, continuous_type = 'predict', individual = T)
  pred_mgcv <- predict(b, newdata = dat, se.fit = TRUE, type = 'response')

  expect_equivalent(pred_gKRLS$individual$est, as.vector(pred_mgcv$fit))
  expect_true(!isTRUE(all.equal(pred_gKRLS$individual$se, as.vector(pred_mgcv$se.fit))))

  # Test when K = 2 (i.e. 3 categories)
  new_dat <- subset(dat, new_y < 3)
  b <- gam(list(new_y ~ s(x0), 
                ~ s(x3, x2, bs = 'gKRLS')), 
           data = new_dat,
           family = multinom(K = 2))
  
  pred_gKRLS <- calculate_effects(model = b, continuous_type = 'predict', individual = T)
  pred_mgcv <- predict(b, newdata = new_dat, se.fit = TRUE, type = 'response')
  
  expect_equivalent(pred_gKRLS$individual$est, as.vector(pred_mgcv$fit))
  expect_equivalent(pred_gKRLS$individual$se, as.vector(pred_mgcv$se.fit))
  
})
