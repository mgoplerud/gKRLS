context("test non-standard standard errors")

test_that("Test Robust for Linear (Unpenalized)", {
  
  N <- 100
  x <- rnorm(N)
  y <- rnorm(N, sin(x)) * 100
  est_linear_gam <- gam(y ~ x)
  est_lm <- lm(y ~ x)
  # Check fitted models are equivalent
  expect_equivalent(sigma(est_lm), sigma(est_linear_gam))
  expect_equivalent(fitted(est_lm), fitted(est_linear_gam))
  expect_equivalent(vcov(est_lm), vcov(est_linear_gam))
  
  hc_lm <- vcovHC(est_lm, type = 'HC0')
  hc_linear_gam <- vcovHC(est_linear_gam, type = 'HC0')
  
  ef_gam <- estfun(est_linear_gam, override_check = TRUE)
  X_gam <- predict(est_linear_gam, newdata = data.frame(x = x), type = 'lpmatrix')
  
  expect_equivalent(ef_gam, as.matrix(Diagonal(x = residuals(est_lm)) %*% X_gam))
  expect_equivalent(ef_gam, estfun(est_lm))
  tXXinv <- solve(crossprod(X_gam))
  direct_estimation <- tXXinv %*%
    t(X_gam) %*% Diagonal(x = residuals(est_linear_gam)^2) %*% X_gam %*% 
    tXXinv
  expect_equivalent(as.matrix(direct_estimation), hc_lm)
  expect_equivalent(hc_linear_gam, hc_lm)
  expect_equivalent(vcovHC(est_lm, type = 'HC1'),
                    vcovHC(est_linear_gam, type = 'HC1'))

  est_bam_linear <- bam(y ~ x)
  expect_equivalent(vcov(est_bam_linear), vcov(est_lm))
  expect_equivalent(vcovHC(est_bam_linear, type = 'HC1'), 
    vcovHC(est_lm, type = 'HC1'))
  expect_equivalent(colSums(ef_gam), rep(0, ncol(ef_gam)))
})

test_that("Test Robust for Linear (Penalized)", {
  
  N <- 1000
  x <- rnorm(N)
  z <- rnorm(N)
  y <- rnorm(N, sin(x * z)) * 2
  est_gam <- gam(y ~ te(x,z))
  hc_gam <- vcovHC(est_gam, type = 'HC0')

  ef_gam <- estfun(est_gam, override_check = TRUE)
  X_gam <- predict(est_gam, newdata = data.frame(x = x, z = z), type = 'lpmatrix')
  
  expect_equivalent(ef_gam, as.matrix(Diagonal(x = residuals(est_gam)) %*% X_gam))
  expect_equivalent(ef_gam, estfun(est_gam, override_check = T))
  S <- est_gam$smooth[[1]]$S[[2]] * est_gam$sp[2] + est_gam$smooth[[1]]$S[[1]] * est_gam$sp[1]
  S <- bdiag(0, S)
  m <- solve(crossprod(X_gam) + S, t(X_gam))
  direct_estimation <- m %*% Diagonal(x = residuals(est_gam)^2) %*% t(m)
  # Check direct formula for HC0
  expect_equivalent(as.matrix(direct_estimation), hc_gam, tol = 1e-5)
  # Check that HC1 works for "naive" correction
  manual_hc1 <- vcovHC(est_gam, type = 'HC0') * N/(N-sum(est_gam$edf))
  expect_equivalent(manual_hc1,
    vcovHC(est_gam, type = 'HC1'))
  
  expect_equivalent(as.vector(colSums(ef_gam) - S %*% coef(est_gam)), 
    rep(0, ncol(S)), tol = 1e-5)
  
})

test_that("Test Robust for GLM", {
  
  N <- 100
  x <- rnorm(N)
  z <- rnorm(N)
  
  # Check for poisson
  y <- round(exp(rnorm(N, sin(x))))
  est_gam <- gam(y ~ s(x,z, bs = 'gKRLS'), family = poisson())
  hc_gam <- vcovHC(est_gam, type = 'HC0')
  
  ef_gam <- estfun(est_gam, override_check = TRUE)
  X_gam <- predict(est_gam, newdata = data.frame(x = x, z = z), type = 'lpmatrix')
  expect_equivalent(ef_gam, as.matrix(Diagonal(x = est_gam$y - fitted(est_gam)) %*% X_gam))
  expect_equivalent(ef_gam, estfun(est_gam, override_check = T))
  S <- est_gam$smooth[[1]]$S[[1]] * est_gam$sp[1] 
  S <- bdiag(0, S)
  tXXinv <- solve(t(X_gam) %*% Diagonal(x = fitted(est_gam)) %*% X_gam + S)
  direct_estimation <- tXXinv %*%
    t(X_gam) %*% Diagonal(x = (y - fitted(est_gam))^2 ) %*% X_gam %*% 
    tXXinv
  # Check direct formula for HC0
  expect_equivalent(as.matrix(direct_estimation), hc_gam)
  # Check that HC1 works for "naive" correction
  manual_hc1 <- vcovHC(est_gam, type = 'HC0') * N/(N-sum(est_gam$edf))
  expect_equivalent(manual_hc1,
                    vcovHC(est_gam, type = 'HC1'))
  
  # Check for binomial with NONSTANDARD LINK
  N <- 1000
  x <- rnorm(N)
  z <- rnorm(N)
  y <- rbinom(N, 1, plogis( exp(x) + cos(z)))
  
  est_gam <- gam(y ~ s(x) + s(z), family = binomial(link = 'cloglog'))
  hc_gam <- vcovHC(est_gam, type = 'HC0')
  
  ef_gam <- estfun(est_gam, override_check = TRUE)
  X_gam <- predict(est_gam, newdata = data.frame(x = x, z = z), type = 'lpmatrix')
  lp_gam <- X_gam %*% coef(est_gam)
  score_weight <- y * exp(-exp(lp_gam) + lp_gam)/(1 - exp(-exp(lp_gam))) +
    (1-y) * -exp(lp_gam)
  hessian_weight <- exp(lp_gam) * (-1 - exp(exp(lp_gam)) * (-1 + exp(lp_gam)))/(-1 + exp(exp(lp_gam)))^2
  e_y <- 1 - exp(-exp(lp_gam))
  if (any(e_y %in% c(0,1))){
    hessian_weight[which(e_y %in% c(0,1))] <- 0
  }
  hessian_weight <- as.vector(e_y * hessian_weight + (1-e_y) * -exp(lp_gam))
  
  expect_equivalent(
    score_weight, residuals(est_gam, 'pearson') * sqrt(weights(est_gam, 'working')),
    tol = 1e-5
  )
  
  expect_equivalent(ef_gam, as.matrix(Diagonal(x = score_weight) %*% X_gam))
  expect_equivalent(ef_gam, estfun(est_gam, override_check = T))
  if (length(est_gam$smooth) == 0){
    S <- Diagonal(x = rep(0, ncol(ef_gam)))
  }else{
    S <- bdiag(est_gam$smooth[[1]]$S[[1]] * est_gam$sp[1] , est_gam$smooth[[2]]$S[[1]] * est_gam$sp[2] )
    S <- bdiag(0, S)
  }
  
  tXXinv <- solve(t(X_gam) %*% Diagonal(x = -hessian_weight) %*% X_gam + S)
  expect_equivalent(as.matrix(tXXinv), vcov(est_gam))
  
  direct_estimation <- tXXinv %*%
    t(X_gam) %*% Diagonal(x = score_weight^2) %*% X_gam %*% 
    tXXinv
  # Check direct formula for HC0
  expect_equivalent(as.matrix(direct_estimation), hc_gam)
  # Check that HC1 works for "naive" correction
  manual_hc1 <- vcovHC(est_gam, type = 'HC0') * N/(N-sum(est_gam$edf))
  expect_equivalent(manual_hc1, vcovHC(est_gam, type = 'HC1'))
  expect_equivalent(as.vector(colSums(ef_gam) - S %*% coef(est_gam)), rep(0, ncol(S)),
                    tol = 1e-6)
})

test_that("Test Robust for bam", {
  
  N <- 1000
  x <- rnorm(N)
  z <- rnorm(N)
  y <- rbinom(N, 1, plogis( exp(x) + cos(z)))
  
  est_gam <- gam(y ~ x + z, family = binomial(link = 'cloglog'), method = 'REML')
  est_bam <- bam(y ~ x + z, family = binomial(link = 'cloglog'), method = 'REML')
  
  expect_equivalent(
    vcovHC(est_gam, type = 'HC1'), 
    vcovHC(est_gam, type = 'HC0'),
    tol = 1e-5)

  est_gam <- gam(y ~ x + z, family = gaussian())
  est_bam <- bam(y ~ x + z, family = gaussian())
  c1 <- sample(1:15, size = N, replace = T)
  expect_equivalent(
    vcovCL(est_gam, cluster = c1), 
    vcovCL(est_gam, cluster = c1),
    tol = 1e-5)

  est_bam <- bam(y ~ te(x,z), family = gaussian())
  mm <- predict(est_bam, newdata = data.frame(x, z), type = 'lpmatrix')
  S <- est_bam$smooth[[1]]$S[[1]] * est_bam$sp[1] + est_bam$smooth[[1]]$S[[2]] * est_bam$sp[2]
  meat <- t(mm) %*% Diagonal(x = residuals(est_bam)^2) %*% mm
  inv <- solve(crossprod(mm) + bdiag(0, S))
  expect_equivalent(as.matrix(inv %*% meat %*% inv), 
    vcovHC(est_bam, type = 'HC0'), tol = 1e-5)
})

test_that("Test Robust for Complex Family", {
  
  N <- 100
  x <- rnorm(N)
  z <- rnorm(N)
  y <- exp(z) + cos(x) + rnorm(N)
  
  est_gam <- gam(y ~ x + z, family = scat(), method = 'REML')
  expect_error(vcovHC(est_gam), regexp = 'Robust SE')  
  
})
