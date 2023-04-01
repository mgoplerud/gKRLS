if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  env_test <- "CRAN"
  set.seed(129)
}else{
  # If on local machine
  env_test <- 'local'
}


test_that("does not run with prior weights", {
  
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
  
  fit <- gam(y ~ s(x1), weights = x2, data = dat)  
  expect_warning(calculate_effects(fit), regexp = 'ignores "weights"')
  fit <- gam(y ~ s(x1), weights = x2, data = dat, family = ocat(R=R))  
  expect_warning(calculate_effects(fit), regexp = 'ignores "weights"')
  dat$new_y <- dat$y - 1
  fit <- gam(list(new_y ~ s(x1), ~ x1, ~ s(x2)), weights = x2, data = dat, family = multinom(K=3))  
  expect_warning(calculate_effects(fit), regexp = 'ignores "weights"')
  
})

test_that("does not run with offset", {
  
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
  
  fit <- gam(y ~ offset(x2) + s(x1), data = dat)  
  expect_error(calculate_effects(fit), regexp = '"offset" not set')
  fit <- gam(y ~ offset(x2) + s(x1), data = dat, family = ocat(R=R))  
  expect_error(calculate_effects(fit), regexp = '"offset" not set')
  dat$new_y <- dat$y - 1
  fit <- gam(list(new_y ~ offset(x2) + s(x1), ~ x1, ~ s(x2)), data = dat, family = multinom(K=3))  
  expect_error(calculate_effects(fit), regexp = '"offset" not set')
  fit <- gam(list(new_y ~ offset(x2) + s(x1), ~ offset(x2) + x1, ~ s(x2)), data = dat, family = multinom(K=3))  
  expect_error(calculate_effects(fit), regexp = '"offset" not set')
  
})
