if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  env_test <- "CRAN"
  set.seed(129) # CRAN SEED
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

test_that("test raw works and fails where expected", {
  
  dat <- gamSim(eg = 1, n =100)
  dat$x0 <- dat$x0 > 0
  dat$x1 <- cut(dat$x1, 3)
  b <- gam(y ~ x0 + x1 + s(x2, by =x1), data = dat)
  
  est_effects <- calculate_effects(b, raw = TRUE)
  # Ordered by effect, raw_1, raw_0
  expect_true(all(sapply(split(est_effects$est, est_effects$variable), FUN=function(i){
    i[1] - (i[2] - i[3])
  }) == 0))
  
  est_effects <- calculate_effects(b, continuous_type = 'derivative', 
       conditional = data.frame(x1 = unique(dat$x1)), 
       variables = 'x2', raw = TRUE)
  
  expect_equivalent(sapply(split(est_effects$est, est_effects$x1), FUN=function(i){
    i[1] - (i[2] - i[3]) / (2 * 1e-7 * max(1, max(abs(dat$x2))))
  }), rep(0, 3))
  
  est_effects <- calculate_effects(b, continuous_type = 'second_derivative', 
       variables = 'x2', raw = TRUE, individual = T)
  # Check formula for second derivative lines up  
  expect_equivalent(
    est_effects$est[1],
    sum(est_effects$est[-1] * c(1, -2, 1))/(1e-7 * max(1, max(abs(dat$x2))))
  )
  individual_raw <- get_individual_effects(est_effects)
  expect_equivalent(sapply(split(individual_raw$est, individual_raw$obs), FUN=function(i){
    i[1] - sum(i[-1] * c(1, -2,1)) / (1e-7 * max(1, max(abs(dat$x2))))
  }), rep(0, length(unique(individual_raw$obs))))

  # Test for errors  
  expect_error(calculate_effects(b, continuous_type = 'predict', raw=T), regexp = 'raw must be FALSE')
  expect_error(calculate_interactions(b, raw = T), regexp = 'raw=T not permitted')
})