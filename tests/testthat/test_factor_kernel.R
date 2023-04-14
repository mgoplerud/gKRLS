if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  env_test <- "CRAN"
  set.seed(128) # CRAN SEED
}else{
  # If on local machine
  env_test <- 'local'
}


test_that("Test everything runs when kernel has categorical/factor variables", {
  
  N <- 50
  X <- cbind(matrix(rnorm(N * 2), ncol = 2), rbinom(N, 1, 0.5))
  y <- X %*% rnorm(ncol(X))
  y <- rpois(length(y), exp(y))
  X4 <- sample(letters[1:5], N, replace = T)
  X5 <- factor(sample(state.region, N, replace = T))
  X6 <- factor(sample(state.abb[1:5], N, replace = T), ordered = TRUE)
  X7 <- sample(as.integer(1:15), N, replace = T)
  
  df <- data.frame(X, X4 = factor(X4), X5, X6, X7, y, stringsAsFactors = F)
  fit_gam <- gam(y ~ 0 +
                   s(X1, X2, X3, X4, X5, X6, X7, bs = "gKRLS",
                     xt = gKRLS(standardize = 'none', sketch_method = 'none')), data = df
  )
  
  contrast_mm <- lapply(df[,c('X4', 'X5', 'X6')],
                        FUN=function(i){contrasts(i, contrasts = FALSE)})
  
  mm <- model.matrix(~ 0 + X1 + X2 + X3 + X4 + X5 + X6 + X7, 
                     data = df,
                     contrasts.arg = contrast_mm)
  # Check the data aligns (no dropped levels), factors expanded correctly
  expect_equivalent(mm, fit_gam$smooth[[1]]$X_train, tolerance = 1e-6)
  
  
  v1 <- predict(fit_gam, 
                newdata = data.frame(X1 = -5:5, X2 = 0:10,
                                     X3 = 2, X4 = 'a', X5 = 'South', X6 = 'CA',
                                     X7 = 1.5, stringsAsFactors = F))
  v2 <- predict(fit_gam, 
                newdata = data.frame(X1 = -5:5, X2 = 0:10,
                                     X3 = 2, X4 = 'a', X5 = 'South', X6 = 'CA',
                                     X7 = 1.5, stringsAsFactors = T))
  expect_equal(v1, v2, tolerance = 1e-6)
  mfx_calc <- calculate_effects(fit_gam)
  expect_error(
    suppressMessages(legacy_marginal_effect(fit_gam)), 
    regexp = 'factors are provided')
  expect_s3_class(mfx_calc, 'gKRLS_mfx')
  
  expect_warning(predict(fit_gam, 
    newdata = data.frame(X1 = -5:5, X2 = 0:10,
      X3 = 2, X4 = 'a', X5 = 'NEW', X6 = 'CA',
      X7 = 1.5, stringsAsFactors = T)), regexp = 'factor levels NEW')
  
  v1 <- suppressMessages(suppressWarnings(predict(fit_gam, 
                                                  newdata = data.frame(X1 = -5:5, X2 = 0:10,
                                                                       X3 = 2, X4 = 'a', X5 = 'NEW', X6 = 'CA',
                                                                       X7 = 1.5, stringsAsFactors = F))))
  v2 <- suppressMessages(suppressWarnings(predict(fit_gam, 
                                                  newdata = data.frame(X1 = -5:5, X2 = 0:10,
                                                                       X3 = 2, X4 = 'a', X5 = 'NEW', X6 = 'CA',
                                                                       X7 = 1.5, stringsAsFactors = T))))
  expect_equal(v1, v2, tolerance = 1e-6)
  
  fit_gam2 <- gam(y ~ 1 +
                    s(X1, X2, X3, X4, X5, X6, X7, bs = "gKRLS"), 
                  data = df
  )
  predict(fit_gam2, newdata = df[1:5,])
  expect_s3_class(calculate_effects(fit_gam2, variables = 'X5'), 'gKRLS_mfx')
  expect_error(legacy_marginal_effect(fit_gam2, data = df[1:5,]),
               regexp = 'factors are provided')
})

test_that("Test that factor vs dummies is equivalent", {
  
  N <- 50
  X <- cbind(matrix(rnorm(N * 2), ncol = 2), rbinom(N, 1, 0.5))
  y <- X %*% rnorm(ncol(X))
  y <- rpois(length(y), exp(y))
  X4 <- sample(letters[1:5], N, replace = T)
  
  df <- data.frame(X, X4 = factor(X4), y, stringsAsFactors = F)
  
  id_nystrom <- sample(1:N, ceiling(5 * N^(1/3)))

  fit_factor <- gam(
    y ~ s(X1, X4, bs = 'gKRLS', xt = gKRLS(sketch_method = id_nystrom)), data = df
  )
  wide_X4 <- model.matrix(~ 0 + X4)
  df <- cbind(df, wide_X4)

  fmla <- as.formula(paste0('y ~ s(X1,', 
        paste0(colnames(wide_X4), collapse=','), 
        ', bs = "gKRLS", xt = gKRLS(sketch_method = id_nystrom))'))
  
  fit_direct <-  gam(
    fmla, data = df
  )
  expect_equivalent(coef(fit_factor), coef(fit_direct), tolerance = 1e-6)
  expect_equivalent(fitted(fit_factor), fitted(fit_direct), tolerance = 1e-6)
  
  mfx_direct <- calculate_effects(fit_direct, variables = 'X1', continuous_type = 'derivative')
  mfx_factor <- calculate_effects(fit_factor, variables = 'X1', continuous_type = 'derivative')
  expect_equivalent(mfx_direct, mfx_factor, tol = 1e-6)
  legacy_direct <- legacy_marginal_effect(fit_direct, variables = 'X1')
  expect_equivalent(
    legacy_direct$AME_pointwise,
    mfx_direct$est,
    tol = 0.01, scale = 1
  )
  expect_equivalent(
    sqrt(legacy_direct$AME_pointwise_var),
    mfx_direct$se,
    tol = 0.01, scale = 1
  )
  
  # Test that "summary" works as expected  
  expect_true(nrow(summary(legacy_direct)) == 1)
  expect_equal(summary(mfx_direct), mfx_direct)
  
  test_custom <- legacy_marginal_effect(fit_direct, data = df[1:5,], variables = 'X1')
  expect_equal(nrow(test_custom$ME_pointwise), 5)
  expect_equal(legacy_direct$ME_pointwise_var[1:5,,drop=F], test_custom$ME_pointwise_var)
})

test_that("test that changing factor reference category works in sensible ways", {
  
  dat <- gamSim(eg = 4, n = 100, verbose = FALSE)
  dat$x0 <- dat$x0 > median(dat$x0)
  dat$fac <- factor(letters[dat$fac])
  dat$y <- rnorm(3)[match(dat$fac, letters)] + dat$y
  
  id <- sample(1:nrow(dat), 5 * ceiling(100^(1/3)))
  b <- gam(y ~ x0 + s(x2, fac, bs = 'gKRLS', xt = gKRLS(sketch_method = id)), data = dat)
  
  dat$fac_new <- relevel(dat$fac, ref = 'c')
  b_refit <- gam(y ~ x0 + s(x2, fac_new, bs = 'gKRLS', xt = gKRLS(sketch_method = id)), data = dat)
  
  fit_default <- calculate_effects(b)
  fit_refit <- calculate_effects(b_refit)
  expect_equivalent(fit_default[1:2,], fit_refit[1:2,], tol = 1e-6)
  
  # Confirm that changing the factor reorders the data.frame
  expect_equal(
    fit_default[fit_default$variable == 'factor(fac)b',]$est - fit_default[fit_default$variable == 'factor(fac)c',]$est,
    fit_refit[fit_refit$variable == 'factor(fac_new)b',]$est,
    tol = 1e-6
  )
  expect_equal(
    -fit_default[fit_default$variable == 'factor(fac)c',]$est ,
    fit_refit[fit_refit$variable == 'factor(fac_new)a',]$est,
    tol = 1e-6
  )
  
  # Confirm that re-leveling the factor in the custom data and using
  # the originally estimated model returns the expected results
  new_dat <- dat
  new_dat$fac <- factor(new_dat$fac, levels = c('c', 'a', 'b'))
  fit_custom <- expect_warning(calculate_effects(b, data = new_dat), regexp  ='For custom argument')
  fit_override <- calculate_effects(b, data = new_dat, use_original = TRUE)
  
  expect_equivalent(fit_override, fit_default, tol = 1e-6)
  expect_equivalent(fit_refit[,-1], fit_custom[,-1], tol = 1e-6)
  
})