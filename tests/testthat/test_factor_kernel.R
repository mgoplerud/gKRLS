test_that("Test everything runs when kernel has categorical/factor variables", {
  
  set.seed(4161)
  
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
    legacy_marginal_effect(fit_gam, newdata = df), 
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
  expect_error(legacy_marginal_effect(fit_gam2, newdata = df[1:5,]),
               regexp = 'factors are provided')
})

test_that("Test that factor vs dummies is equivalent", {
  
<<<<<<< HEAD
  set.seed(4999)
=======

>>>>>>> 48cd6b5a4963fc0e14db0e248ee6ae2707a59ebe
  N <- 50
  X <- cbind(matrix(rnorm(N * 2), ncol = 2), rbinom(N, 1, 0.5))
  y <- X %*% rnorm(ncol(X))
  y <- rpois(length(y), exp(y))
  X4 <- sample(letters[1:5], N, replace = T)
  
  df <- data.frame(X, X4 = factor(X4), y, stringsAsFactors = F)
  
  set.seed(54321)
  
  fit_factor <- gam(
    y ~ s(X1, X4, bs = 'gKRLS'), data = df
  )
  wide_X4 <- model.matrix(~ 0 + X4)
  df <- cbind(df, wide_X4)
  
  set.seed(54321)
  
  fmla <- as.formula(paste0('y ~ s(X1,', paste0(colnames(wide_X4), collapse=','), ', bs = "gKRLS")'))
  fit_direct <-  gam(
    fmla, data = df
  )
  expect_equivalent(coef(fit_factor), coef(fit_direct), tolerance = 1e-6)
  expect_equivalent(fitted(fit_factor), fitted(fit_direct), tolerance = 1e-6)
  
  mfx_direct <- calculate_effects(fit_direct, variables = 'X1', continuous_type = 'derivative')
  mfx_factor <- calculate_effects(fit_factor, variables = 'X1', continuous_type = 'derivative')
  expect_equivalent(mfx_direct, mfx_factor, tol = 1e-6, scale = 1)
  legacy_direct <- legacy_marginal_effect(fit_direct, newdata = df, keep = 'X1')
  expect_equivalent(
    legacy_direct$AME_pointwise,
    mfx_direct$marginal_effects$est,
    tol = 0.01, scale = 1
  )
  expect_equivalent(
    sqrt(legacy_direct$AME_pointwise_var),
    mfx_direct$marginal_effects$se,
    tol = 0.01, scale = 1
  )
})
