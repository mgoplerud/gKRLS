
context("Marginal Effects Tests")

test_that("mfx with degenerate Mahalanobis", {
  
  N <- 100
  x1 <- rnorm(N)
  x2 <- rbinom(N,size=1,prob=.2)
  y <- x1^3 - 0.5 *x2 + rnorm(N,0, 1)
  y <- y * 10
  X <- cbind(x1, x2, x1 + x2 * 3)
  colnames(X) <- paste0('x', 1:ncol(X))
  fit_gKRLS <- gam(y ~ 0 + s(x1,x2,x3, bs = 'gKRLS', 
                            xt = gKRLS(force_base = TRUE, standardize = 'Mahalanobis')),
                   family = gaussian(), method = 'REML', data = data.frame(y, X)
  )
  expect_equivalent(fit_gKRLS$smooth[[1]]$term, c('x1', 'x2', 'x3'))
  expect_equivalent(ncol(fit_gKRLS$smooth[[1]]$X_train), 2)
  
  mfx_num <- calculate_effects(fit_gKRLS, data = data.frame(X), continuous_type = 'deriv')
  mfx_legacy <- legacy_marginal_effect(fit_gKRLS, newdata = data.frame(X),
                         keep = c('x1', 'x2', 'x3'))
  
  expect_equivalent(mfx_num$marginal_effects$est, mfx_legacy$AME_pointwise, tol = 1e-3)
  expect_equivalent(mfx_num$marginal_effects$se, sqrt(mfx_legacy$AME_pointwise_var), tol = 1e-3)
  
})

test_that('Test MFX', {

  N <- 100
  x1 <- rnorm(N)
  x2 <- rbinom(N,size=1,prob=.2)
  y <- x1^3 - 0.5 *x2 + rnorm(N,0, 1)
  y <- y * 10
  X <- cbind(x1, x2, x1 + x2 * 3)
  colnames(X) <- paste0('x', 1:ncol(X))
  fit_gKRLS <- gam(y ~ 0 + s(x1,x2,x3, bs = 'gKRLS', 
                             xt = gKRLS(force_base = TRUE, standardize = 'scaled', sketch_method = 'gaussian')),
                   family = gaussian(), method = 'REML', data = data.frame(y, X)
  )
  expect_equivalent(fit_gKRLS$smooth[[1]]$term, c('x1', 'x2', 'x3'))
  expect_equivalent(ncol(fit_gKRLS$smooth[[1]]$X_train), 3)
  
  mfx_num <- calculate_effects(fit_gKRLS, data = data.frame(X), continuous_type = 'deriv')
  mfx_legacy <- legacy_marginal_effect(fit_gKRLS, newdata = data.frame(X),
                                       keep = c('x1', 'x2', 'x3'))
  
  expect_equivalent(mfx_num$marginal_effects$est, mfx_legacy$AME_pointwise, tol = 1e-3)
  expect_equivalent(mfx_num$marginal_effects$se, sqrt(mfx_legacy$AME_pointwise_var), tol = 1e-3)
  
})

test_that("test 'calculate_effects'", {

  N <- 100
  x1 <- rnorm(N)
  x2 <- rbinom(N,size=1,prob=.2)
  y <- x1^3 - 0.5 *x2 + rnorm(N,0, 1)
  y <- y * 10
  s <- sample(letters[1:5], N, replace = T)
  X <- data.frame(x1, x2, s = factor(s), stringsAsFactors = F)
  
  colnames(X) <- paste0('x', 1:ncol(X))
  fit_gKRLS <- gam(y ~ s(x1,x2,x3, bs = 'gKRLS'),
   family = gaussian(), method = 'REML', data = data.frame(y, X)
  )
  
  factor_test <- calculate_effects(model = fit_gKRLS, variables = 'x3')
  expect_null(factor_test$individual)
  factor_test <- calculate_effects(model = fit_gKRLS, variables = 'x3',
      individual = TRUE, conditional = data.frame(x1 = c(-1, 1)))
  expect_false(is.null(factor_test$individual))
  expect_equal(nrow(factor_test$marginal_effects), 8)
  
  fit_inter <- kernel_interactions(
    fit_gKRLS, variables = list(c('x1', 'x3')))
  # TO DO: tests for this
})
