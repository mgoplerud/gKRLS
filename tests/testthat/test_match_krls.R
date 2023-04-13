if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  env_test <- "CRAN"
  set.seed(130) # CRAN SEED
}else{
  # If on local machine
  env_test <- 'local'
}


test_that("gKRLS agrees with direct solution", {
  
  N <- 200
  x1 <- rnorm(N)
  x2 <- rbinom(N, size = 1, prob = .2)
  y <- x1^3 - 0.5 * x2 + rnorm(N, 0, 1)
  y <- y * 10
  X <- cbind(x1, x2)
  X <- cbind(X, model.matrix(~ 0 + sample(letters[1:5], N, replace = T)))
  colnames(X) <- paste0("x", 1:ncol(X))
  X_copy <- cbind(X, X, X)
  colnames(X_copy) <- paste0("x", 1:ncol(X_copy))

  est_gKRLS_single <- gam(y ~ s(x1, bs = "gKRLS"),
    family = gaussian(), data = data.frame(y, X)
  )

  est_gKRLS <- gam(y ~ s(x1, x2, x3, x4, x5, x6, x7,
    bs = "gKRLS",
    xt = gKRLS(sketch_method = "none")
  ),
  family = gaussian(), data = data.frame(y, X)
  )

  est_copy <- gam(y ~ s(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13,
    x14, x15, x16, x17, x18, x19, x20, x21,
    bs = "gKRLS",
    xt = gKRLS(sketch_method = "none")
  ),
  family = gaussian(), data = data.frame(y, X_copy)
  )

  expect_equal(fitted(est_gKRLS), fitted(est_copy), tol = 1e-5)

  fit_gKRLS <- gam(y ~ 0 + s(x1, x2, x3, x4, x5, x6, x7,
    bs = "gKRLS",
    xt = gKRLS(
      standardize = "scaled", 
      sketch_method = "none",
    )
  ),
  family = gaussian(), method = "REML", data = data.frame(y, X)
  )

  K_manual <- exp(-as.matrix(dist(apply(X, MARGIN = 2, scale))^2) / fit_gKRLS$smooth[[1]]$bandwidth)
  pred_kern <- predict(fit_gKRLS, newdata = data.frame(X), type = "lpmatrix")

  # if (requireNamespace('KRLS', quietly = TRUE)){
  #
  #   fit_krls <- suppressWarnings(KRLS::krls(X = X, y = y))
  #
  #   direct_fitted_KRLS <- as.vector(K_manual %*%
  #    solve(K_manual + as.numeric(fit_krls$lambda) * Diagonal(n = nrow(K_manual))) %*%
  #    scale(y)) * sd(y) + mean(y)
  #
  #   expect_equivalent(as.vector(fitted(fit_krls)), as.vector(direct_fitted_KRLS))
  #   rm(direct_fitted_KRLS, fit_krls)
  # }
  #
  direct_fit <- as.vector(pred_kern %*%
    solve(crossprod(pred_kern) + fit_gKRLS$sp * fit_gKRLS$smooth[[1]]$S[[1]], t(pred_kern) %*% y))
  expect_equivalent(direct_fit, fitted(fit_gKRLS), tol = 1e-5)
})


test_that("Legacy Agrees with Numerical", {
  
  N <- 200
  x1 <- rnorm(N)
  x2 <- rbinom(N, size = 1, prob = .2)
  y <- x1^3 - 0.5 * x2 + rnorm(N, 0, 1)
  y <- y * 10
  X <- cbind(x1, x2)
  X <- cbind(X, model.matrix(~ 0 + sample(letters[1:5], N, replace = T)))
  colnames(X) <- paste0("x", 1:ncol(X))
  fit_gKRLS <- gam(y ~ 0 + s(x1, x2, x3, x4, x5, x6, x7, bs = "gKRLS", xt = gKRLS(standardize = "scaled")),
    family = gaussian(), method = "REML", data = data.frame(y, X)
  )

  mfx_gKRLS <- legacy_marginal_effect(fit_gKRLS, data = data.frame(X))
  mfx_numerical <- calculate_effects(
    model = fit_gKRLS, individual = TRUE,
    data = data.frame(X), continuous_type = "deriv"
  )

  expect_equivalent(mfx_gKRLS$AME_pointwise, mfx_numerical$est, tol = 1e-2)
  expect_equivalent(sqrt(mfx_gKRLS$AME_pointwise_var), mfx_numerical$se, tol = 1e-2)

  ind_est <- get_individual_effects(mfx_numerical)
  expect_equivalent(
    do.call("rbind", split(ind_est$est, ind_est$obs)),
    mfx_gKRLS$ME_pointwise,
    tol = 1e-4
  )

  expect_equivalent(
    do.call("rbind", split(ind_est$se^2, ind_est$obs)),
    mfx_gKRLS$ME_pointwise_var,
    tol = 1e-4
  )
})


test_that("Logistic KRLS Tests", {
  
  N <- 200
  x1 <- rnorm(N)
  x2 <- rbinom(N, size = 1, prob = .2)
  x2 <- x2 + sample(c(1, 0, -1), size = N, replace = T) * 1e-6
  b1 <- 1
  b2 <- -3
  y <- b1 * x1^3 + b2 * x2 + rnorm(N, 0, .15)
  X <- cbind(x1, x2)

  bin_y <- rbinom(nrow(X), 1, plogis(scale(y)))

  fit_binary_gKRLS <- gam(bin_y ~ s(x1, x2,
    xt = gKRLS(sketch_method = "gaussian"),
    bs = "gKRLS"
  ),
  data = data.frame(X, bin_y),
  family = binomial()
  )

  mfx_logit <- legacy_marginal_effect(fit_binary_gKRLS, data = data.frame(X))
  mfx_logit_num <- calculate_effects(fit_binary_gKRLS,
    data = data.frame(X),
    continuous_type = "deriv", individual = TRUE
  )

  expect_equivalent(mfx_logit$AME_pointwise[-1], mfx_logit_num$est, tol = 1e-3)
  expect_equivalent(mfx_logit$AME_pointwise_var[-1], mfx_logit_num$se^2, tol = 1e-3)

  ind_est <- get_individual_effects(mfx_logit_num)
  expect_equivalent(
    do.call("rbind", split(ind_est$est, ind_est$obs)),
    mfx_logit$ME_pointwise[, -1],
    tol = 1e-3
  )

  expect_equivalent(
    do.call("rbind", split(ind_est$se^2, ind_est$obs)),
    mfx_logit$ME_pointwise_var[, -1],
    tol = 1e-3
  )

  test_print <- print(mfx_logit)
  test_print <- print(mfx_logit_num)
  
})
