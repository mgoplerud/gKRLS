if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  env_test <- "CRAN"
  set.seed(150) # CRAN SEED
}else{
  # If on local machine
  env_test <- 'local'
}

context("Test Nystrom approximation")

# Shared data for tests
N <- 200
set.seed(42)
x1 <- rnorm(N)
x2 <- rnorm(N)
x3 <- rnorm(N)
y_cont <- sin(x1) + 0.5 * x2 + rnorm(N, sd = 0.5)
y_bin <- rbinom(N, 1, plogis(0.5 * x1 + 0.3 * x2))
dat <- data.frame(y_cont, y_bin, x1, x2, x3)

test_that("Nystrom model fits and produces valid output", {

  fit_nys <- gam(y_cont ~ s(x1, x2, x3, bs = "gKRLS",
    xt = gKRLS(sketch_method = "nystrom", sketch_multiplier = 5)),
    data = dat)

  expect_s3_class(fit_nys, "gam")
  expect_true(inherits(fit_nys$smooth[[1]], "gKRLS.smooth"))
  expect_true(summary(fit_nys)$s.table[, "edf"] > 1)
  # Landmarks are not rows of X, so subsampling_id should be NULL
  expect_null(fit_nys$smooth[[1]]$subsampling_id)
})

test_that("Nystrom fit is broadly similar to subsampling", {

  fit_sub <- gam(y_cont ~ s(x1, x2, x3, bs = "gKRLS",
    xt = gKRLS(sketch_method = "subsampling", sketch_multiplier = 5)),
    data = dat)
  fit_nys <- gam(y_cont ~ s(x1, x2, x3, bs = "gKRLS",
    xt = gKRLS(sketch_method = "nystrom", sketch_multiplier = 5)),
    data = dat)

  # R-squared should be in the same ballpark
  r2_sub <- 1 - sum(residuals(fit_sub)^2) / sum((dat$y_cont - mean(dat$y_cont))^2)
  r2_nys <- 1 - sum(residuals(fit_nys)^2) / sum((dat$y_cont - mean(dat$y_cont))^2)
  expect_true(abs(r2_sub - r2_nys) < 0.15)

  # In-sample predictions should correlate highly
  expect_true(cor(fitted(fit_sub), fitted(fit_nys)) > 0.85)
})

test_that("Nystrom out-of-sample prediction works", {

  fit_nys <- gam(y_cont ~ s(x1, x2, x3, bs = "gKRLS",
    xt = gKRLS(sketch_method = "nystrom", sketch_multiplier = 5)),
    data = dat)

  newdata <- data.frame(x1 = rnorm(20), x2 = rnorm(20), x3 = rnorm(20))
  pred <- predict(fit_nys, newdata = newdata, se.fit = TRUE)

  expect_length(pred$fit, 20)
  expect_length(pred$se.fit, 20)
  expect_true(all(is.finite(pred$fit)))
  expect_true(all(pred$se.fit > 0))
})

test_that("Nystrom prediction aligns with subsampling on held-out data", {

  dat_test <- data.frame(x1 = rnorm(50), x2 = rnorm(50), x3 = rnorm(50))

  fit_sub <- gam(y_cont ~ s(x1, x2, x3, bs = "gKRLS",
    xt = gKRLS(sketch_method = "subsampling", sketch_multiplier = 5)),
    data = dat)
  fit_nys <- gam(y_cont ~ s(x1, x2, x3, bs = "gKRLS",
    xt = gKRLS(sketch_method = "nystrom", sketch_multiplier = 5)),
    data = dat)

  pred_sub <- predict(fit_sub, newdata = dat_test)
  pred_nys <- predict(fit_nys, newdata = dat_test)

  expect_true(cor(pred_sub, pred_nys) > 0.85)
})

test_that("Nystrom marginal effects (AME) are close to subsampling", {

  fit_sub <- gam(y_cont ~ s(x1, x2, x3, bs = "gKRLS",
    xt = gKRLS(sketch_method = "subsampling", sketch_multiplier = 5)),
    data = dat)
  fit_nys <- gam(y_cont ~ s(x1, x2, x3, bs = "gKRLS",
    xt = gKRLS(sketch_method = "nystrom", sketch_multiplier = 5)),
    data = dat)

  mfx_sub <- calculate_effects(fit_sub, variables = c("x1", "x2", "x3"))
  mfx_nys <- calculate_effects(fit_nys, variables = c("x1", "x2", "x3"))

  expect_equal(nrow(mfx_sub), 3)
  expect_equal(nrow(mfx_nys), 3)

  # Estimates should be in the same direction with similar magnitude
  for (v in c("x1", "x2", "x3")) {
    est_sub <- mfx_sub$est[mfx_sub$variable == v]
    est_nys <- mfx_nys$est[mfx_nys$variable == v]
    se_sub  <- mfx_sub$se[mfx_sub$variable == v]
    se_nys  <- mfx_nys$se[mfx_nys$variable == v]
    # Estimates should not diverge by more than a few pooled SEs
    pooled_se <- sqrt(se_sub^2 + se_nys^2)
    expect_true(abs(est_sub - est_nys) < 5 * pooled_se,
      label = paste("AME for", v, "diverges beyond 5 pooled SEs"))
    # SEs should be in the same order of magnitude
    expect_true(se_nys / se_sub > 0.3 && se_nys / se_sub < 3,
      label = paste("SE ratio for", v, "is extreme"))
  }
})

test_that("Nystrom individual marginal effects correlate with subsampling", {

  fit_sub <- gam(y_cont ~ s(x1, x2, x3, bs = "gKRLS",
    xt = gKRLS(sketch_method = "subsampling", sketch_multiplier = 5)),
    data = dat)
  fit_nys <- gam(y_cont ~ s(x1, x2, x3, bs = "gKRLS",
    xt = gKRLS(sketch_method = "nystrom", sketch_multiplier = 5)),
    data = dat)

  ime_sub <- calculate_effects(fit_sub, variables = "x1", individual = TRUE)
  ime_nys <- calculate_effects(fit_nys, variables = "x1", individual = TRUE)

  ind_sub <- attr(ime_sub, "individual")$est
  ind_nys <- attr(ime_nys, "individual")$est

  expect_length(ind_sub, N)
  expect_length(ind_nys, N)
  expect_true(cor(ind_sub, ind_nys) > 0.4)
})

test_that("Nystrom interaction effects are computed without error", {

  fit_nys <- gam(y_cont ~ s(x1, x2, x3, bs = "gKRLS",
    xt = gKRLS(sketch_method = "nystrom", sketch_multiplier = 5)),
    data = dat)

  int_nys <- calculate_interactions(fit_nys, variables = list(c("x1", "x2")))

  expect_true(nrow(int_nys) > 0)
  expect_true("AMIE" %in% int_nys$QOI)
  expect_true(all(is.finite(int_nys$est)))
  expect_true(all(is.finite(int_nys$se)))
})

test_that("Nystrom robust SEs (HC1) work", {

  fit_nys <- gam(y_cont ~ s(x1, x2, x3, bs = "gKRLS",
    xt = gKRLS(sketch_method = "nystrom", sketch_multiplier = 5)),
    data = dat)

  vcov_robust <- vcovHC(fit_nys, type = "HC1")

  expect_true(is.matrix(vcov_robust))
  expect_true(all(diag(vcov_robust) > 0))

  # AME with robust SEs should also work
  mfx_robust <- calculate_effects(fit_nys, variables = "x1", vcov = vcov_robust)
  expect_true(mfx_robust$se > 0)
  expect_true(is.finite(mfx_robust$est))
})

test_that("Nystrom works with logistic regression", {

  fit_nys <- gam(y_bin ~ s(x1, x2, x3, bs = "gKRLS",
    xt = gKRLS(sketch_method = "nystrom", sketch_multiplier = 5)),
    data = dat, family = binomial(link = "logit"))

  expect_s3_class(fit_nys, "gam")

  # Predictions should be probabilities
  pred <- predict(fit_nys, type = "response")
  expect_true(all(pred > 0 & pred < 1))

  # Marginal effects should work
  mfx <- calculate_effects(fit_nys, variables = "x1")
  expect_true(is.finite(mfx$est))
  expect_true(mfx$se > 0)
})

test_that("Nystrom bandwidth calibration works", {

  fit_cal <- gam(y_cont ~ s(x1, x2, x3, bs = "gKRLS",
    xt = gKRLS(sketch_method = "nystrom", bandwidth = "calibrate",
      sketch_multiplier = 5)),
    data = dat)

  cal_info <- get_calibration_information(fit_cal)
  expect_true(nrow(cal_info) == 1)
  expect_true(cal_info$bandwidth > 0)
  expect_true(is.finite(cal_info$time))
})

test_that("Nystrom custom nystrom_args are respected", {

  fit_custom <- gam(y_cont ~ s(x1, x2, x3, bs = "gKRLS",
    xt = gKRLS(sketch_method = "nystrom",
      nystrom_args = list(batch_size = 50, max_iter = 200, tol = 1e-8))),
    data = dat)

  expect_s3_class(fit_custom, "gam")
  expect_true(all(is.finite(fitted(fit_custom))))
})

test_that("gKRLS validates nystrom_args", {

  expect_error(
    gKRLS(sketch_method = "nystrom", nystrom_args = list(bad_arg = 1)),
    "Unknown nystrom_args"
  )
})

test_that("Nystrom with sketch_size_raw works", {

  fit_raw <- gam(y_cont ~ s(x1, x2, bs = "gKRLS",
    xt = gKRLS(sketch_method = "nystrom",
      sketch_multiplier = NULL, sketch_size_raw = 20)),
    data = dat)

  expect_s3_class(fit_raw, "gam")
  expect_true(all(is.finite(fitted(fit_raw))))
})

test_that("Nystrom with demean_kernel works", {

  fit_dm <- gam(y_cont ~ s(x1, x2, bs = "gKRLS",
    xt = gKRLS(sketch_method = "nystrom", demean_kernel = TRUE)),
    data = dat)

  expect_s3_class(fit_dm, "gam")
  expect_true(all(is.finite(fitted(fit_dm))))
})

test_that("Nystrom errors when sketch_size > N", {

  expect_error(
    suppressWarnings(
      gam(y_cont ~ s(x1, x2, bs = "gKRLS",
        xt = gKRLS(sketch_method = "nystrom",
          sketch_multiplier = NULL, sketch_size_raw = N + 10)),
        data = dat)
    ),
    "sketch_size must be less than N"
  )
})

test_that("minibatch_kmeans_cpp returns correct dimensions", {

  X <- matrix(rnorm(100 * 3), ncol = 3)
  centers <- minibatch_kmeans_cpp(X, k = 5, batch_size = 20,
    max_iter = 50, tol = 1e-6)

  expect_equal(nrow(centers), 5)
  expect_equal(ncol(centers), 3)
  expect_true(all(is.finite(centers)))
})

test_that("minibatch_kmeans_cpp returns X when k >= N", {

  X <- matrix(rnorm(10 * 2), ncol = 2)
  centers <- minibatch_kmeans_cpp(X, k = 15, batch_size = 5,
    max_iter = 10, tol = 1e-6)

  expect_equal(nrow(centers), 10)
  expect_equal(ncol(centers), 2)
})

if (env_test != "CRAN") {

  test_that("Nystrom Monte Carlo: estimates are stable across replications", {

    n_rep <- 5
    ame_sub <- ame_nys <- numeric(n_rep)

    for (r in seq_len(n_rep)) {
      set.seed(1000 + r)
      y_r <- sin(x1) + 0.5 * x2 + rnorm(N, sd = 0.5)
      dat_r <- data.frame(y_cont = y_r, x1, x2, x3)

      fit_sub_r <- gam(y_cont ~ s(x1, x2, x3, bs = "gKRLS",
        xt = gKRLS(sketch_method = "subsampling", sketch_multiplier = 5)),
        data = dat_r)
      fit_nys_r <- gam(y_cont ~ s(x1, x2, x3, bs = "gKRLS",
        xt = gKRLS(sketch_method = "nystrom", sketch_multiplier = 5)),
        data = dat_r)

      mfx_sub_r <- calculate_effects(fit_sub_r, variables = "x1")
      mfx_nys_r <- calculate_effects(fit_nys_r, variables = "x1")

      ame_sub[r] <- mfx_sub_r$est
      ame_nys[r] <- mfx_nys_r$est
    }

    # Means should be close
    expect_true(abs(mean(ame_sub) - mean(ame_nys)) < 0.3)
    # SDs (across replications) should be similar order of magnitude
    expect_true(sd(ame_nys) / sd(ame_sub) > 0.2 && sd(ame_nys) / sd(ame_sub) < 5)
  })

}
