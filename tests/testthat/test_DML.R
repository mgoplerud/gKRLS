
if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  env_test <- "CRAN"
  set.seed(126)
}else{
  # If on local machine
  env_test <- 'local'
}

test_that("Test DoubleML", {
  N <- 100
  x1 <- rnorm(N)
  x2 <- rbinom(N, size = 1, prob = .2)
  y <- x1^3 - 0.5 * x2 + rnorm(N, 0, 1)
  y <- y * 10
  X <- cbind(x1, x2, x1 + x2 * 3)
  X <- cbind(X, "x3" = rexp(nrow(X)))

  if (requireNamespace("DoubleML", quietly = TRUE)) {
    require(DoubleML)
    double_bam_1 <- LearnerRegrBam$new()
    double_bam_1$param_set$values$formula <- ~ s(x1, x3, bs = "gKRLS", xt = gKRLS(sketch_multiplier = NULL, sketch_size_raw = 2))
    double_bam_2 <- LearnerClassifBam$new()
    double_bam_2$param_set$values$formula <- ~ s(x1, x3, bs = "gKRLS", xt = gKRLS(sketch_multiplier = NULL, sketch_size_raw = 2))

    dml_data <- DoubleMLData$new(
      data = data.frame(X, y),
      x_cols = c("x1", "x3"), y_col = "y",
      d_cols = "x2"
    )

    dml_est <- DoubleMLIRM$new(
      data = dml_data,
      n_folds = 2,
      ml_g = double_bam_1,
      ml_m = double_bam_2
    )$fit()

    expect_s3_class(dml_est, "DoubleML")
  }
})

test_that("Test SuperLearner", {
  N <- 100
  x1 <- rnorm(N)
  x2 <- rbinom(N, size = 1, prob = .2)
  y <- x1^3 - 0.5 * x2 + rnorm(N, 0, 1)
  y <- y * 10
  X <- cbind(x1, x2, x1 + x2 * 3)

  if (requireNamespace("SuperLearner", quietly = TRUE)) {
    require(SuperLearner)
    sl_m <- function(...) {
      SL.mgcv(formula = ~ x1 + x2, ...)
    }
    fit_SL <- SuperLearner::SuperLearner(
      Y = y, obsWeights = rep(1, nrow(X)),
      X = data.frame(X),
      SL.library = "sl_m"
    )
    expect_s3_class(fit_SL, "SuperLearner")

    sl_m <- function(...) {
      SL.mgcv(bam = TRUE, formula = ~ x1 + x2, ...)
    }
    fit_SL <- SuperLearner::SuperLearner(
      Y = as.numeric(y > mean(y)), obsWeights = rep(1, nrow(X)),
      X = data.frame(X), family = "binomial",
      SL.library = "sl_m"
    )
    expect_s3_class(fit_SL, "SuperLearner")

    pred <- predict(fit_SL, newdata = data.frame(X))
    expect_length(pred, n = 2)
  }
})
