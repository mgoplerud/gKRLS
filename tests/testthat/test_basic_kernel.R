
context("Test basic kernel operations and edge cases")

test_that("Check kernel CPP aligns with direct", {
  N <- 50
  X <- cbind(matrix(rnorm(N * 2), ncol = 2), rbinom(N, 1, 0.5))
  y <- X %*% rnorm(ncol(X))

  K_cpp <- create_sketched_kernel(
    X_test = X,
    X_train = X[1:15, ], tS = diag(15),
    bandwidth = 2
  )

  K_base <- exp(-as.matrix(dist(X)^2) / 2)[, 1:15]

  K_outer <- exp(-base_kernel(X, X[1:15, ]) / 2)

  expect_equivalent(K_cpp, K_outer, tol = 1e-6)
  expect_equivalent(K_cpp, K_base, tol = 1e-6)

  expect_error(create_sketched_kernel(
    X_test = X, X_train = X, bandwidth = 2,
    tS = diag(3)
  ))

  S <- create_sketch_matrix(N = nrow(X), sketch_size = 3, sketch_method = "gaussian")
  K_cpp <- create_sketched_kernel(
    X_test = X,
    X_train = X, tS = t(S),
    bandwidth = 2
  )
  expect_equivalent(K_cpp, exp(-base_kernel(X, X) / 2) %*% S)
  expect_equivalent(K_cpp, exp(-as.matrix(dist(X)^2) / 2) %*% S)
})

test_that("Test everything works when kernel has limited columns", {
  N <- 50
  X <- cbind(matrix(rnorm(N * 2), ncol = 2), rbinom(N, 1, 0.5))
  y <- X %*% rnorm(ncol(X))
  y <- rpois(length(y), exp(y))

  fit_gam <- suppressWarnings(gam(y ~ s(X1, X2, X3,
    bs = "gKRLS",
    xt = gKRLS(
      standardize = "Mahalanobis", truncate.eigen.tol = 1e-4,
      sketch_prob = 0.2, sketch_multiplier = NULL,
      sketch_size_raw = N,
      sketch_method = "bernoulli"
    )
  ),
  data = data.frame(X, y),
  family = poisson()
  ))

  expect_equal(length(coef(fit_gam)) - 1, ncol(fit_gam$smooth[[1]]$sketch_matrix))

  evalues <- eigen(exp(-as.matrix(dist(X)^2) / ncol(X)))$values
  fit_single <- gam(y ~ s(X1, X2, X3,
    bs = "gKRLS",
    xt = gKRLS(
      standardize = "none",
      truncate.eigen.tol = max(evalues) - 0.001, sketch_method = "none"
    )
  ),
  data = data.frame(X, y),
  family = poisson()
  )
  expect_equal(length(coef(fit_single)), 2)
})


test_that("Test sketch size options work as anticipated", {
  N <- 50
  X <- cbind(matrix(rnorm(N * 2), ncol = 2), rbinom(N, 1, 0.5))
  y <- X %*% rnorm(ncol(X))
  fit_one <- gam(y ~ s(X1, X2, X3,
    bs = "gKRLS",
    xt = gKRLS(sketch_multiplier = 3.255, remove_instability = FALSE)
  ), data = data.frame(X, y))
  expect_equal(floor(ceiling(N^(1 / 3)) * 3.255), nrow(fit_one$smooth[[1]]$sketch_matrix))
  expect_s3_class(fit_one, "gam")
  expect_error(
    gam(y ~ s(X1, X2, X3,
      bs = "gKRLS",
      xt = gKRLS(sketch_size_raw = 1)
    ), data = data.frame(X, y))
  )
  fit_two <- gam(y ~ s(X1, X2, X3,
    bs = "gKRLS",
    xt = gKRLS(sketch_size_raw = 3, sketch_multiplier = NULL)
  ), data = data.frame(X, y))
  expect_equal(nrow(fit_two$smooth[[1]]$sketch_matrix), 3)
})


test_that("Test everything works when kernel has categorical/factor variables", {

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
  expect_equivalent(mm, fit_gam$smooth[[1]]$X_train)
  
  
  v1 <- predict(fit_gam, 
          newdata = data.frame(X1 = -5:5, X2 = 0:10,
                     X3 = 2, X4 = 'a', X5 = 'South', X6 = 'CA',
                     X7 = 1.5, stringsAsFactors = F))
  v2 <- predict(fit_gam, 
                newdata = data.frame(X1 = -5:5, X2 = 0:10,
                                     X3 = 2, X4 = 'a', X5 = 'South', X6 = 'CA',
                                     X7 = 1.5, stringsAsFactors = T))
  expect_equal(v1, v2)
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
  expect_equal(v1, v2)
  
  fit_gam2 <- gam(y ~ 1 +
    s(X1, X2, X3, X4, X5, X6, X7, bs = "gKRLS"), 
    data = df
  )
  predict(fit_gam2, newdata = df[1:5,])
  expect_s3_class(calculate_effects(fit_gam2, variables = 'X5'), 'gKRLS_mfx')
  expect_error(legacy_marginal_effect(fit_gam2, newdata = df[1:5,]),
               regexp = 'factors are provided')
})
  

