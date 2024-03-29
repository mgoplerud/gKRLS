if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  env_test <- "CRAN"
  set.seed(124) # CRAN SEED
}else{
  # If on local machine
  env_test <- 'local'
}


context("Test basic kernel operations and edge cases")

test_that("Check kernel CPP aligns with direct", {
  
  N <- 50
  X <- cbind(matrix(rnorm(N * 2), ncol = 2), rbinom(N, 1, 0.5))
  y <- X %*% rnorm(ncol(X))

  K_cpp <- create_sketched_kernel(
    X_test = X,
    X_train = X[1:15, ], S = diag(15),
    bandwidth = 2
  )

  K_base <- exp(-as.matrix(dist(X)^2) / 2)[, 1:15]

  K_outer <- exp(-base_kernel(X, X[1:15, ]) / 2)

  expect_equivalent(K_cpp, K_outer, tol = 1e-6)
  expect_equivalent(K_cpp, K_base, tol = 1e-6)

  expect_error(create_sketched_kernel(
    X_test = X, X_train = X, bandwidth = 2,
    S = diag(3)
  ))

  S <- create_sketch_matrix(N = nrow(X), sketch_size = 3, sketch_method = "gaussian")
  K_cpp <- create_sketched_kernel(
    X_test = X,
    X_train = X, S = S,
    bandwidth = 2
  )
  expect_equivalent(K_cpp, exp(-base_kernel(X, X) / 2) %*% t(S))
  expect_equivalent(K_cpp, exp(-as.matrix(dist(X)^2) / 2) %*% t(S))

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
  
  # Confirm that Predict.matrix.gKRLS.smooth works as expected
  expect_equivalent(
    as.vector(
      cbind(1, Predict.matrix.gKRLS.smooth(fit_gam$smooth[[1]], data = data.frame(X))) %*% coef(fit_gam)
    ), predict(fit_gam, newdata = data.frame(X))
  )
  
  expect_equal(length(coef(fit_gam)) - 1, nrow(fit_gam$smooth[[1]]$sketch_matrix))
  evalues <- eigen(exp(-as.matrix(dist(X)^2) / ncol(X)))$values
  fit_single <- gam(y ~ s(X1, X2, X3,
    bs = "gKRLS",
    xt = gKRLS(
      standardize = "none",
      truncate.eigen.tol = max(evalues) - 0.001, sketch_method = "none"
    )
  ), data = data.frame(X, y), family = poisson())
  
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

test_that("Test custom vector", {
  
  N <- 50
  X <- cbind(matrix(rnorm(N * 2), ncol = 2), rbinom(N, 1, 0.5))
  y <- X %*% rnorm(ncol(X))
  cust_vector <- sample(1:50, 5)
  fit_one <- gam(y ~ s(X1, X2, X3,
   bs = "gKRLS",
   xt = gKRLS(sketch_method = cust_vector, truncate.eigen.tol = 0.001)
  ), data = data.frame(X, y))
  expect_equal(length(fit_one$smooth[[1]]$subsampling_id), 5)
  
  # Check that prediction works as expected with eigentruncation
  predict_mat <- Predict.matrix.gKRLS.smooth(object = fit_one$smooth[[1]], 
                                             data = data.frame(X))
  eig_X <- eigen(cov(X))
  std_X <- sweep(X, MARGIN = 2, STATS = colMeans(X), FUN = '-')
  std_X <- as.matrix(std_X %*% eig_X$vectors %*% Diagonal(x=1/sqrt(eig_X$values)))
  K <- exp(-as.matrix(dist(std_X)^2) / ncol(std_X))
  P <- K[cust_vector, cust_vector] * N/5
  eig_P <- eigen(P, symmetric = TRUE)
  eig_P$values[eig_P$values < 1e-4] <- 0
  eig_nonzero <- which(eig_P$values != 0)
  TMat <- eig_P$vectors[,eig_nonzero] %*% Diagonal(x = 1/sqrt(eig_P$values[eig_nonzero]))
  custom_pred <- sqrt(N/5) * as.matrix(K[,cust_vector] %*% TMat)
  # Address issue of eigen sign-flip
  custom_pred <- custom_pred %*% Diagonal(x=sapply(1:ncol(custom_pred), FUN=function(i){cor(custom_pred[,i], predict_mat[,i])}))
  expect_equivalent(predict_mat, as.matrix(custom_pred))
  
  expect_s3_class(fit_one, "gam")
  fit_one <- gam(y ~ s(X1, X2, X3,
                       bs = "gKRLS",
                       xt = gKRLS(sketch_method = rep(1,4), remove_instability = FALSE)
  ), data = data.frame(X, y))
  expect_equal(length(fit_one$smooth[[1]]$subsampling_id), 4)
  expect_true(all(fit_one$smooth[[1]]$subsampling_id == 1))
  expect_s3_class(fit_one, "gam")
  v <- predict(fit_one, newdata = data.frame(X)[1:5,])
  expect_vector(v, size = 5)
  
  
})

test_that("Test polynomial works", {
  
  N <- 50
  X <- cbind(matrix(rnorm(N * 2), ncol = 2), rbinom(N, 1, 0.5))
  y <- X %*% rnorm(ncol(X))
  X <- data.frame(X)
  
  fit_three <- gam(y ~ s(X1, bs = 'unregpoly', xt = list(degree = 3)), data = X)  
  fit_three_lm <- lm(y ~ poly(X1, 3), data = X)
  expect_equivalent(coef(fit_three), coef(fit_three_lm))
  expect_equivalent(vcov(fit_three), vcov(fit_three_lm))

  df <- data.frame(X1 = rnorm(100))
  pred_three <- predict(fit_three, newdata = df, se.fit = TRUE)
  pred_three <- lapply(pred_three, as.numeric)
  pred_three_lm <- predict(fit_three_lm, newdata = df, se.fit = TRUE)
  expect_equivalent(pred_three$fit, pred_three_lm$fit)
  expect_equivalent(pred_three$se.fit, pred_three_lm$se.fit)
  
  fit_three <- gam(y ~ X2 + s(X1, bs = 'unregpoly', xt = list(raw = T, degree = 3)), data = X)  
  fit_three_lm <- lm(y ~ X2 + poly(X1, 3, raw = T), data = X)
  expect_equivalent(coef(fit_three), coef(fit_three_lm))
  expect_equivalent(vcov(fit_three), vcov(fit_three_lm))
  
  
})
