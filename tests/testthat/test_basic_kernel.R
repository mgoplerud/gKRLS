
context('Test basic kernel operations and edge cases')

test_that("Check kernel CPP aligns with direct", {
  
  N <- 50
  X <- cbind(matrix(rnorm(N * 2), ncol = 2), rbinom(N, 1, 0.5))
  y <- X %*% rnorm(ncol(X))
  
  K_cpp <- create_sketched_kernel(X_test = X, 
    X_train = X[1:15,], tS = diag(15),
    bandwidth = 2)

  K_base <- exp(-as.matrix(dist(X)^2)/2)[,1:15]
  
  K_outer <- exp(-base_kernel(X, X[1:15,])/2)
  
  expect_equivalent(K_cpp, K_outer, tol = 1e-6)
  expect_equivalent(K_cpp, K_base, tol = 1e-6)
  
  expect_error(create_sketched_kernel(X_test = X, X_train = X, bandwidth = 2,
                                      tS = diag(3)))
  
  S <- create_sketch_matrix(N = nrow(X), sketch_size = 3, sketch_method = 'gaussian')
  K_cpp <- create_sketched_kernel(X_test = X, 
                                  X_train = X, tS = t(S),
                                  bandwidth = 2)
  expect_equivalent(K_cpp, exp(-base_kernel(X, X)/2) %*% S)
  expect_equivalent(K_cpp, exp(-as.matrix(dist(X)^2)/2) %*% S)

})

test_that("Test everything works when kernel has limited columns", {
  
  N <- 50
  X <- cbind(matrix(rnorm(N * 2), ncol = 2), rbinom(N, 1, 0.5))
  y <- X %*% rnorm(ncol(X))
  y <- rpois(length(y), exp(y))
  
  fit_gam <- suppressWarnings(gam(y ~ s(X1, X2, X3, bs = 'gKRLS',
    xt = gKRLS(standardize = 'Mahalanobis', truncate.eigen.tol = 1e-4, sketch_prob = 0.2, sketch_size = N, sketch_method = 'bernoulli')), data = data.frame(X, y),
    family = poisson()))
  
  expect_equal(length(coef(fit_gam)) - 1, ncol(fit_gam$smooth[[1]]$sketch_matrix))
  
  evalues <- eigen(exp(-as.matrix(dist(X)^2)/ncol(X)))$values
  fit_single <- gam(y ~ s(X1, X2, X3, bs = 'gKRLS',
    xt = gKRLS(standardize = 'none', 
               truncate.eigen.tol = max(evalues) - 0.001, sketch_method = 'none')), data = data.frame(X, y),
    family = poisson())
  expect_equal(length(coef(fit_single)), 2)
  
})