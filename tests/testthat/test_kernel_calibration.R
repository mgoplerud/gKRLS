if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  env_test <- "CRAN"
  set.seed(555) # CRAN SEED
}else{
  # If on local machine
  env_test <- 'local'
}

context("Test calibration of kernel works as expected")

test_that("Check calibration for unsketched", {
  
  N <- 50
  X <- cbind(matrix(rnorm(N * 2), ncol = 2), rbinom(N, 1, 0.5))
  y <- X %*% rnorm(ncol(X))
  
  K <- create_sketched_kernel(
    X_test = X,
    X_train = X, S = diag(N),
    bandwidth = 1, raw = TRUE
  )
  K_base <- as.matrix(dist(X)^2)
  expect_equivalent(K, K_base)
  
  f <- function(ln_b){
    K_b <- exp(-K/exp(ln_b))
    return(var(as.vector(K_b)))
  }
  direct_opt <- optimize(f = f, 
   interval = c(-10, 10), 
   maximum = TRUE)
  opt_1 <- calibrate_bandwidth(X = X, S = diag(N))  
  opt_2 <- calibrate_bandwidth(X = X, id_S = 1:N)
  expect_equal(opt_1, opt_2)  
  expect_equal(exp(direct_opt$maximum), opt_1)  
  
})
  

test_that("Check calibration for sub-sampling sketching", {
  
  N <- 50
  X <- cbind(matrix(rnorm(N * 2), ncol = 2), rbinom(N, 1, 0.5))
  y <- X %*% rnorm(ncol(X))
  
  id <- sample(1:N, 5)
  raw_KSt <- create_sketched_kernel(
    X_test = X,
    X_train = X[id,], S = diag(length(id)),
    bandwidth = 1, raw = TRUE
  )
  raw_KSt_base <- as.matrix(dist(X)^2)[,id]
  expect_equivalent(raw_KSt, raw_KSt_base)
  
  f <- function(ln_b){
    KSt <- exp(-raw_KSt_base/exp(ln_b))
    recons_K <- KSt %*% solve(KSt[id,]) %*% t(KSt)
    return(var(as.vector(recons_K)))
  }
  direct_opt <- optimize(f = f, 
                         interval = c(-10, 10), 
                         maximum = TRUE)
  
  S <- as.matrix(sparseMatrix(i = 1:length(id), j = id, x = 1, dims = c(length(id), N)))
  opt_1 <- calibrate_bandwidth(X = X, S =S)
  opt_2 <- calibrate_bandwidth(X = X, id_S = id)
  expect_equal(opt_1, opt_2)  
  expect_equal(exp(direct_opt$maximum), opt_1)  
  
})
  
test_that("Check calibration for Gaussian sketching", {
  
  N <- 50
  X <- cbind(matrix(rnorm(N * 2), ncol = 2), rbinom(N, 1, 0.5))
  y <- X %*% rnorm(ncol(X))
  S <- matrix(rnorm(5 * N), nrow = 5)
  raw_K <- as.matrix(dist(X)^2)

  f <- function(ln_b){
    KSt <- exp(-raw_K/exp(ln_b)) %*% t(S)
    recons_K <- KSt %*% solve(S %*% KSt) %*% t(KSt)
    return(var(as.vector(recons_K)))
  }
  direct_opt <- optimize(f = f, 
                         interval = c(-10, 10), 
                         maximum = TRUE)
  
  opt_1 <- calibrate_bandwidth(X = X, S =S)
  expect_equal(exp(direct_opt$maximum), opt_1)  
  
})

test_that("Using calibration with mgcv", {
  
  N <- 50
  X <- cbind(matrix(rnorm(N * 2), ncol = 2), rbinom(N, 1, 0.5))
  y <- X %*% rnorm(ncol(X))
  dat <- data.frame(X,y)  
  fit_dat <- gam(y ~ s(X2, bs = 'gKRLS') +
        s(X1,X3,bs = 'gKRLS', xt=gKRLS(bandwidth = 'calibrate')),
      data = dat)
  fit_dat_2 <- bam(y ~ s(X2, bs = 'gKRLS', xt =gKRLS(bandwidth = 'calibrate')) +
                   s(X1,X3,bs = 'gKRLS', xt=gKRLS(bandwidth = 'calibrate')),
                 data = dat)
  fit_simple <- gam(y ~ s(X1) + s(X2) + X3, data = dat)
  calib_1 <- get_calibration_information(fit_dat)
  expect_true(all(is.na(calib_1$time) == c(TRUE, FALSE)))
  expect_true(calib_1$bandwidth[1] == 1 & calib_1$bandwidth[2] != 2)
  calib_2 <- get_calibration_information(fit_dat_2)
  expect_false(any(is.na(calib_2$time)))
  expect_false(any(calib_2$bandwidth == c(1,2)))
  calib_3 <- get_calibration_information(fit_simple)  
  expect_identical(calib_3, data.frame())
  
})
