test_that("Test that kernel bandwidth calibration works", {
  
  N <- 50
  X <- cbind(matrix(rnorm(N * 2), ncol = 2), rbinom(N, 1, 0.5))
  y <- X %*% rnorm(ncol(X))
  # Check that it runs with calibration
  fit_one <- gam(y ~ s(X1, X2, X3,
                       bs = "gKRLS", xt = gKRLS(bandwidth = "calibrate")), data = data.frame(X, y))
  # Extract the calibration information
  calib_info <- get_calibration_information(fit_one)
  
  fit_two <- gam(y ~ 
                   s(X2) +
      s(X1, X3, bs = "gKRLS", xt = gKRLS(bandwidth = "calibrate")) +
      s(X1, X2, X3, bs = "gKRLS", xt = gKRLS(bandwidth = "calibrate")), 
    data = data.frame(X, y), method = "GCV.Cp")
  calib_info <- get_calibration_information(fit_two)
  # Check that calibration works correctly
  expect_equal(calib_info$bandwidth, sapply(calib_info$smooth, FUN=function(i){
    fit_two$smooth[[i]]$bandwidth
  }))
  
  # Check the manually reconstructed kernel uses the correct, calibrated,
  # bandwidths
  X_test <- data.frame(cbind(matrix(rnorm(N * 2 * 2), ncol = 2), rbinom(N * 2, 1, 0.5)))
  pred_test <- predict(fit_two, newdata = X_test, type = 'terms')
  kern_two <- Predict.matrix(fit_two$smooth[[2]], data = X_test)
  kern_three <- Predict.matrix(fit_two$smooth[[3]], data = X_test)
 
  input_X <- sweep(as.matrix(X_test[,c('X1','X3')]), MARGIN = 2, 
        FUN='-', STATS = fit_two$smooth[[2]]$std_train$mean)
  input_X <- input_X %*% fit_two$smooth[[2]]$std_train$whiten
  test_K <- gKRLS:::create_sketched_kernel(
    X_test = as.matrix(input_X),
    X_train = fit_two$smooth[[2]]$X_train, 
    bandwidth = fit_two$smooth[[2]]$bandwidth,
    S = fit_two$smooth[[2]]$sketch_matrix)
  mgcv_K <- Predict.matrix(fit_two$smooth[[2]], data = X_test)  
  expect_equal(test_K, mgcv_K)

  input_X <- sweep(as.matrix(X_test[,c('X1','X2', 'X3')]), MARGIN = 2, 
                   FUN='-', STATS = fit_two$smooth[[3]]$std_train$mean)
  input_X <- input_X %*% fit_two$smooth[[3]]$std_train$whiten
  test_K <- gKRLS:::create_sketched_kernel(
    X_test = as.matrix(input_X),
    X_train = fit_two$smooth[[3]]$X_train, 
    bandwidth = fit_two$smooth[[3]]$bandwidth,
    S = fit_two$smooth[[3]]$sketch_matrix)
  mgcv_K <- Predict.matrix(fit_two$smooth[[3]], data = X_test)  
  expect_equal(test_K, mgcv_K)
  
})
