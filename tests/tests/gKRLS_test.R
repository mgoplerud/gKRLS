library(testthat)
est_gKRLS <- gKRLS(y = y, X = X, standardize = 'Mahalanobis',
                   sketch_method = 'gaussian',
                   family = gaussian())

pred_gKRLS <- predict(est_gKRLS, newdata = X)


fitted_internal <- fitted(est_gKRLS)
manual_predict <- predict(est_gKRLS, newdata = X)
expect_equivalent(manual_predict$fitted, fitted_internal)


N <- 200
x1 <- rnorm(N)
x2 <- rbinom(N,size=1,prob=.2)
y <- x1^3 + .5*x2 + rnorm(N,0,.15)
X <- cbind(x1, x2)
X_copy <- cbind(X, X, X)

est_gKRLS <- gKRLS(y = y, X = X, 
                   sketch_method = 'none',
                   family = gaussian())

est_copy <- gKRLS(y = y, X = X_copy, 
                   sketch_method = 'none',
                   family = gaussian())

expect_equal(est_gKRLS$fitted, est_copy$fitted)

plot(predict(est_copy, newdata = X_copy)$fitted, est_copy$fitted)
