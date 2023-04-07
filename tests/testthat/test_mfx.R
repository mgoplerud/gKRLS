if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  env_test <- "CRAN"
  set.seed(135)
}else{
  # If on local machine
  env_test <- 'local'
}


test_that("mfx with degenerate Mahalanobis", {
  
  N <- 100
  x1 <- rnorm(N)
  x2 <- rbinom(N, size = 1, prob = .2)
  y <- x1^3 - 0.5 * x2 + rnorm(N, 0, 1)
  y <- y * 10
  X <- cbind(x1, x2, x1 + x2 * 3)
  colnames(X) <- paste0("x", 1:ncol(X))
  
  for (std in c('Mahalanobis', 'scaled')){
    fit_gKRLS <- gam(y ~ 0 + s(x1, x2, x3,
                               bs = "gKRLS",
                               xt = gKRLS(standardize = std)
    ),
    family = gaussian(), method = "REML", data = data.frame(y, X)
    )
    expect_equivalent(fit_gKRLS$smooth[[1]]$term, c("x1", "x2", "x3"))
    if (std != 'scaled'){
      expect_equivalent(ncol(fit_gKRLS$smooth[[1]]$X_train), 2)
    }
    mfx_num <- calculate_effects(fit_gKRLS, data = data.frame(X), continuous_type = "deriv")
    mfx_legacy <- legacy_marginal_effect(fit_gKRLS,
                                         newdata = data.frame(X),
                                         keep = c("x1", "x2", "x3")
    )
    
    expect_equivalent(mfx_num$est, mfx_legacy$AME_pointwise, tol = 1e-3)
    expect_equivalent(mfx_num$se, sqrt(mfx_legacy$AME_pointwise_var), tol = 1e-3)
  }
})

test_that("Test MFX", {
  N <- 100
  x1 <- rnorm(N)
  x2 <- rbinom(N, size = 1, prob = .2)
  y <- x1^3 - 0.5 * x2 + rnorm(N, 0, 1)
  y <- y * 10
  X <- cbind(x1, x2, x1 + x2 * 3)
  colnames(X) <- paste0("x", 1:ncol(X))
  fit_gKRLS <- gam(y ~ 0 + s(x1, x2, x3,
    bs = "gKRLS",
    xt = gKRLS(standardize = "scaled", sketch_method = "gaussian")
  ),
  family = gaussian(), method = "REML", data = data.frame(y, X)
  )
  expect_equivalent(fit_gKRLS$smooth[[1]]$term, c("x1", "x2", "x3"))
  expect_equivalent(ncol(fit_gKRLS$smooth[[1]]$X_train), 3)

  mfx_num <- calculate_effects(fit_gKRLS, data = data.frame(X), continuous_type = "deriv")
  mfx_legacy <- legacy_marginal_effect(fit_gKRLS,
    newdata = data.frame(X),
    keep = c("x1", "x2", "x3")
  )

  expect_equivalent(mfx_num$est, mfx_legacy$AME_pointwise, tol = 1e-3)
  expect_equivalent(mfx_num$se, sqrt(mfx_legacy$AME_pointwise_var), tol = 1e-3)
})

test_that("test 'calculate_effects'", {
  
  N <- 100
  x1 <- rnorm(N, sd = 1.6)
  x2 <- rbinom(N, size = 1, prob = .2)
  s <- sample(letters[1:5], N, replace = T)
  X <- data.frame(x1, x2, s = factor(s), 
    l = sample(c(TRUE,FALSE), N, replace = T), stringsAsFactors = F)
  colnames(X) <- paste0("x", 1:ncol(X))
  X$x4 <- factor(X$x4)
  
  y <- x1^3 - 0.5 * x2 + 1/5 * match(X$x3, letters) + rnorm(N, 0, 1)
  y <- y * 10
  
  fit_gKRLS <- gam(y ~ s(x1, x2, x3, x4, bs = "gKRLS"),
    family = gaussian(), method = "REML", data = data.frame(y, X)
  )

  factor_test <- calculate_effects(model = fit_gKRLS, variables = "x3")
  expect_null(factor_test$individual)
  factor_test <- calculate_effects(
    model = fit_gKRLS, variables = "x3",
    individual = TRUE, conditional = data.frame(x1 = c(-1, 1))
  )
  expect_false(is.null(get_individual_effects(factor_test)))
  expect_equal(nrow(factor_test), 8)

  logical_test <- calculate_effects(model = fit_gKRLS, variables = "x4")
  custom_cont_test <- calculate_effects(
    model = fit_gKRLS, variables = "x1",
      continuous_type = list('x1' = c(-1, 1)))
  
  # Test alignment with derivative and second derivative
  
  fit_first_deriv <- calculate_effects(model = fit_gKRLS, 
      variables = 'x1', continuous_type = 'derivative', individual = TRUE)  
  fit_second_deriv <- calculate_effects(model = fit_gKRLS, 
      variables = 'x1', continuous_type = 'second_derivative', individual = TRUE)  
  
  
  obj <- fit_gKRLS$smooth[[1]]
  S <- obj$sketch_matrix
  
  X_test <- create_data_gKRLS(
    term_levels = obj$term_levels, 
    term_class = obj$term_class, 
    data_term = data.frame(X)[obj$term], 
    terms = obj$term, 
    allow.missing.levels = TRUE)
  if (identical(obj$xt$sketch_method, 'gaussian')){
    X_nystrom <- X_test
  }else{
    X_nystrom <- X_test[obj$subsampling_id,, drop = F]
  }
  
  obj$xt$return_raw <- TRUE
  K <- create_sketched_kernel(
    X_test = Predict.matrix.gKRLS.smooth(object = obj, data = data.frame(X)), 
    X_train = as.matrix(obj$X_train),
    S = diag(nrow(obj$X_train)),
    bandwidth = obj$bandwidth
  )
  
  W <- obj$std_train$whiten
  W2 <- tcrossprod(W)
  man_first_deriv <- -2/obj$bandwidth * as.vector( (K * outer( (X_test %*% W2)[,1], (X_nystrom %*% W2)[,1], FUN=function(x,y){x-y})) %*% t(S) %*% coef(fit_gKRLS)[-1] )

  man_second_deriv_a <- 4/obj$bandwidth^2 * as.vector( (K * outer( (X_test %*% W2)[,1], (X_nystrom %*% W2)[,1], FUN=function(x,y){x-y})^2 ) %*% t(S) %*% coef(fit_gKRLS)[-1] )
  man_second_deriv_b <- -2/obj$bandwidth * as.vector( (K * diag(W2)[1]) %*% t(S) %*% coef(fit_gKRLS)[-1] )
  man_second_deriv <- man_second_deriv_a + man_second_deriv_b
  
  expect_equivalent(get_individual_effects(fit_first_deriv)$est, man_first_deriv, tol = 1e-5)  
  expect_equivalent(get_individual_effects(fit_second_deriv)$est, man_second_deriv, tol = 1e-5)  
  
  SE_MAT <- -2/obj$bandwidth * ( (K * outer( (X_test %*% W2)[,1], (X_nystrom %*% W2)[,1], FUN=function(x,y){x-y})) %*% t(S))
  man_first_deriv_se <- sqrt(rowSums( (SE_MAT %*% vcov(fit_gKRLS)[-1,-1]) * SE_MAT ))
  expect_equivalent(get_individual_effects(fit_first_deriv)$se, man_first_deriv_se, tol = 1e-5)  

  
  SE_MAT2 <- 4/obj$bandwidth^2 * (K * outer( (X_test %*% W2)[,1], (X_nystrom %*% W2)[,1], FUN=function(x,y){x-y})^2 ) %*% t(S) +
    -2/obj$bandwidth * (K * diag(W2)[1]) %*% t(S)
  man_second_deriv_se <- sqrt(rowSums( (SE_MAT2 %*% vcov(fit_gKRLS)[-1,-1]) * SE_MAT2 ))
  expect_equivalent(get_individual_effects(fit_second_deriv)$se, man_second_deriv_se, tol = 1e-5)  
  
})

test_that("test logical and binary", {

  N <- 100
  x1 <- rnorm(N, sd = 1.6)
  x2 <- rbinom(N, size = 1, prob = .2)
  s <- sample(letters[1:5], N, replace = T)
  X <- data.frame(x1, x2, s = factor(s), 
                  l = sample(c(TRUE,FALSE), N, replace = T), stringsAsFactors = F)
  colnames(X) <- paste0("x", 1:ncol(X))

  y <- x1^3 - 0.5 * x2 + 1/5 * match(X$x3, letters) + rnorm(N, 0, 1)
  y <- y * 10
  
  fit_gKRLS <- suppressWarnings(gam(list(y ~ x4 * x3 + s(x1, x2, x3, bs = "gKRLS"), ~ x4 * x3 + s(x1,x2)),
    family = gaulss(), method = "REML", data = data.frame(y, X)
  ))
  
  est_effects <- calculate_effects(fit_gKRLS, 
    variables = list(c('x4', 'x3')))
  
  copy_X <- data.frame(X)  
  copy_X$x4 <- FALSE  
  copy_X$x3 <- 'a'
  pred_base <- predict(fit_gKRLS, newdata = copy_X, type = 'response')
  pred_base <- colMeans(pred_base)
  for (v in c('b', 'c', 'd')){
    copy_X$x3 <- v
    copy_X$x4 <- TRUE  
    pred_alt <- colMeans(predict(fit_gKRLS, newdata = copy_X, type = 'response'))
    expect_equal(pred_alt - pred_base, subset(est_effects, grepl(variable, pattern=paste0('factor.x3.', v)))$est)
  }
})

