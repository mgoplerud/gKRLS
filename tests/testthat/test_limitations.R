if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  env_test <- "CRAN"
  set.seed(129) # CRAN SEED
}else{
  # If on local machine
  env_test <- 'local'
}


test_that("does not run with prior weights", {
  
  # From mgcv
  n <- 400
  dat <- gamSim(1,n=n, verbose = FALSE)
  dat$f <- dat$f - mean(dat$f)
  
  alpha <- c(-Inf,-1,0,5,Inf)
  R <- length(alpha)-1
  y <- dat$f
  u <- runif(n)
  u <- dat$f + log(u/(1-u)) 
  for (i in 1:R) {
    y[u > alpha[i]&u <= alpha[i+1]] <- i
  }
  dat$y <- y
  
  fit <- gam(y ~ s(x1), weights = x2, data = dat)  
  expect_warning(calculate_effects(fit), regexp = 'ignores "weights"')
  fit <- gam(y ~ s(x1), weights = x2, data = dat, family = ocat(R=R))  
  expect_warning(calculate_effects(fit), regexp = 'ignores "weights"')
  dat$new_y <- dat$y - 1
  fit <- gam(list(new_y ~ s(x1), ~ x1, ~ s(x2)), weights = x2, data = dat, family = multinom(K=3))  
  expect_warning(calculate_effects(fit), regexp = 'ignores "weights"')
  
})

test_that("does not run with offset", {
  
  # From mgcv
  n <- 400
  dat <- gamSim(1,n=n, verbose=FALSE)
  dat$f <- dat$f - mean(dat$f)
  
  alpha <- c(-Inf,-1,0,5,Inf)
  R <- length(alpha)-1
  y <- dat$f
  u <- runif(n)
  u <- dat$f + log(u/(1-u)) 
  for (i in 1:R) {
    y[u > alpha[i]&u <= alpha[i+1]] <- i
  }
  dat$y <- y
  
  fit <- gam(y ~ offset(x2) + s(x1), data = dat)  
  expect_error(calculate_effects(fit), regexp = '"offset" not set')
  fit <- gam(y ~ offset(x2) + s(x1), data = dat, family = ocat(R=R))  
  expect_error(calculate_effects(fit), regexp = '"offset" not set')
  dat$new_y <- dat$y - 1
  fit <- gam(list(new_y ~ offset(x2) + s(x1), ~ x1, ~ s(x2)), data = dat, family = multinom(K=3))  
  expect_error(calculate_effects(fit), regexp = '"offset" not set')
  fit <- gam(list(new_y ~ offset(x2) + s(x1), ~ offset(x2) + x1, ~ s(x2)), data = dat, family = multinom(K=3))  
  expect_error(calculate_effects(fit), regexp = '"offset" not set')
  
})

test_that("test raw works and fails where expected", {
  
  dat <- gamSim(eg = 1, n =100, verbose=FALSE)
  dat$x0 <- dat$x0 > 0
  dat$x1 <- cut(dat$x1, 3)
  b <- gam(y ~ x0 + x1 + s(x2, by =x1), data = dat)
  
  est_effects <- calculate_effects(b, raw = TRUE)
  # Ordered by effect, raw_1, raw_0
  expect_true(all(sapply(split(est_effects$est, est_effects$variable), FUN=function(i){
    i[1] - (i[2] - i[3])
  }) == 0))
  
  est_effects <- calculate_effects(b, continuous_type = 'derivative', 
       conditional = data.frame(x1 = unique(dat$x1)), 
       variables = 'x2', raw = TRUE)
  
  expect_equivalent(sapply(split(est_effects$est, est_effects$x1), FUN=function(i){
    i[1] - (i[2] - i[3]) / (2 * 1e-7 * max(1, max(abs(dat$x2))))
  }), rep(0, 3))
  
  est_effects <- calculate_effects(b, continuous_type = 'second_derivative', 
       variables = 'x2', raw = TRUE, individual = T)
  # Check formula for second derivative lines up  
  expect_equivalent(
    est_effects$est[1],
    sum(est_effects$est[-1] * c(1, -2, 1))/(1e-7 * max(1, max(abs(dat$x2))))
  )
  individual_raw <- get_individual_effects(est_effects)
  expect_equivalent(sapply(split(individual_raw$est, individual_raw$obs), FUN=function(i){
    i[1] - sum(i[-1] * c(1, -2,1)) / (1e-7 * max(1, max(abs(dat$x2))))
  }), rep(0, length(unique(individual_raw$obs))))

  # Test for errors  
  expect_error(calculate_effects(b, continuous_type = 'predict', raw=T), regexp = 'raw must be FALSE')
  expect_error(calculate_interactions(b, raw = T), regexp = 'raw=T not permitted')
})

test_that("test sketching dimension and bam", {
  
  dat <- gamSim(eg = 1, n = 500, verbose = FALSE)
  fit_bam_default <- bam(y ~ s(x1, x2, bs = 'gKRLS'), data = dat)  
  fit_bam <- bam(y ~ s(x1, x2, bs = 'gKRLS'), data = dat, chunk.size = 100)  
  fit_bam_custom <- bam(y ~ s(x1, x2, bs = 'gKRLS', xt = gKRLS(sketch_size_raw = 5, sketch_multiplier = NULL)), data = dat, chunk.size = 100)  
  # Check that sketching dimension is governed by "chunk.size" to bam
  expect_equal(length(fit_bam_default$smooth[[1]]$subsampling_id), ceiling(nrow(dat)^(1/3)) * 5)
  expect_equal(length(fit_bam$smooth[[1]]$subsampling_id), ceiling(100^(1/3)) * 5)
  expect_equal(length(fit_bam_custom$smooth[[1]]$subsampling_id), 5)
  
  if (env_test != "CRAN"){
    dat <- gamSim(eg = 1, n = 15000, verbose = FALSE)
    fit_gam <- gam(y ~ s(x1, x2, bs = 'gKRLS'), data = dat)  
    fit_bam <- bam(y ~ s(x1, x2, bs = 'gKRLS'), data = dat)  
    fit_bam_large <- bam(y ~ s(x1, x2, bs = 'gKRLS'), data = dat, chunk.size = Inf)  
    expect_equal(length(fit_bam$smooth[[1]]$subsampling_id), ceiling(10^(4/3)) * 5)
    expect_equal(length(fit_gam$smooth[[1]]$subsampling_id), ceiling(nrow(dat)^(1/3)) * 5)
    expect_equal(length(fit_bam_large$smooth[[1]]$subsampling_id), ceiling(nrow(dat)^(1/3)) * 5)
  }
  
})

test_that("test behavior when 'data' is specified", {
  
  dat <- gamSim(eg = 4, n = 100, verbose = FALSE)
  dat$x0 <- dat$x0 > median(dat$x0)
  dat$fac <- factor(letters[dat$fac])
  dat$y <- rnorm(3)[match(dat$fac, letters)] + dat$y
  b <- gam(y ~ x0 + s(fac, bs = 're') + s(x2, fac, bs = 'gKRLS'), data = dat)
  
  fit_default <- calculate_effects(b)
  
  new_dat <- dat[1:5,]
  new_dat$x0 <- FALSE
  new_dat$fac <- c('e', 'd', 'd', 'e', 'e')
  fit_limited_data <- expect_message(
    expect_warning(calculate_effects(b, data = new_dat), regexp = 'not in original fit'),
    regexp = 'New levels found in factor fac.'
  )
  expect_equal(sum(grepl(fit_limited_data$variable, pattern='factor')), 1) 
  expect_equal(fit_limited_data$est[3], 0)
  expect_equal(fit_limited_data$t[3], NA_real_)
  
  grid_x <- quantile(new_dat$x2, c(0.25, 0.75))
  limited_x2 <- diff(sapply(grid_x, FUN=function(i){
    copy <- new_dat
    copy$x2 <- i
    suppressWarnings(mean(predict(b, newdata = copy)))
  }))
  limited_x0 <- diff(sapply(c(FALSE, TRUE), FUN=function(i){
    copy <- new_dat
    copy$x0 <- i
    suppressWarnings(mean(predict(b, newdata = copy)))
  }))
  # Check values align for using data but *not* manually
  expect_equivalent(fit_limited_data$est[1:2], c(limited_x2, limited_x0))
  
  fit_data_override <- suppressWarnings(suppressMessages(calculate_effects(b, 
    data = new_dat, use_original = TRUE))) 
  expect_equal(sum(grepl(fit_data_override$variable, pattern='factor')), 2) 
  

  override_x2 <- diff(sapply(quantile(dat$x2, c(0.25, 0.75)), FUN=function(i){
    copy <- new_dat
    copy$x2 <- i
    suppressWarnings(mean(predict(b, newdata = copy)))
  }))
  override_x0 <- diff(sapply(c(FALSE, TRUE), FUN=function(i){
    copy <- new_dat
    copy$x0 <- i
    suppressWarnings(mean(predict(b, newdata = copy)))
  }))
  expect_equivalent(fit_data_override$est[1:2], c(override_x2, override_x0))
  
  expect_false(isTRUE(all.equal(fit_limited_data$est[1:2], fit_default$est[1:2])))
  expect_false(isTRUE(all.equal(fit_data_override$est[1:2], fit_default$est[1:2])))
  
  expect_identical(attr(fit_limited_data, 'N_eff'), attr(fit_default, 'N_eff'))
  expect_identical(attr(fit_limited_data, 'N'), attr(fit_default, 'N'))
  
  # Check that it works as expected with conditional:
  fit_cond <- expect_message(
    expect_warning(calculate_effects(b, variables = 'x2', conditional = data.frame(x0 = FALSE), data = new_dat), regexp = 'not in original fit'),
    regexp = 'New levels found in factor fac.'
  )
  cond_x2 <- diff(sapply(quantile(new_dat$x2, c(0.25, 0.75)), FUN=function(i){
    copy <- new_dat
    copy$x2 <- i
    suppressWarnings(mean(predict(b, newdata = copy)))
  }))
  expect_equivalent(fit_cond$est, cond_x2)
  
})

test_that("test other cases where 'data' is provided", {
  
  dat <- gamSim(eg = 4, n = 100, verbose = FALSE)
  dat$x2 <- dat$x2
  dat$x0 <- dat$x0 > median(dat$x0)
  dat$fac <- factor(letters[dat$fac])
  dat$y <- rnorm(3)[match(dat$fac, letters)] + dat$y
  b <- gam(y ~ x0 + s(x2, fac, bs = 'gKRLS'), data = dat)
  
  if (env_test == 'CRAN'){
    loop_list <- c('derivative', 'second_derivative')
  }else{
    loop_list <- c('derivative', 'second_derivative', 'IQR', 'minmax', 'onesd')
  }
  
  for (v in loop_list){
    
    print(v)
    fit_default <- calculate_effects(b, continuous_type = v, variables = 'x2')
    
    new_dat <- dat[1:5,]
    new_dat$x0 <- FALSE
    new_dat$fac <- c('e', 'd', 'd', 'e', 'e')
    fit_limited_data <- expect_message(
      expect_warning(calculate_effects(b, variables = 'x2', data = new_dat, continuous_type = v), regexp = 'not in original fit'),
      regexp = 'New levels found in factor fac.'
    )
    
    if (v == 'derivative'){
      f <- function(s, x, i){
        step <- max(1, max(abs(x))) * 1e-7 * c(-1, 1)[i]
        return(s + step)
      }
    }else if (v == 'second_derivative'){
      f <- function(s, x, i){
        step <- sqrt(max(1,max(abs(x))) * 1e-7) * c(1, 0, -1)[i]
        return(s + step)
      }
    }else if (v == 'onesd'){
      f <- function(s, x, i){mean(x) + c(-1, 1)[i] * sd(x)}
    }else if (v == 'minmax'){
      f <- function(s, x, i){range(x)[i]}
    }else if (v == 'IQR'){
      f <- function(s, x, i){quantile(x, c(0.25, 0.75))[i]}
    }else{stop('...')}
    
    if (v == 'second_derivative'){
      grid_x <- 1:3
    }else{
      grid_x <- 1:2
    }
    
    if (v == 'second_derivative'){
      g <- function(x){sum(x * c(1,-2,1))}
    }else{
      g <- function(x){diff(x)}
    }
    
    limited_x2 <- g(sapply(grid_x, FUN=function(i){
      copy <- new_dat
      copy$x2 <- f(copy$x2, copy$x2,i)
      out <- suppressWarnings(mean(predict(b, newdata = copy)))
      if (v == 'derivative'){
        out <- 1/(2 * max(1, max(abs(new_dat$x2))) * 1e-7) * out
      }else if (v == 'second_derivative'){
        out <- 1/(max(1, max(abs(new_dat$x2))) * 1e-7) * out
      }else{}
      return(out)
    }))
    expect_equivalent(fit_limited_data$est[1], limited_x2)
    
    fit_data_override <- suppressWarnings(suppressMessages(calculate_effects(b, 
                                                                             data = new_dat, variables = 'x2', continuous_type = v, use_original = TRUE))) 
    
    override_x2 <- g(sapply(grid_x, FUN=function(i){
      copy <- new_dat
      copy$x2 <- f(copy$x2, dat$x2,i)
      out <- suppressWarnings(mean(predict(b, newdata = copy)))
      if (v == 'derivative'){
        out <- 1/(2 * max(1, max(abs(dat$x2))) * 1e-7) * out
      }else if (v == 'second_derivative'){
        out <- 1/(max(1, max(abs(dat$x2))) * 1e-7) * out
      }else{}
      return(out)
    }))
    expect_equivalent(fit_data_override$est[1], override_x2)
    
  }
  
})