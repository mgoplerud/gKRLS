if (isTRUE(as.logical(Sys.getenv("CI")))){
  # If on CI
  env_test <- "CI"
}else if (!identical(Sys.getenv("NOT_CRAN"), "true")){
  # If on CRAN
  env_test <- "CRAN"
  set.seed(1351)
}else{
  # If on local machine
  env_test <- 'local'
}

test_that("test calculate_interactions function", {

  # from mgcv example
  n <- 100
  dat <- gamSim(2, n=n)
  dat <- dat$data
  dat$y <- round(exp(dat$y))
  b <- gam(y~s(x, z, bs = 'gKRLS'),
          family=poisson, data=dat, method="REML")
  
  fit_inter <- calculate_interactions(
    b, variables = list(c("x", "z")),
    continuous_type = 'onesd'
  )
  fit_main <- calculate_effects(b, variables = c('x', 'z', list(c('x', 'z'))), 
    individual = TRUE, continuous_type = 'onesd')
  
  range_x <- mean(dat$x) + c(-1, 1) * sd(dat$x)
  range_z <- mean(dat$z) + c(-1, 1) * sd(dat$z)
  avg_x <- sapply(range_x, FUN=function(i){
    mean(predict(b, newdata = data.frame(x = i, z = dat$z), type = 'response'))
  })
  avg_z <- sapply(range_z, FUN=function(i){
    mean(predict(b, newdata = data.frame(z = i, x = dat$x), type = 'response'))
  })

  expect_equal(fit_main$est[1], diff(avg_x))
  expect_equal(fit_main$est[2], diff(avg_z))

  grid_test <- expand.grid(x = range_x, z = range_z)
  grid_test$predict <- predict(b, grid_test, type = 'response')
  expect_equal(grid_test$predict[4] - grid_test$predict[1], fit_main$est[3])
  expect_equal(fit_inter[1:3, c('est', 'se', 't')], fit_main[, c('est', 'se', 't')])
  expect_equal(fit_inter$est[4], grid_test$predict[4] - grid_test$predict[1] - diff(avg_x) - diff(avg_z))
  
})

test_that("test calculate_interactions function for complex families", {
  
  # from mgcv example
  n <- 100
  dat <- gamSim(eg = 2, n = n, scale = 0.05)
  dat <- dat$data
  dat$w <- rexp(nrow(dat))
  dat$y <- dat$y * dat$w + dat$x * dat$z
  dat$w <- as.logical(dat$w > median(dat$w))
  dat$y <- cut(dat$y, 4, labels = FALSE)
  dat$fake <- rnorm(nrow(dat))
  b <- gam(y ~ w + s(x, z, fake, bs = 'gKRLS'),
           family=ocat(R=4), data=dat, method="REML")
  
  fake_cond <- data.frame(fake = rnorm(2))
  grid_type <- list(x = c(0, 0.5), z = c(0, 0.25))
  
  fit_inter <- calculate_interactions(
    b, conditional = fake_cond,
    variables = list(c("x", "z"), c("x", "w")),
    continuous_type = grid_type
  )
  
  fit_main <- calculate_effects(b, 
    variables = c('x', 'z', 'w'), 
    conditional = fake_cond,
    continuous_type = grid_type,
    individual = TRUE)
  
  orig_copy_dat <- dat[, c('x', 'z', 'w', 'fake')]
  for (f in fake_cond$fake){
    copy_dat <- orig_copy_dat
    copy_dat$fake <- f
    fake_grid <- c(grid_type, list(w = c(FALSE, TRUE)))
    est_AME <- do.call('rbind', lapply(names(fake_grid), FUN=function(i){
      copy_dat[[i]] <- fake_grid[[i]][1]
      copy_low <- colMeans(predict(b, newdata = copy_dat, type = 'response'))
      copy_dat[[i]] <- fake_grid[[i]][2]
      copy_high <- colMeans(predict(b, newdata = copy_dat, type = 'response'))
      return(data.frame(variable = i, AME = copy_high - copy_low, response = 1:length(copy_high)))
    }))
    expect_equal(est_AME$AME, fit_main$est[fit_main$fake == f])
    
    est_ACE <- do.call('rbind', lapply(list(c('x', 'z'), c('x', 'w')), FUN=function(i){
      copy_dat[[i[1]]] <- fake_grid[[i[1]]][1]
      copy_dat[[i[2]]] <- fake_grid[[i[2]]][1]
      copy_low <- colMeans(predict(b, newdata = copy_dat, type = 'response'))
      copy_dat[[i[1]]] <- fake_grid[[i[1]]][2]
      copy_dat[[i[2]]] <- fake_grid[[i[2]]][2]
      copy_high <- colMeans(predict(b, newdata = copy_dat, type = 'response'))
      return(data.frame(variable = paste(i, collapse=':'), AME = copy_high - copy_low, response = 1:length(copy_high)))
    }))
    expect_equal(est_ACE$AME, subset(fit_inter, QOI == 'ACE' & (fake == f))$est)
    
    combine_manual <- merge(
      subset(est_ACE, variable == 'x:w')[,-1],
      merge(
        subset(est_AME, variable == 'x')[,-1], 
        subset(est_AME, variable == 'w')[,-1], by = 'response'),
      by = 'response'
    )
    combine_manual$AMIE <- combine_manual$AME - combine_manual$AME.x - combine_manual$AME.y
    expect_equal(combine_manual$AMIE, subset(fit_inter, fake == f & variable == 'x:w' & QOI == 'AMIE')$est)
    
    combine_manual <- merge(
      subset(est_ACE, variable == 'x:z')[,-1],
      merge(
        subset(est_AME, variable == 'x')[,-1], 
        subset(est_AME, variable == 'z')[,-1], by = 'response'),
      by = 'response'
    )
    combine_manual$AMIE <- combine_manual$AME - combine_manual$AME.x - combine_manual$AME.y
    expect_equal(combine_manual$AMIE, subset(fit_inter, fake == f & variable == 'x:z' & QOI == 'AMIE')$est)
  }
  
})
