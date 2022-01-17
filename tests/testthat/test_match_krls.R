
test_that("gKRLS agrees with direct solution", {

  N <- 200
  x1 <- rnorm(N)
  x2 <- rbinom(N,size=1,prob=.2)
  y <- x1^3 - 0.5 *x2 + rnorm(N,0,.15)
  y <- y * 10
  X <- cbind(x1, x2)
  X <- cbind(X, model.matrix(~ 0 + sample(letters[1:5], N, replace  = T)))
  X_copy <- cbind(X, X, X)
  
  est_gKRLS <- gKRLS(formula = y ~ 1, kernel_X = X, 
                     sketch_method = 'none', data = NULL,
                     family = gaussian())
  
  est_copy <- gKRLS(formula = y ~ 1, kernel_X = X_copy, data = NULL,
                    sketch_method = 'none',
                    family = gaussian())
  
  expect_equal(est_gKRLS$fitted, est_copy$fitted)
  
  fit_gKRLS <- gKRLS(formula = y ~ 0, kernel_X = X, 
                     sketch_method = 'none', data = NULL, 
                     truncate_eigen_penalty = .Machine$double.eps,
                     family = gaussian(), prior_stabilize = F)
  lambda <- fit_gKRLS$internal$KS_scale^2 * 
    1/fit_gKRLS$fmt_varcorr[[1]][1,1] * 
    fit_gKRLS$sigma^2
  K_manual <- exp(-as.matrix(dist(apply(X, MARGIN = 2, scale))^2) / fit_gKRLS$internal$bandwidth)
  
  direct_fitted <- as.vector(K_manual %*% 
                               solve(K_manual + lambda * Diagonal(n = nrow(K_manual))) %*% 
                               scale(y)) * sd(y) + mean(y)  
  fit_krls <- KRLS::krls(X = X, y = y, sigma = fit_gKRLS$internal$bandwidth, 
                         lambda = lambda)
  plot(as.vector(fitted(fit_krls)) - as.vector(fitted(fit_gKRLS)))
  
  expect_equivalent(fitted(fit_krls), fitted(fit_gKRLS))
  expect_equivalent(direct_fitted, fitted(fit_gKRLS))
  expect_equal(direct_fitted, as.vector(fitted(fit_krls)))
  
  
})


fit_krls$K

fake_K <- function(object, raw_X){
  
  # Standardize the incoming new data.
  std_newkernel_X <- sweep(raw_X, 2, object$internal$std_train$mean, FUN = "-")
  
  std_newkernel_X <- std_newkernel_X %*% 
    object$internal$std_train$whiten
  std_newkernel_X <- as.matrix(std_newkernel_X)  
  # Standardize the saved training kernel
  std_kernel_X <- object$internal$kernel_X_train
  std_kernel_X <- sweep(std_kernel_X, 2, object$internal$std_train$mean, FUN = "-")
  std_kernel_X <- std_kernel_X %*% 
    object$internal$std_train$whiten
  std_kernel_X <- as.matrix(std_kernel_X)  
  
  # Get the "test" data after using same sketch matrix
  newdataKS <- create_sketched_kernel(X_test = raw_X, 
                                      X_train = raw_X, 
                                      tS = t(object$internal$sketch), 
                                      bandwidth = object$internal$bandwidth)
  
}
plot(fitted(fit_gKRLS), fitted(fit_krls))
abline(a=0,b=1)

mfx_gKRLS_base <- marginal_effect_base(object = fit_gKRLS, newdata = NULL, newkernel_X = X, method = 'base')
mfx_gKRLS <- marginal_effect_base(object = fit_gKRLS, newdata = NULL, newkernel_X = X)
mfx_gKRLS$AME_pointwise
mfx_gKRLS_base$AME_pointwise
mfx_gKRLS$AME_pointwise - fit_krls$avgderivatives

expect_equal(as.vector(fit_krls$derivatives), 
             as.vector(mfx_gKRLS_base$ME_pointwise), tol = 1e-5)


expect_equal(as.vector(mfx_gKRLS$ME_pointwise), 
             as.vector(mfx_gKRLS_base$ME_pointwise), tol = 1e-5)


fit_with_FE <- gKRLS(formula = y ~ 1 + w + factor(id), kernel_X =  X, 
      standardize = 'scaled',
      sketch_method = 'gaussian', data = data.frame(w = w, id = sample(letters, length(w), replace = T)),
      family = gaussian())

predict(fit_with_FE, newdata = data.frame(w = rnorm(3), id = letters[c(3,5,10)]),
        newkernel_X = X[c(14,5,2),])


mfx_gKRLS_base <- marginal_effect_base(object = fit_with_FE, keep = c('x_1', 'w'),
                                       newdata = data.frame(w = rnorm(3), id = letters[c(3,5,10)]), 
                                       newkernel_X = X[c(3,10,6),])

stop()



N <- 200
x1 <- rnorm(N)
x2 <- rbinom(N,size=1,prob=.2)
x2 <- x2 + sample(c(1,0,-1), size = N, replace = T) * 1e-6
b1 <- 1
b2 <- -3
y <- b1 * x1^3 + b2 * x2 + rnorm(N,0,.15)
X <- cbind(x1, x2)

bin_y <- rbinom(nrow(X), 1, plogis(scale(y)))

fit_binary_gKRLS <- gKRLS(formula = bin_y ~ 1, kernel_X = X, 
                   sketch_method = 'none', data = NULL, bandwidth = 2,
                   family = binomial(), prior_stabilize = T)

fit_krls_logit <- KRLS2::krls(X = X, y = bin_y, b = 2, 
            loss = 'logistic', truncate = 1 - sqrt(.Machine$double.eps),
            lambda = 1/2 * 1/nrow(X) * 1/fit_binary_gKRLS$fmt_varcorr[[1]][1,1])

plot(fitted(fit_krls_logit), fitted(fit_binary_gKRLS))

mfx_logit <- marginal_effect_base(fit_binary_gKRLS, newdata = NULL, newkernel_X = X)
mfx_KRLS_2 <- KRLS2::summary.krls2(fit_krls_logit)

base_glm <- glm(bin_y ~ I(X[,1]^3) + X[,2], family = binomial())
coef(base_glm)
plot(as.vector(mfx_KRLS_2$derivatives), as.vector(mfx_logit$ME_pointwise[,-1]))

cbind(mfx_logit$AME_pointwise[-1], mfx_KRLS_2$avgderivatives)

ggplot() + geom_density(aes(x=value, group = Var2), data = reshape2::melt(mfx_KRLS_2$derivatives)) +
  facet_wrap(~Var2) +
  geom_density(aes(x=value,group=Var2), data = reshape2::melt(mfx_logit$ME_pointwise[,-1]), col = 'red') 


# Linear Case
fit_extend_gKRLS <- gKRLS(formula = y ~ 1, kernel_X = X[,1,drop=F], bandwidth = 10^6,
                   sketch_method = 'none', data = data.frame(X),
                   family = gaussian(), prior_stabilize = T)
extrapolate_X <- as.matrix(expand.grid(seq(3 * min(X[,1]), 3 * max(X[,2]), length.out = 1000), c(0)))
colnames(extrapolate_X) <- c('x1','x2')
plot(extrapolate_X[,1] , predict(fit_extend_gKRLS, newdata = data.frame(as.matrix(extrapolate_X)), newkernel_X = extrapolate_X[,1,drop=F])$fitted)
abline(lm(y ~ X[,1]), col = 'red')

      
g <- sample(1:5, N, replace = T)
verify_glmer <- gKRLS(formula = bin_y ~ 1 + (1 | g), kernel_X = X[,1,drop=F], bandwidth = 10^6,
                          sketch_method = 'none', data = data.frame(X),
                          family = binomial(), prior_stabilize = T)
basic_glmer <- glmer(bin_y ~ 1 + (1 | g), family = binomial())

#


gKRLS_fixed <- gKRLS(formula = out_y ~ 1 + x1 + x2,
                     kernel_X = kern_data,  prior_stabilize = F,
                     data = reg_data %>% dplyr::mutate(out_y = (100 * as.numeric(outcome_y > 0))), 
                     family = gaussian())

gKRLS_fixed2 <- gKRLS(formula = out_y ~ 1 + x1 + x2, standardize_outcome = F,
                      kernel_X = kern_data,  prior_stabilize = F,
                      data = reg_data %>% dplyr::mutate(out_y = (100 * as.numeric(outcome_y > 0))), 
                      family = gaussian())

b <- glm(out_y ~ 1 + x1 + x2,
         data = reg_data %>% dplyr::mutate(out_y = (100 * as.numeric(outcome_y > 0))),
         family = gaussian())

