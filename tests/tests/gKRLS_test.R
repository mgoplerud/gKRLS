library(testthat)
est_gKRLS <- gKRLS(formula = y ~ 1, kernel_X =  X, standardize = 'Mahalanobis',
                   sketch_method = 'gaussian',
                   family = gaussian())

pred_gKRLS <- predict(est_gKRLS, newdata = X)


fitted_internal <- fitted(est_gKRLS)
manual_predict <- predict(est_gKRLS, newdata = X)
expect_equivalent(manual_predict$fitted, fitted_internal)


N <- 200
x1 <- rnorm(N)
x2 <- rbinom(N,size=1,prob=.2)
y <- x1^3 - 0.5 *x2 + rnorm(N,0,.15)
X <- cbind(x1, x2)
X_copy <- cbind(X, X, X)

est_gKRLS <- gKRLS(formula = y ~ 1, kernel_X = X, 
                   sketch_method = 'none',
                   family = gaussian())

est_copy <- gKRLS(y = y, X = X_copy, 
                   sketch_method = 'none',
                   family = gaussian())

expect_equal(est_gKRLS$fitted, est_copy$fitted)

plot(predict(est_copy, newdata = X_copy)$fitted, est_copy$fitted)



fit_gKRLS <- gKRLS(formula = y ~ 0, kernel_X = X, 
      sketch_method = 'none', data = NULL,
      family = gaussian(), prior_stabilize = F)
fit_krls <- krls(X = X, y = y, sigma = fit_gKRLS$internal$bandwidth, 
                 lambda = 1/fit_gKRLS$fmt_varcorr[[1]][1,1] * fit_gKRLS$sigma^2)

plot(fitted(fit_gKRLS), fitted(fit_krls))
abline(a=0,b=1)

mfx_gKRLS <- marginal_effect_base(object = fit_gKRLS, newdata = NULL, newkernel_X = X)
mfx_gKRLS_base <- marginal_effect_base(object = fit_gKRLS, newdata = NULL, newkernel_X = X, method = 'base')
mfx_gKRLS$AME_pointwise
mfx_gKRLS_base$AME_pointwise
mfx_gKRLS$AME_pointwise - fit_krls$avgderivatives

plot(as.vector(fit_krls$derivatives), as.vector(mfx_gKRLS$ME_pointwise))

plot(fit_krls$coeffs, fit_gKRLS$)


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

