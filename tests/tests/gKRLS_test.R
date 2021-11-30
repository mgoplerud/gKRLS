rm(list=ls())
library(testthat)
library(gKRLS)

X <- matrix(rnorm(500 * 2), ncol = 2)
X <- apply(X, MARGIN = 2, scale)
X <- cbind(X, model.matrix(~0 + sample(letters[1:7], nrow(X), replace = T)))
w <- runif(nrow(X))
colnames(X) <- paste0('x_', 1:ncol(X))
y <- rnorm(nrow(X), mean = X %*% rnorm(ncol(X), sd = 5))


est_gKRLS <- gKRLS(formula = y ~ 1 + w, kernel_X =  X, 
                   standardize = 'scaled',
                   sketch_method = 'gaussian', data = data.frame(w = w),
                   family = gaussian())

pred_gKRLS <- predict(est_gKRLS, newkernel_X = X, newdata = data.frame(w = w))


fitted_internal <- fitted(est_gKRLS)
expect_equivalent(pred_gKRLS$fitted, fitted_internal)


MFX_SIZE <- sample(1:nrow(X), nrow(X), replace = T)

mfx_gKRLS_base <- marginal_effect_base(object = est_gKRLS,
    newdata = data.frame(w = w[MFX_SIZE]), 
    newkernel_X = X[MFX_SIZE,], method = 'base')

mfx_gKRLS <- marginal_effect_base(object = est_gKRLS,
    newdata = data.frame(w = w[MFX_SIZE]), newkernel_X = X[MFX_SIZE,])

# mfx calculated both ways aligns
expect_true(all(sapply(names(mfx_gKRLS), FUN=function(i){
  all.equal(mfx_gKRLS[[i]], mfx_gKRLS_base[[i]])
})))


mfx_gKRLS_base_select <- marginal_effect_base(object = est_gKRLS, keep = c('w', 'x_1', 'x_6', 'x_2'),
                                       newdata = data.frame(w = w[MFX_SIZE]), 
                                       newkernel_X = X[MFX_SIZE,], method = 'base')

mfx_gKRLS_select <- marginal_effect_base(object = est_gKRLS, keep = c('w', 'x_1', 'x_6', 'x_2'),
        newdata = data.frame(w = w[MFX_SIZE]), newkernel_X = X[MFX_SIZE,])

expect_true(all(sapply(names(mfx_gKRLS), FUN=function(i){
  all.equal(mfx_gKRLS_select[[i]], mfx_gKRLS_base_select[[i]])
})))
# expect_identical(colnames(mfx_gKRLS_select$ME_pointwise), c('w', 'x_1', 'x_6', 'x_2'))

expect_equal(
  mfx_gKRLS$ME_pointwise_var[,colnames(mfx_gKRLS_select$ME_pointwise_var)] ,
  mfx_gKRLS_select$ME_pointwise_var
)



stop()



stop()

N <- 200
x1 <- rnorm(N)
x2 <- rbinom(N,size=1,prob=.2)
y <- x1^3 - 0.5 *x2 + rnorm(N,0,.15)
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
      family = gaussian(), prior_stabilize = F)
fit_krls <- KRLS::krls(X = X, y = y, sigma = fit_gKRLS$internal$bandwidth, 
                 lambda = 1/fit_gKRLS$fmt_varcorr[[1]][1,1] * fit_gKRLS$sigma^2)

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

