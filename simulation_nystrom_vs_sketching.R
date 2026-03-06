## Simulation: Sketching vs Nystrom Approximation
##
## Compares subsampling sketching and Nystrom approximation across:
##   1. Coefficient recovery / model fit (in-sample R^2)
##   2. Out-of-sample prediction (RMSE)
##   3. Marginal effects (AME estimates and SEs)
##   4. Interactions
##   5. Robust standard errors
##   6. Binary outcome (logistic)
##   7. Bandwidth calibration
##
## If both methods approximate the full kernel well, results should
## be broadly similar.

devtools::load_all(".")
library(sandwich)

set.seed(2026)

# ---------------------------------------------------------------------------
# Helper: root mean squared error
# ---------------------------------------------------------------------------
rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2))

# ---------------------------------------------------------------------------
# 1. Data-generating process
# ---------------------------------------------------------------------------
N <- 800
N_test <- 200
P <- 4

X <- matrix(rnorm(N * P), ncol = P)
colnames(X) <- paste0("x", 1:P)

# Non-linear truth: sin + interaction + quadratic
f_true <- function(x) {
  sin(2 * x[, 1]) + 0.6 * x[, 2] * x[, 3] + 0.3 * x[, 4]^2
}

y <- f_true(X) + rnorm(N, sd = 0.5)
dat <- data.frame(y = y, X)

# Hold-out test set
X_test <- matrix(rnorm(N_test * P), ncol = P)
colnames(X_test) <- paste0("x", 1:P)
y_test <- f_true(X_test) + rnorm(N_test, sd = 0.5)
dat_test <- data.frame(y = y_test, X_test)

cat("=== Data generated: N =", N, ", N_test =", N_test, ", P =", P, "===\n\n")

# ---------------------------------------------------------------------------
# 2. Fit models
# ---------------------------------------------------------------------------
formula <- y ~ s(x1, x2, x3, x4, bs = "gKRLS", xt = gKRLS(
  sketch_method = METHOD, sketch_multiplier = 5
))

fit_models <- function() {
  cat("--- Fitting subsampling model ---\n")
  t_sub <- system.time(
    fit_sub <- gam(y ~ s(x1, x2, x3, x4, bs = "gKRLS",
      xt = gKRLS(sketch_method = "subsampling", sketch_multiplier = 5)),
      data = dat)
  )

  cat("--- Fitting nystrom model ---\n")
  t_nys <- system.time(
    fit_nys <- gam(y ~ s(x1, x2, x3, x4, bs = "gKRLS",
      xt = gKRLS(sketch_method = "nystrom", sketch_multiplier = 5)),
      data = dat)
  )

  cat("--- Fitting no-sketch (full kernel) model ---\n")
  t_full <- system.time(
    fit_full <- gam(y ~ s(x1, x2, x3, x4, bs = "gKRLS",
      xt = gKRLS(sketch_method = "none")),
      data = dat)
  )

  list(
    sub  = list(fit = fit_sub, time = t_sub["elapsed"]),
    nys  = list(fit = fit_nys, time = t_nys["elapsed"]),
    full = list(fit = fit_full, time = t_full["elapsed"])
  )
}

models <- fit_models()

cat("\n====================================================\n")
cat("SECTION 1: Model Fit & Timing\n")
cat("====================================================\n\n")

for (nm in names(models)) {
  m <- models[[nm]]
  r2 <- 1 - sum(residuals(m$fit)^2) / sum((dat$y - mean(dat$y))^2)
  cat(sprintf("  %-12s  R2 = %.4f   edf = %6.1f   time = %.2fs\n",
    nm, r2, sum(summary(m$fit)$s.table[, "edf"]), m$time))
}

# ---------------------------------------------------------------------------
# 3. Out-of-sample prediction
# ---------------------------------------------------------------------------
cat("\n====================================================\n")
cat("SECTION 2: Out-of-Sample Prediction (RMSE)\n")
cat("====================================================\n\n")

for (nm in names(models)) {
  pred <- predict(models[[nm]]$fit, newdata = dat_test, se.fit = TRUE)
  cat(sprintf("  %-12s  RMSE = %.4f   mean(SE) = %.4f\n",
    nm, rmse(dat_test$y, pred$fit), mean(pred$se.fit)))
}

# Correlation of predictions between methods
pred_sub  <- predict(models$sub$fit, newdata = dat_test)
pred_nys  <- predict(models$nys$fit, newdata = dat_test)
pred_full <- predict(models$full$fit, newdata = dat_test)

cat(sprintf("\n  Prediction correlation (sub vs nys):   %.6f\n", cor(pred_sub, pred_nys)))
cat(sprintf("  Prediction correlation (sub vs full):  %.6f\n", cor(pred_sub, pred_full)))
cat(sprintf("  Prediction correlation (nys vs full):  %.6f\n", cor(pred_nys, pred_full)))

# ---------------------------------------------------------------------------
# 4. Marginal effects (AME)
# ---------------------------------------------------------------------------
cat("\n====================================================\n")
cat("SECTION 3: Average Marginal Effects\n")
cat("====================================================\n\n")

vars <- paste0("x", 1:P)

mfx_sub  <- calculate_effects(models$sub$fit,  variables = vars)
mfx_nys  <- calculate_effects(models$nys$fit,  variables = vars)
mfx_full <- calculate_effects(models$full$fit, variables = vars)

mfx_all <- merge(
  merge(
    data.frame(variable = mfx_sub$variable, type = mfx_sub$type,
      est_sub = mfx_sub$est, se_sub = mfx_sub$se),
    data.frame(variable = mfx_nys$variable, type = mfx_nys$type,
      est_nys = mfx_nys$est, se_nys = mfx_nys$se),
    by = c("variable", "type")
  ),
  data.frame(variable = mfx_full$variable, type = mfx_full$type,
    est_full = mfx_full$est, se_full = mfx_full$se),
  by = c("variable", "type")
)

cat("  AME Estimates:\n")
print(mfx_all[, c("variable", "type", "est_sub", "est_nys", "est_full")], row.names = FALSE)
cat("\n  AME Standard Errors:\n")
print(mfx_all[, c("variable", "type", "se_sub", "se_nys", "se_full")], row.names = FALSE)

cat(sprintf("\n  Max |est_sub - est_nys|:  %.6f\n", max(abs(mfx_all$est_sub - mfx_all$est_nys))))
cat(sprintf("  Max |est_sub - est_full|: %.6f\n", max(abs(mfx_all$est_sub - mfx_all$est_full))))
cat(sprintf("  Max |est_nys - est_full|: %.6f\n", max(abs(mfx_all$est_nys - mfx_all$est_full))))

# ---------------------------------------------------------------------------
# 5. Individual-level marginal effects
# ---------------------------------------------------------------------------
cat("\n====================================================\n")
cat("SECTION 4: Individual Marginal Effects (first 5 obs)\n")
cat("====================================================\n\n")

ime_sub  <- calculate_effects(models$sub$fit,  variables = "x1", individual = TRUE)
ime_nys  <- calculate_effects(models$nys$fit,  variables = "x1", individual = TRUE)
ime_full <- calculate_effects(models$full$fit, variables = "x1", individual = TRUE)

ind_sub  <- attr(ime_sub, "individual")$est
ind_nys  <- attr(ime_nys, "individual")$est
ind_full <- attr(ime_full, "individual")$est

cat("  Correlation of individual effects (x1):\n")
cat(sprintf("    sub vs nys:  %.6f\n", cor(ind_sub, ind_nys)))
cat(sprintf("    sub vs full: %.6f\n", cor(ind_sub, ind_full)))
cat(sprintf("    nys vs full: %.6f\n", cor(ind_nys, ind_full)))

# ---------------------------------------------------------------------------
# 6. Interaction effects
# ---------------------------------------------------------------------------
cat("\n====================================================\n")
cat("SECTION 5: Interaction Effects (x2:x3)\n")
cat("====================================================\n\n")

int_sub  <- calculate_interactions(models$sub$fit,  variables = list(c("x2", "x3")))
int_nys  <- calculate_interactions(models$nys$fit,  variables = list(c("x2", "x3")))
int_full <- calculate_interactions(models$full$fit, variables = list(c("x2", "x3")))

int_compare <- data.frame(
  QOI = int_sub$QOI,
  variable = int_sub$variable,
  est_sub  = int_sub$est,
  est_nys  = int_nys$est,
  est_full = int_full$est
)
print(int_compare, row.names = FALSE)

# ---------------------------------------------------------------------------
# 7. Robust / clustered standard errors
# ---------------------------------------------------------------------------
cat("\n====================================================\n")
cat("SECTION 6: Robust Standard Errors (HC1)\n")
cat("====================================================\n\n")

vcov_sub  <- vcovHC(models$sub$fit,  type = "HC1")
vcov_nys  <- vcovHC(models$nys$fit,  type = "HC1")
vcov_full <- vcovHC(models$full$fit, type = "HC1")

# Compare diagonal elements (variances of kernel coefficients)
se_sub  <- sqrt(diag(vcov_sub))
se_nys  <- sqrt(diag(vcov_nys))
se_full <- sqrt(diag(vcov_full))

# The number of coefficients may differ, so compare the intercept
cat(sprintf("  Intercept robust SE:  sub = %.6f,  nys = %.6f,  full = %.6f\n",
  se_sub[1], se_nys[1], se_full[1]))
cat(sprintf("  Mean kernel coef SE:  sub = %.6f,  nys = %.6f,  full = %.6f\n",
  mean(se_sub[-1]), mean(se_nys[-1]), mean(se_full[-1])))

# Marginal effects with robust SEs
cat("\n  AME with robust SEs (x1):\n")
mfx_r_sub  <- calculate_effects(models$sub$fit,  variables = "x1", vcov = vcov_sub)
mfx_r_nys  <- calculate_effects(models$nys$fit,  variables = "x1", vcov = vcov_nys)
mfx_r_full <- calculate_effects(models$full$fit, variables = "x1", vcov = vcov_full)

cat(sprintf("    sub:  est = %.4f, robust_se = %.4f\n", mfx_r_sub$est, mfx_r_sub$se))
cat(sprintf("    nys:  est = %.4f, robust_se = %.4f\n", mfx_r_nys$est, mfx_r_nys$se))
cat(sprintf("    full: est = %.4f, robust_se = %.4f\n", mfx_r_full$est, mfx_r_full$se))

# ---------------------------------------------------------------------------
# 8. Binary outcome (logistic regression)
# ---------------------------------------------------------------------------
cat("\n====================================================\n")
cat("SECTION 7: Binary Outcome (Logistic)\n")
cat("====================================================\n\n")

prob <- plogis(f_true(X))
y_bin <- rbinom(N, 1, prob)
dat_bin <- data.frame(y = y_bin, X)

prob_test <- plogis(f_true(X_test))
y_bin_test <- rbinom(N_test, 1, prob_test)
dat_bin_test <- data.frame(y = y_bin_test, X_test)

cat("--- Fitting logistic models ---\n")
fit_bin_sub <- gam(y ~ s(x1, x2, x3, x4, bs = "gKRLS",
    xt = gKRLS(sketch_method = "subsampling", sketch_multiplier = 5)),
  data = dat_bin, family = binomial(link = "logit"))

fit_bin_nys <- gam(y ~ s(x1, x2, x3, x4, bs = "gKRLS",
    xt = gKRLS(sketch_method = "nystrom", sketch_multiplier = 5)),
  data = dat_bin, family = binomial(link = "logit"))

# Prediction
p_sub  <- predict(fit_bin_sub, newdata = dat_bin_test, type = "response")
p_nys  <- predict(fit_bin_nys, newdata = dat_bin_test, type = "response")

cat(sprintf("  Prediction correlation (sub vs nys): %.6f\n", cor(p_sub, p_nys)))
cat(sprintf("  Brier score sub:  %.6f\n", mean((dat_bin_test$y - p_sub)^2)))
cat(sprintf("  Brier score nys:  %.6f\n", mean((dat_bin_test$y - p_nys)^2)))

# Marginal effects (binary)
mfx_bin_sub <- calculate_effects(fit_bin_sub, variables = "x1")
mfx_bin_nys <- calculate_effects(fit_bin_nys, variables = "x1")

cat(sprintf("\n  AME x1 (logistic):  sub = %.4f (se %.4f),  nys = %.4f (se %.4f)\n",
  mfx_bin_sub$est, mfx_bin_sub$se, mfx_bin_nys$est, mfx_bin_nys$se))

# ---------------------------------------------------------------------------
# 9. Bandwidth calibration
# ---------------------------------------------------------------------------
cat("\n====================================================\n")
cat("SECTION 8: Bandwidth Calibration\n")
cat("====================================================\n\n")

fit_cal_sub <- gam(y ~ s(x1, x2, x3, x4, bs = "gKRLS",
    xt = gKRLS(sketch_method = "subsampling", bandwidth = "calibrate",
      sketch_multiplier = 5)),
  data = dat)

fit_cal_nys <- gam(y ~ s(x1, x2, x3, x4, bs = "gKRLS",
    xt = gKRLS(sketch_method = "nystrom", bandwidth = "calibrate",
      sketch_multiplier = 5)),
  data = dat)

cal_sub <- get_calibration_information(fit_cal_sub)
cal_nys <- get_calibration_information(fit_cal_nys)

cat(sprintf("  Calibrated bandwidth:  sub = %.4f,  nys = %.4f  (default = %d)\n",
  cal_sub$bandwidth, cal_nys$bandwidth, P))

# Compare fit with calibrated bandwidth
r2_cal_sub <- 1 - sum(residuals(fit_cal_sub)^2) / sum((dat$y - mean(dat$y))^2)
r2_cal_nys <- 1 - sum(residuals(fit_cal_nys)^2) / sum((dat$y - mean(dat$y))^2)
pred_cal_sub <- predict(fit_cal_sub, newdata = dat_test)
pred_cal_nys <- predict(fit_cal_nys, newdata = dat_test)

cat(sprintf("  Calibrated R2:         sub = %.4f,  nys = %.4f\n", r2_cal_sub, r2_cal_nys))
cat(sprintf("  Calibrated RMSE:       sub = %.4f,  nys = %.4f\n",
  rmse(dat_test$y, pred_cal_sub), rmse(dat_test$y, pred_cal_nys)))

# ---------------------------------------------------------------------------
# 10. Monte Carlo: repeated draws to assess variability
# ---------------------------------------------------------------------------
cat("\n====================================================\n")
cat("SECTION 9: Monte Carlo (10 replications)\n")
cat("====================================================\n\n")

n_rep <- 10
mc_results <- data.frame(
  rep = integer(), method = character(),
  r2 = numeric(), rmse_oos = numeric(),
  ame_x1 = numeric(), se_x1 = numeric(),
  stringsAsFactors = FALSE
)

for (r in seq_len(n_rep)) {
  set.seed(2025 + r)

  # Resample y with new noise
  y_r <- f_true(X) + rnorm(N, sd = 0.5)
  dat_r <- data.frame(y = y_r, X)

  for (method in c("subsampling", "nystrom")) {
    fit_r <- gam(y ~ s(x1, x2, x3, x4, bs = "gKRLS",
        xt = gKRLS(sketch_method = method, sketch_multiplier = 5)),
      data = dat_r)

    r2_r <- 1 - sum(residuals(fit_r)^2) / sum((dat_r$y - mean(dat_r$y))^2)
    pred_r <- predict(fit_r, newdata = dat_test)
    rmse_r <- rmse(dat_test$y, pred_r)
    mfx_r <- calculate_effects(fit_r, variables = "x1")

    mc_results <- rbind(mc_results, data.frame(
      rep = r, method = method,
      r2 = r2_r, rmse_oos = rmse_r,
      ame_x1 = mfx_r$est, se_x1 = mfx_r$se,
      stringsAsFactors = FALSE
    ))
  }
  cat(sprintf("  Replication %d/%d complete\n", r, n_rep))
}

cat("\n  Summary across replications:\n\n")
agg <- aggregate(cbind(r2, rmse_oos, ame_x1, se_x1) ~ method, data = mc_results,
  FUN = function(x) c(mean = mean(x), sd = sd(x)))

for (m in c("subsampling", "nystrom")) {
  row <- agg[agg$method == m, ]
  cat(sprintf("  %-12s  R2 = %.4f (%.4f)  RMSE = %.4f (%.4f)  AME_x1 = %.4f (%.4f)  SE_x1 = %.4f (%.4f)\n",
    m,
    row$r2[1], row$r2[2],
    row$rmse_oos[1], row$rmse_oos[2],
    row$ame_x1[1], row$ame_x1[2],
    row$se_x1[1], row$se_x1[2]))
}

cat("\n====================================================\n")
cat("SIMULATION COMPLETE\n")
cat("====================================================\n")
