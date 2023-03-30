# License/copyright notice: internal_ocat_jacobian(), internal_ziP_jacobian(),
# and internal_multinom_jacobian() are slight modifications of the corresponding
# ocat()$predict, ziP()$predict, multinom()$predict from mgcv version 1.8-42
# (see <https://CRAN.R-project.org/package=mgcv>)
#
# The copyright statement for mgcv version 1.8-42 is:
# Copyright (C) 2003-2023 Simon Wood (see
# <https://CRAN.R-project.org/package=mgcv>).
#
# The license of mgcv version 1.8-42 is:
# "GPL (>=2)" (see <https://CRAN.R-project.org/package=mgcv>).

#' Internal Functions for Predict
#'
#' These functions are taken from "mgcv" with slight modifications, mostly to
#' export the derivative of the predicted value with respect to the linear
#' predictor so they can be used to calculate the jacobian for the delta method.
#' Please see \code{mgcv} for more details on the argument.
#' 
#' @param family Name of family
#' @param se Calculate standard errors
#' @param eta Not used internally for \code{gKRLS}
#' @param y Not used internally for \code{gKRLS}
#' @param X Design matrix
#' @param beta Coefficients
#' @param off Not used internally for \code{gKRLS}
#' @param Vb Variance matrix
#' @name internal_functions
#' @keywords internal
internal_ocat_jacobian <- function (family, se = FALSE, eta = NULL, y = NULL, X = NULL, 
                                    beta = NULL, off = NULL, Vb = NULL) 
{
  ocat.prob <- function(theta, lp, se = NULL) {
    R <- length(theta)
    dp <- prob <- matrix(0, length(lp), R + 2)
    prob[, R + 2] <- 1
    for (i in 1:R) {
      x <- theta[i] - lp
      ind <- x > 0
      prob[ind, i + 1] <- 1/(1 + exp(-x[ind]))
      ex <- exp(x[!ind])
      prob[!ind, i + 1] <- ex/(1 + ex)
      dp[, i + 1] <- prob[, i + 1] * (prob[, i + 1] - 1)
    }
    prob <- t(diff(t(prob)))
    dp <- t(diff(t(dp)))
    if (!is.null(se)) 
      se <- as.numeric(se) * abs(dp)
    list(prob, se, dp)
  }
  theta <- family$getTheta(TRUE)
  if (is.null(eta)) {
    discrete <- is.list(X)
    mu <- off + if (discrete) 
      Xbd(X$Xd, beta, k = X$kd, ks = X$ks, ts = X$ts, dt = X$dt, 
          v = X$v, qc = X$qc, drop = X$drop)
    else drop(X %*% beta)
    if (se) {
      se <- if (discrete) 
        sqrt(pmax(0, diagXVXd(X$Xd, Vb, k = X$kd, ks = X$ks, 
                              ts = X$ts, dt = X$dt, v = X$v, qc = X$qc, drop = X$drop, 
                              nthreads = 1)))
      else sqrt(pmax(0, rowSums((X %*% Vb) * X)))
    }
    else se <- NULL
    p <- ocat.prob(theta, mu, se)
    if (is.null(se)) 
      return(p)
    else {
      names(p) <- c("fit", "se.fit", "jacobian")
      return(p)
    }
  }
}

#'
#' @rdname internal_functions
#' @keywords internal
internal_ziP_jacobian <- function (family, se = FALSE, eta = NULL, y = NULL, X = NULL, 
                                   beta = NULL, off = NULL, Vb = NULL) 
{
  theta <- family$getTheta()
  .b <- environment(family$getTheta)$.b
  .Theta_copy <- environment(family$getTheta)$.Theta
  if (!isTRUE(identical(theta, .Theta_copy))){
    stop('Unable to extract theta from ziP environment.')
  }
  if (is.null(eta)) {
    discrete <- is.list(X)
    gamma <- off + if (discrete) 
      Xbd(X$Xd, beta, k = X$kd, ks = X$ks, ts = X$ts, dt = X$dt, 
          v = X$v, qc = X$qc, drop = X$drop)
    else drop(X %*% beta)
    if (se) {
      se <- if (discrete) 
        sqrt(pmax(0, diagXVXd(X$Xd, Vb, k = X$kd, ks = X$ks, 
                              ts = X$ts, dt = X$dt, v = X$v, qc = X$qc, drop = X$drop, 
                              nthreads = 1)))
      else sqrt(pmax(0, rowSums((X %*% Vb) * X)))
    }
    else se <- NULL
  }
  else {
    se <- NULL
    gamma <- eta
  }
  b <- get(".b")
  eta <- theta[1] + (b + exp(theta[2])) * gamma
  et <- exp(eta)
  mu <- p <- 1 - exp(-et)
  fv <- lambda <- exp(gamma)
  ind <- gamma < log(.Machine$double.eps)/2
  mu[!ind] <- lambda[!ind]/(1 - exp(-lambda[!ind]))
  mu[ind] <- 1
  fv <- list(p * mu)
  if (is.null(se)) 
    return(fv)
  else {
    dp.dg <- p
    dp.dg <- exp(-et) * et * (b + exp(theta[2]))
    dmu.dg <- (lambda + 1) * mu - mu^2
    fv[[2]] <- abs(dp.dg * mu + dmu.dg * p) * se
    fv[[3]] <- dp.dg * mu + dmu.dg * p
    names(fv) <- c("fit", "se.fit", "jacobian")
    return(fv)
  }
}

#' @rdname internal_functions
#' @keywords internal
internal_multinom_jacobian <- function (family, se = FALSE, eta = NULL, y = NULL, X = NULL, 
                                        beta = NULL, off = NULL, Vb = NULL) 
{
  if (is.null(eta)) {
    discrete <- is.list(X)
    lpi <- attr(X, "lpi")
    if (is.null(lpi)) {
      lpi <- list(1:ncol(X))
    }
    K <- length(lpi)
    nobs <- if (discrete) 
      nrow(X$kd)
    else nrow(X)
    eta <- matrix(0, nobs, K)
    if (se) {
      ve <- matrix(0, nobs, K)
      ce <- matrix(0, nobs, K * (K - 1)/2)
    }
    ii <- 0
    for (i in 1:K) {
      if (discrete) {
        eta[, i] <- Xbd(X$Xd, beta, k = X$kd, ks = X$ks, 
                        ts = X$ts, dt = X$dt, v = X$v, qc = X$qc, drop = X$drop, 
                        lt = X$lpid[[i]])
      }
      else {
        Xi <- X[, lpi[[i]], drop = FALSE]
        eta[, i] <- Xi %*% beta[lpi[[i]]]
      }
      if (!is.null(off[[i]])) 
        eta[, i] <- eta[, i] + off[[i]]
      if (se) {
        ve[, i] <- if (discrete) 
          diagXVXd(X$Xd, Vb, k = X$kd, ks = X$ks, ts = X$ts, 
                   dt = X$dt, v = X$v, qc = X$qc, drop = X$drop, 
                   nthreads = 1, lt = X$lpid[[i]], rt = X$lpid[[i]])
        else drop(pmax(0, rowSums((Xi %*% Vb[lpi[[i]], 
                                             lpi[[i]]]) * Xi)))
        # ii <- 0
        if (i < K) 
          for (j in (i + 1):K) {
            ii <- ii + 1
            ce[, ii] <- if (discrete) 
              diagXVXd(X$Xd, Vb, k = X$kd, ks = X$ks, 
                       ts = X$ts, dt = X$dt, v = X$v, qc = X$qc, 
                       drop = X$drop, nthreads = 1, lt = X$lpid[[i]], 
                       rt = X$lpid[[j]])
            else drop(pmax(0, rowSums((Xi %*% Vb[lpi[[i]], 
                                                 lpi[[j]]]) * X[, lpi[[j]]])))
          }
      }
    }
  }
  else {
    se <- FALSE
  }
  gamma <- cbind(1, exp(eta))
  beta <- rowSums(gamma)
  gamma <- gamma/beta
  vp <- gamma * 0
  store_dp <- as.list(rep(NA, K+1))
  if (se) {
    for (j in 1:(K + 1)) {
      if (j == 1) 
        dp <- -gamma[, -1, drop = FALSE]/beta
      else {
        dp <- -gamma[, j] * gamma[, -1, drop = FALSE]
        dp[, j - 1] <- gamma[, j] * (1 - gamma[, j])
      }
      store_dp[[j]] <- dp
      vp[, j] <- rowSums(dp^2 * ve)
      ii <- 0
      for (i in 1:K) if (i < K) 
        for (k in (i + 1):K) {
          ii <- ii + 1
          vp[, j] <- vp[, j] + 2 * dp[, i] * dp[, k] * ce[, ii]
        }
      vp[, j] <- sqrt(pmax(0, vp[, j]))
    }
    return(list(fit = gamma, se.fit = vp, ce = ce,
                gamma = gamma, beta = beta,
                jacobian = store_dp, ve = ve))
  }
  list(fit = gamma)
}