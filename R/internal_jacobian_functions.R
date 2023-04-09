# Functions to calculate the relevant derivative of the
# outcomes with respect to the linear predictor for extended
# families that have custom predict functions

ocat_jacobian <- function(lp_i, family_object){
  cumprob <- cbind(0, 
                   outer(lp_i, family_object$getTheta(TRUE), 
                         FUN=function(l, th){plogis(th-l)}),
                   1)
  dprob <-  cumprob^2 - cumprob # - p * (1 - p)
  fit <- t(diff(t(cumprob)))
  jacobian <- t(diff(t(dprob)))
  return(
    list(fit = fit, jacobian = jacobian)
  )
}

zip_jacobian <- function(lp_i, family_object){
  transf_theta <- family_object$getTheta(trans = TRUE)
  lambda <- exp(lp_i)
  eta <- transf_theta[1] + transf_theta[2] * lp_i
  # Get the probabiltiy of non-zero outcome
  p_nonzero <- 1 - exp(-exp(eta))
  # Get expected outcome if non-zero, i.e. E[Y | X > 0]
  E_nonzero <- lambda/(1-exp(-lambda))
  fit <- p_nonzero * E_nonzero
  # Get the derivative of fit w.r.t. lp
  d1 <- exp(eta - exp(eta)) * transf_theta[2] * E_nonzero
  d2 <- E_nonzero * (1 - E_nonzero * exp(-lambda)) * p_nonzero
  jacobian <- d1 + d2
  return(list(fit = fit, jacobian = jacobian))
}

multinom_jacobian <- function(lp_i, family_object){
  
  aug_lpi <- cbind(0, lp_i)
  lp_max <- do.call(pmax, data.frame(aug_lpi))
  lp_transf <- sweep(aug_lpi, MARGIN = 1, STATS = lp_max, FUN = '-')
  lp_transf <- exp(lp_transf)
  row_transf <- rowSums(lp_transf)
  lp_out <- lp_transf/row_transf
  lp_denom <- exp(lp_max) * row_transf  
  
  lp_jacob <- lapply(1:ncol(lp_out), FUN=function(k){
    if (k == 1){
      fmt_k <- -lp_out[,-1, drop = F] / lp_denom
    }else{
      fmt_k <- -lp_out[,-1, drop = F] * lp_out[,k]
      fmt_k[,k-1] <- fmt_k[, k-1, drop = F] + lp_out[,k]
    }
    return(fmt_k)
  })
  return(
    list(fit = lp_out, jacobian = lp_jacob)
  )
}
