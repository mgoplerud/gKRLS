
#' @export
discrete_mfx <- function(model, data_list, vcov,
                         weights, individual = FALSE){
  
  model_coef <- coef(model)
  
  if (missing(vcov)){
    vcov <- vcov(model)
  }else if (identical(vcov, 'none')){
    vcov <- NULL
  }else{
    if (nrow(vcov) != length(coef(model))){
      stop('If vcov is provided manually, it must be the same size as the coefficient vector.')
    }
  }
  
  if (missing(weights)){
    if (length(data_list) != 2){
      stop('weights may only be null if data_list contains exactly two elements.')
    }else{
      message('Reporting data_list[[2]] - data_list[[1]]')
      weights <- c(-1, 1)
    }
  }else{
    if (length(weights) != length(data_list)){
      stop('"weights" must be the same length as "data_list".')
    }
  }
  
  raw_predictions <- lapply(data_list, FUN=function(data_i){
    # Get the design
    matrix_i <- predict(model, newdata = data_i, na.action = na.pass, type = 'lpmatrix')
    lp_i <- as.vector(matrix_i %*% model_coef)
    e_i <- model$family$linkinv(lp_i)
    ex <- mean(e_i, na.rm=T)
    if (individual){
      jacob_i <- Diagonal(x = model$family$mu.eta(lp_i)) %*% matrix_i
      jacob <- colMeans(jacob_i, na.rm=T)
      nrow_valid <- sum(!is.na(lp_i))
    }else{
      se_i <- NULL
      e_i <- NULL
      jacob_i <- NULL
      nrow_valid <- NULL
      jacob <- colMeans(Diagonal(x = model$family$mu.eta(lp_i)) %*% matrix_i, na.rm=T)
    }
    return(list(
      expectation = ex,
      expectation_i = e_i,
      jacobian_i = jacob_i,
      jacobian = jacob,
      nrow_valid = nrow_valid
    ))
  })
  
  jacobian_net <- sapply(raw_predictions, FUN=function(i){i$jacobian}) %*% weights
  out_se <- sqrt(as.numeric(t(jacobian_net) %*% vcov %*% jacobian_net))
  out_est <- as.numeric(sapply(raw_predictions, FUN=function(i){i$expectation}) %*% weights)
  
  out_aggregate <- data.frame(est = out_est, se = out_se)  
  if (individual){
    checksum_ind <- length(unique(sapply(raw_predictions, FUN=function(i){i$nrow_valid})))
    if (checksum_ind != 1){
      stop('individual=TRUE requires same pattern of missing data across all elements of data_list.')
    }
    jacob_net_i <- Reduce('+', mapply(sapply(raw_predictions, FUN=function(i){i$jacobian_i}), weights, FUN=function(i,j){i * j}))
    out_se_i <- apply(jacob_net_i, MARGIN = 1, FUN=function(i){sqrt(as.numeric(t(i) %*% vcov %*% i))})
    out_est <- as.vector(sapply(raw_predictions, FUN=function(i){i$expectation_i}) %*% weights)
    out_individual <- data.frame(est = out_est, se = out_se)
  }else{
    out_individual <- NULL
  }
  return(list(
    aggregate = out_aggregate,
    individual = out_individual
  ))
}


