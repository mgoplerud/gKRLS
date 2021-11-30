
#' @export
print.gKRLS <- function(object){
  print(object$fmt_varcorr)
}

#' @export
plot.gKRLS_ME <- function(object){
  
  list_pointwise <- as.list(data.frame(object$ME_pointwise, check.names = F))
  list_pointwise_var <- as.list(data.frame(object$ME_pointwise_var, check.names = F))
  
  pointwise_flat <- do.call('rbind', mapply(names(list_pointwise), list_pointwise, list_pointwise_var, SIMPLIFY = FALSE, 
         FUN=function(name, est, var){
           data.frame(variable = name, id = 1:length(est), est = est, se = sqrt(var), check.names = F, stringsAsFactors = F)
         }))
  rownames(pointwise_flat) <- NULL
  
  out <- data.frame(est = object$AME_pointwise, se = sqrt(object$AME_pointwise_var))
  
  thresh <- abs(qt(0.025, df = nrow(object$ME_pointwise)))
  
  out$t.stat <- out$est/out$se
  out$p.value <- 2 * pt(-abs(out$t.stat), df = nrow(object$ME_pointwise))
  out$variable <- rownames(out)
  
  out$variable <- factor(out$variable, levels = out$variable[order(out$est)])
  pointwise_flat$variable <- factor(pointwise_flat$variable, levels = levels(out$variable))

  g <- ggplot(pointwise_flat) +
    geom_boxplot(aes(x=est,y=variable, col = 'black'), linetype = 'dashed', alpha = 0.1) +
    geom_point(aes(x=est,y=variable), col = 'red', data = out) +
    geom_errorbarh(aes(xmin=est - thresh * se,xmax=est + thresh * se,y=variable), col = 'red',
                    data = out) +
    theme_bw() +
    scale_color_manual(values = rgb(0,0,0, alpha = 0.3), guide = FALSE) +
    geom_vline(aes(xintercept=0)) + xlab('Marginal Effect') +
    ylab('Varible')
  plot(g)
  invisible(g)
}

#' @export
print.gKRLS_ME <- function(object){
  
  summary_pointwise <- apply(object$ME_pointwise, MARGIN = 2, FUN=function(i){quantile(i, c(0.25, 0.5, 0.75))})
  cat(paste0('Distribution of Pointwise Marginal Effects: N = ', nrow(object$ME_pointwise), '\n'))
  print(summary_pointwise)
  
  out <- data.frame(est = object$AME_pointwise, se = sqrt(object$AME_pointwise_var))
  
  thresh <- abs(qt(0.025, df = nrow(object$ME_pointwise)))
  
  out$t.stat <- out$est/out$se
  out$p.value <- 2 * pt(-abs(out$t.stat), df = nrow(object$ME_pointwise))
  cat('\nSummary of Average Marginal Effects\n')
  print(out)
  invisible()  
}

prepare_predict_data <- function(object, newdata, newkernel_X, allow_missing_levels){

  # Standardize the incoming new data.
  std_newkernel_X <- sweep(newkernel_X, 2, object$internal$std_train$mean, FUN = "-")
  
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
  newdataKS <- create_sketched_kernel(X_test = std_newkernel_X, 
                                      X_train = std_kernel_X, 
                                      tS = t(object$internal$sketch), 
                                      bandwidth = object$internal$bandwidth)
  rownames(newdataKS) <- as.character(1:nrow(newdataKS))
  
  if (is.null(newdata)){
    newdata <- data.frame(matrix(nrow = nrow(newkernel_X), ncol = 0))
  }else{
    rownames(newdata) <- as.character(1:nrow(newdata))
  }
  
  fmla <- formula(object)
  fmla <- update.formula(fmla, '. ~ . + (1 | kernel_RE)')
  # Adapted from "predict.lm"
  tt <- object$internal$fe_options$terms
  Terms <- delete.response(tt)
  m <- model.frame(Terms, newdata, na.action = object$internal$fe_options$na.action, 
                   xlev = object$internal$fe_options$xlevels)
  if (!is.null(cl <- attr(Terms, "dataClasses"))){
    .checkMFClasses(cl, m)
  } 
  
  newdata_FE <- model.matrix(Terms, m, contrasts.arg = contrasts)

  if (!identical(colnames(newdata_FE), names(object$fe$mean))) {
    print(all.equal(colnames(newdata_FE), names(object$fe$mean)))
    stop("Misaligned Fixed Effects")
  }
  
  # Adapted from vglmer (github.com/mgoplerud/vglmer)
  # Extract the Z (Random Effect) design matrix.
  newdata$kernel_RE <- 1
  mk_Z <- model.frame(delete.response(terms(subbars(fmla))), data = newdata)
  if (nrow(newdataKS) != nrow(mk_Z)){stop('Kernel data size misaligned with RE/FE.')}
  if (nrow(mk_Z) != nrow(newdata_FE)){stop('Data size for RE and FE are misaligned.')}
  rownames_Z <- rownames(mk_Z)
  mk_Z <- lme4::mkReTrms(lme4::findbars(fmla), mk_Z, reorder.terms = FALSE, reorder.vars = FALSE)
  mk_Z$Zt <- NULL
  
  Z <- mapply(names(mk_Z$Ztlist), mk_Z$Ztlist, SIMPLIFY = FALSE, FUN=function(name_zj, zj){
    if (name_zj == "1 | kernel_RE"){return(NA)}
    levels_training <- object$levels_of_RE[[name_zj]]
    levels_test <- rownames(zj)
    
    not_in_train <- setdiff(levels_test, levels_training)
    not_in_test <- setdiff(levels_training, levels_test)
    
    if (length(not_in_train) > 0) {
      if (!allow_missing_levels) {
        stop("New levels not allowed unless allow_missing_levels = TRUE")
      }
    }
    
    in_both <- intersect(levels_training, levels_test)
    recons_Z <- drop0(sparseMatrix(i = 1, j = 1, x = 0, dims = c(ncol(zj), length(levels_training))))
    colnames(recons_Z) <- levels_training
    rownames(recons_Z) <- colnames(zj)
    
    recons_Z[, match(in_both, levels_training)] <- t(zj[match(in_both, rownames(zj)),])
    return(recons_Z)      
  })
  
  
  Z[["1 | kernel_RE"]] <- newdataKS
  n_pos <- sapply(Z, ncol)
  Z <- do.call('cbind', Z)
  
  gc()
  
  total_obs <- rownames(newdata)
  obs_in_both <- intersect(rownames(newdata_FE), rownames(Z))

  newdata_FE <- newdata_FE[match(obs_in_both, rownames(newdata_FE)), , drop = F]
  Z <- Z[match(obs_in_both, rownames(Z)), , drop = F]
  
  
  return(list(
    Z = Z,
    n_pos = n_pos,
    newdata_FE = newdata_FE,
    total_obs = total_obs,
    obs_in_both = obs_in_both,
    newdataKS = newdataKS,
    std_newkernel_X = std_newkernel_X,
    std_kernel_X = std_kernel_X
  ))
}

#' Prediction after gKRLS
#' Pass "newdata"
#' @importFrom lme4 nobars findbars mkReTrms
#' @export
predict.gKRLS <- function(object, newdata, newkernel_X,
                          type = c('link', 'response'),
                          return_full = FALSE,
                          allow_missing_levels = FALSE){

    type <- match.arg(type)
    
    prepped_data <- prepare_predict_data(object = object, newdata = newdata, 
                         newkernel_X = newkernel_X,
                         allow_missing_levels = allow_missing_levels)

    Z <- prepped_data$Z
    newdata_FE <- prepped_data$newdata_FE
    total_obs <- prepped_data$total_obs
    obs_in_both <- prepped_data$obs_in_both
    newdataKS <- prepped_data$newdataKS
    
    rm(prepped_data); gc()
    
    fit_fe <- fixef(object)
    fit_re <- ranef(object)
    fit_vcov_ridge <- object$vcov_ridge
    
    
    yfitted <- as.vector(newdata_FE %*% fit_fe + Z %*% fit_re)
    

    # Decompose the variance term into L * L^T to get the 
    # standard error for each predicted value.
    
    decomp_re <- eigen(fit_vcov_ridge)
    decomp_re <- with(decomp_re, vectors %*% 
            Diagonal(x = sqrt(ifelse(values < 0, 0, values))))
    yfitted_se <- sqrt(rowSums(as.matrix(cbind(newdata_FE, Z) %*% decomp_re)^2 ))

    # yfitted_se_man <- apply(newdataKS, MARGIN = 1, 
    #   FUN=function(i){sqrt(as.numeric(t(i) %*% fit_re_var %*% i))})
    
    if (type == 'response' & family(object)$family != 'gaussian'){
      yfitted_se <- yfitted_se * abs(family(object)$mu.eta(yfitted))
      yfitted <- family(object)$linkinv(yfitted)
    }
    
    if (family(object)$family == 'gaussian'){
      yfitted <- (yfitted * object$internal$sd_y) + object$internal$mean_y
      yfitted_se <- yfitted_se * object$internal$sd_y
    }
    
    # Reinject missing data if needed
    
    yfitted <- yfitted[match(total_obs, obs_in_both)]
    yfitted_se <- yfitted_se[match(total_obs, obs_in_both)]
    
    out <- list(fitted = yfitted, se = yfitted_se)
    
    if (return_full){
      out$full_vcov <- (cbind(newdata_FE, Z) %*% fit_vcov_ridge %*% t(cbind(newdata_FE, Z))) * 
        object$internal$sd_y^2
      if (length(total_obs) != length(obs_in_both)){stop('"return_full" cannot be used with missing data.')}
      
    }
    
    return(out)
}

fitted.gKRLS <- function(object){
  object$fitted
}

family.gKRLS <- function(object){
  object$internal$family
}

coef.gKRLS <- function(object){object$fe$mean}
vcov.gKRLS <- function(object){as.matrix(object$fe$vcov)}