
#' @export
print.gKRLS <- function(object){
  print(object$fmt_varcorr)
}

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
#' Prediction after gKRLS
#' Pass "newdata"
#' @export
predict.gKRLS <- function(object, newdata, newkernel_X,
                          type = c('link', 'response'),
                          return_full = FALSE){

    type <- match.arg(type)
  
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
    
    if (is.null(newdata)){
      newdata <- data.frame(matrix(nrow = nrow(newkernel_X), ncol = 0))
    }
    newdata_FE <- model.matrix(delete.response(terms(lme4::nobars(formula(object)))), data = newdata)
    
    orig_X_names <- names(fixef(object))
    if (!identical(colnames(newdata_FE), orig_X_names)) {
      print(all.equal(colnames(newdata_FE), orig_X_names))
      stop("Misaligned Fixed Effects")
    }
    
    fit_fe <- fixef(object)
    fit_re <- ranef(object)

    fit_vcov_ridge <- object$vcov_ridge
    yfitted <- as.vector(newdata_FE %*% fit_fe + newdataKS %*% fit_re)
    

    # Decompose the variance term into L * L^T to get the 
    # standard error for each predicted value.
    
    decomp_re <- eigen(fit_vcov_ridge)
    decomp_re <- with(decomp_re, vectors %*% 
            Diagonal(x = sqrt(ifelse(values < 0, 0, values))))
    yfitted_se <- sqrt(rowSums(as.matrix(cbind(newdata_FE, newdataKS) %*% decomp_re)^2 ))

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
    out <- list(fitted = yfitted, se = yfitted_se)
    if (return_full){
      out$full_vcov <- (newdataKS %*% fit_re_var %*% t(newdataKS)) * 
        object$internal$sd_y^2
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