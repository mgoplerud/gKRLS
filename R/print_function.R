
print.gKRLS <- function(object){
  print(object$fmt_varcorr)
}

#' Prediction after gKRLS
#' Pass "newdata"
#' @export
predict.gKRLS <- function(object, newdata, 
                          type = c('link', 'response'),
                          return_full = FALSE){

    type <- match.arg(type)
  
    # Standardize the incoming new data.
    newdata <- sweep(newdata, 2, object$internal$std_train$mean, FUN = "-")
    newdata <- newdata %*% 
      Diagonal(x = 1/sqrt(object$internal$std_train$var)) %*% 
      object$internal$std_train$whiten
    newdata <- as.matrix(newdata)  
  
    # Get the "test" data after using same sketch matrix
    newdataKS <- create_sketched_kernel(X_test = newdata, 
        X_train = object$internal$X_train, 
        tS = t(object$internal$sketch), 
        bandwidth = object$internal$bandwidth)

    fit_fe <- fixef(object)
    fit_re <- ranef(object)
    fit_re_var <- ranef(object, type = 'variance')
    if (length(fit_fe) != 0){stop('Not set up for FE yet...')}

    yfitted <- as.vector(newdataKS %*% fit_re)
    # Decompose the variance term into L * L^T to get the 
    # standard error for each predicted value.
    
    decomp_re <- eigen(fit_re_var)
    decomp_re <- with(decomp_re, vectors %*% 
            Diagonal(x = sqrt(ifelse(values < 0, 0, values))))
    yfitted_se <- sqrt(rowSums(as.matrix(newdataKS %*% decomp_re)^2 ))

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
