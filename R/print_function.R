
print.gKRLS <- function(object){
  print(object$fmt_varcorr)
}

#' Prediction after gKRLS
#' Pass "newdata"
#' @export
predict.gKRLS <- function(object, newdata, newkernel_X,
                          type = c('link', 'response'),
                          return_full = FALSE){

    type <- match.arg(type)
  
    # Standardize the incoming new data.
    newkernel_X <- sweep(newkernel_X, 2, object$internal$std_train$mean, FUN = "-")
    newkernel_X <- newkernel_X %*% 
      Diagonal(x = 1/sqrt(object$internal$std_train$var)) %*% 
      object$internal$std_train$whiten
    newkernel_X <- as.matrix(newkernel_X)  
  
    # Get the "test" data after using same sketch matrix
    newdataKS <- create_sketched_kernel(X_test = newkernel_X, 
        X_train = object$internal$kernel_X_train, 
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