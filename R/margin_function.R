

#' @import prediction
prediction.gam <- function (model, data = find_data(model, parent.frame()), at = NULL, 
          type = c("response", "link"), vcov = stats::vcov(model), 
          calculate_se = TRUE, ...) 
{
  
  mdf <- utils::getFromNamespace("make_data_frame", "prediction")
  model_coef <- coef(model)

  type <- match.arg(type)
  data <- data
  
  if (missing(data) || is.null(data)) {
    if (isTRUE(calculate_se)) {
      pred <- predict(model, type = type, se.fit = TRUE, 
                      ...)
      pred <- mdf(fitted = pred[["fit"]], 
                              se.fitted = pred[["se.fit"]])
    }
    else {
      pred <- predict(model, type = type, se.fit = FALSE, 
                      ...)
      pred <- mdf(fitted = pred, se.fitted = rep(NA_real_, length(pred)))
    }
  }
  else {
    model[["model"]] <- NULL
    out <- build_datalist(data, at = at, as.data.frame = TRUE)
    at_specification <- attr(out, "at_specification")
    if (isTRUE(calculate_se)) {
      tmp <- predict(model, newdata = out, type = type, 
                     se.fit = TRUE, ...)
      pred <- mdf(out, fitted = tmp[["fit"]], 
                              se.fitted = tmp[["se.fit"]])
    }
    else {
      tmp <- predict(model, newdata = out, type = type, 
                     se.fit = FALSE, ...)
      pred <- mdf(out, fitted = tmp, se.fitted = rep(NA_real_, 
                                                                 nrow(out)))
    }
  }
  
  if (isTRUE(calculate_se)) {
    model_terms <- delete.response(terms(model))
    if (is.null(at)) {
      
      model_mat <- predict(object = model, na.action = na.pass,
                           newdata = data, type = 'lpmatrix')
      if (type == "link") {
        means_for_prediction <- colMeans(model_mat)
      }
      else if (type == "response") {
        predictions_link <- model_mat %*% model_coef
        v <- as.vector(model$family$mu.eta(predictions_link))
        means_for_prediction <- colMeans(v * model_mat)
        
      }
      J <- matrix(means_for_prediction, nrow = 1L)
    }
    else {
      datalist <- build_datalist(data, at = at, as.data.frame = FALSE)
      jacobian_list <- lapply(datalist, function(one) {
        
        model_mat <- predict(object = model, newdata = one, 
            na.action = na.pass, type = 'lpmatrix')
        
        if (type == "link") {
          means_for_prediction <- colMeans(model_mat)
        }
        else if (type == "response") {
          
          predictions_link <- model_mat %*% model_coef
          v <- as.vector(model$family$mu.eta(predictions_link))
          means_for_prediction <- colMeans(v * model_mat)
        }
        means_for_prediction
      })
      J <- do.call("rbind", jacobian_list)
    }

    vc <- diag(J %*% vcov %*% t(J))
  }
  else {
    J <- NULL
    if (length(at)) {
      vc <- rep(NA_real_, nrow(at_specification))
    }
    else {
      vc <- NA_real_
    }
  }
  structure(pred, class = c("prediction", "data.frame"), 
            at = if (is.null(at)) 
              at
            else at_specification, type = type, call = if ("call" %in% 
                                                           names(model)) 
              model[["call"]]
            else NULL, model_class = class(model), row.names = seq_len(nrow(pred)), 
            vcov = vc, jacobian = J, weighted = FALSE)
}

#' @import margins
margins.gam <- function(
  model, data = find_data(model, parent.frame()), variables = NULL, 
  at = NULL, type = c("response", "link"), vcov = stats::vcov(model), 
  vce = c("delta", "simulation", "bootstrap", "none"), iterations = 50L, 
  unit_ses = FALSE, eps = 1e-07, 
  ...
){
  
  type <- match.arg(type)
  vce <- match.arg(vce)
  data_list <- build_datalist(data, at = at)
  if (is.null(names(data_list))) {
    names(data_list) <- NA_character_
  }
  at_specification <- attr(data_list, "at_specification")
  varslist <- find_terms_in_model.gam(model, variables = variables)
  model[["model"]] <- NULL
  attr(model[["terms"]], ".Environment") <- NULL
  out <- list()
  bmarg <- utils::getFromNamespace("build_margins", "margins")
  
  for (i in seq_along(data_list)) {
    out[[i]] <- bmarg(model = model, data = data_list[[i]], 
                                     variables = variables, type = type, vcov = vcov, 
                                     vce = vce, iterations = iterations, unit_ses = unit_ses, 
                                     eps = eps, varslist = varslist, ...)
    out[[i]][["_at_number"]] <- i
  }
  if (vce == "delta") {
    jac <- do.call("rbind", lapply(out, attr, "jacobian"))
    rownames(jac) <- paste0(rownames(jac), ".", rep(seq_len(length(out)), 
                                                    each = length(unique(rownames(jac)))))
    vc <- jac %*% vcov %*% t(jac)
  }
  else {
    jac <- NULL
    vc <- NULL
  }
  structure(do.call("rbind", out), class = c("margins", 
                                             "data.frame"), at = if (is.null(at)) 
                                               NULL
            else at_specification, type = type, call = if ("call" %in% 
                                                           names(model)) 
              model[["call"]]
            else NULL, model_class = class(model), vce = vce, vcov = vc, 
            jacobian = jac, weighted = FALSE, iterations = if (vce == 
                                                               "bootstrap") 
              iterations
            else NULL)
}

find_terms_in_model.gam <- function(model, variables = NULL){
  if (!is.null(attributes(terms(model))[["dataClasses"]])) {
    classes <- attributes(terms(model))[["dataClasses"]][-1L]
  }
  else if ("model" %in% names(model) && !is.null(attributes(terms(model$model))[["dataClasses"]])) {
    classes <- attributes(terms(model$model))[["dataClasses"]][-1L]
  }
  else {
    att <- attributes(terms(find_data(model)))
    if ("dataClasses" %in% names(att)) {
      classes <- att[["dataClasses"]][-1L]
      rm(att)
    }
    else {
      stop("No variable classes found in model.")
    }
  }
  classes <- classes[!names(classes) %in% "(weights)"]
  classes[classes == "character"] <- "factor"
  
  cterms <- utils::getFromNamespace('clean_terms', 'margins')
  names(classes) <- cterms(names(classes))
  vars <- list(nnames = unique(names(classes)[!classes %in% 
                                                c("factor", "ordered", "logical")]), 
               lnames = unique(names(classes)[classes == "logical"]), 
               fnames = unique(names(classes)[classes %in% c("factor", 
                                                             "ordered")]))
  if (!is.null(variables)) {
    tmp <- c(vars$nnames, vars$lnames, vars$fnames)
    if (any(!variables %in% tmp)) {
      stop("Some values in 'variables' are not in the model terms.")
    }
    vars$nnames <- vars$nnames[vars$nnames %in% variables]
    vars$lnames <- vars$lnames[vars$lnames %in% variables]
    vars$fnames <- vars$fnames[vars$fnames %in% variables]
  }
  if (is.null(unlist(vars))) {
    stop("No variables found in model.")
  }
  return(vars)
}


