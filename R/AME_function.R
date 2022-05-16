#' Marginal Effects by Numerical Derivatives
#' 
#' This function calculates the marginal effects using numeric approximations of 
#' the partial derivatives. For a description of this method, please see Thomas Leeper's 
#' article "Interpreting Regression Results using Average Marginal Effects with R’s margins."
#' 
#' @name calculate_effects
#' @param model A fitted gKRLS model.
#' @param object A fitted gKRLS model.
#' @param data A new data frame that used to calculate the marginal effect, or 
#' set to ``NULL'', which the data used to estimate the model will be used. The 
#' default is ``NULL.''
#' @param variables Specify the variable names that need to calculate marginal 
#' effect. The default is ``NULL'', which means calculate marginal effect for all variables.
#' @param vcov Specify the covariance matrix.It accepts a user-defined covariance 
#' matrix or clustered covariance matrices using functions from sandwich package. 
#' @param raw Placeholder.
#' @param individual Calculate individual effects (i.e. an effect for each observation in the provided data).
#' @param conditional  This is an analogue of Stata's ``at()'' option and ``at'' argument 
#' in ``margins'' package. Specify the values at which to calculate the marginal 
#' effect in a named data from. See an example below.
#' @param epsilon A numerical value to define the step when calculating numerical 
#' derivatives. See Leeper's articel for details. 
#' @param verbose A logical value indicates whether to report the current stage 
#' when calculating the marginal effects.
#' @param continuous_type A character string indicating the type of marginal effects 
#' to estimate. Options are ``IQR'': variable values change from 25% to 75%, ``minmax'': 
#' variable values changes from minimum to maximum, ``derivative'': variable values 
#' change by epsilon defined in the epsilon argument, ``onesd'': variable values 
#' change by one standard deviation.
#' 
#' @return  \code{calculate_effects} return a list which contain a data frame for marginal effects. 
#' The \code{variable} represents the variable names used to calculate marginal effects. \code{type} 
#' indicates the method for calculating marginal effects. It is one of the IQR, minmax, derivative,
#' and onesd. The \code{est} represents the marginal effects. The \code{se}, \code{t}, \code{p.value} are 
#' standard error, t-value, and p-value for marginal effect. This function also return several information 
#' such as Jacobian matrix (\code{jacobian}), how many variables specified to calculate marginal effect (\code{counter}),
#' and effective sample size (\code{N_eff}) and sample size (N). 
#' 
#' @examples
#' n <- 50
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' state <- sample(letters[1:5], n, replace = TRUE)
#' y = 0.3*x1 + 0.4*x2 +0.5*x3 + rnorm(n)
#' data <- data.frame(y, x1, x2, x3, state)
#' 
#' # Make character variables into factors for mgcv
#' data$state <- factor(data$state)
#' 
#' # A gKRLS model 
#' gkrls_est <- mgcv::gam(y ~ state + s(x1,x2,x3, bs="gKRLS"), data = data)
#' 
#' # calculate marginal effect using derivative
#' calculate_effects(gkrls_est, variables = "x1", continuous_type = 'derivative') 
#' 
#' # calculate marginal effect by specifying conditional variables
#' calculate_effects(gkrls_est, variables = "x1", 
#' conditional = data.frame(x2 = c(0.6, 0.8), x3=0.3))
#' 
#' # calculate marginal effect by specifying a factor conditional variable
#' calculate_effects(gkrls_est, variables = "x1", 
#' conditional = data.frame(state = c("a", "b", "c")), continuous_type = 'derivative')
#' 
#' @importFrom stats model.frame sd
#' @export
calculate_effects <- function(model, data = NULL, 
  variables = NULL, vcov = NULL, raw = FALSE, individual = FALSE,
  conditional = NULL, epsilon = 1e-7, verbose = FALSE,
  continuous_type = c('IQR', 'minmax', 'derivative', 'onesd')){
  
  if (!is.list(continuous_type)){
    continuous_type <- match.arg(continuous_type)
  }

  N_eff <- length(model$y) - sum(model$edf)
  N <- length(model$y)

  if (is.null(data)){
    raw_data <- model.frame(model)
  }else{
    raw_data <- data
    rm(data)
  }
  
  if (is.null(vcov)){
    vcov <- stats::vcov(model)
  }
  
  all_variables <- gam_terms(model = model, variables = NULL)
  variable_list <- unlist(all_variables)
  variable_type <- rep(names(all_variables), lengths(all_variables))
  names(variable_type) <- unlist(all_variables)
  
  # Address the conditional variables
  if (!is.null(conditional)){
    
    names_cond <- names(conditional)
    # Check the variables in conditional are valid
    if ( (!all(names_cond %in% variable_list)) & (length(setdiff(names_cond, variable_list)) != 0) ){
      stop('all variables in "conditional" must be in the original model.')
    }
    # Remove those variables from the ones to calculate the AME over
    # all_variables <- find_terms_in_model.gam(model = model, 
    #   variables = setdiff(names(variable_type), names_cond))
    # if (length(unlist(all_variables)) == 0){
    #   stop('After removing variables in conditional, none are left.')
    # }
    # variable_list <- unlist(all_variables)
    # variable_type <- rep(names(all_variables), lengths(all_variables))
    # names(variable_type) <- unlist(all_variables)
    ncond <- nrow(conditional)
  }else{
    ncond <- 1
  }

  any_binary <- apply(raw_data[,names(variable_type)], 
                      MARGIN = 2, 
                      FUN=function(i){all(i %in% c(0,1))})
  variable_type[which(any_binary)] <- 'bnames'
  
  if (is.null(variables)){
    variable_list <- as.list(variable_list)
  }else if (is.list(variables)){
    check_valid <- all(sapply(variables, FUN=function(i){
      all(i %in% variable_list)
    }))
    if (!check_valid){
      stop('All elements of list of "variables" must be in model.')
    }
    check_valid <- all(sapply(variables, anyDuplicated) == 0)
    if (!check_valid){
      stop('All elements of a list of "variables" must be distinct.')
    }
    variable_list <- variables
  }else if (is.vector(variables)){
    if (!all(variables %in% variable_list)){
      stop('All elements of "variables" must be in model.')
    }
    if (anyDuplicated(variables)){stop('All elements of "variables" must be unique.')}
    variable_list <- as.list(variables)
  }else{stop('"variables" must be NULL, list or a vector.')}
  
  out_mfx <- data.frame()
  if (individual){
    out_mfx_individual <- data.frame()
  }
  out_jacobian <- matrix(nrow = length(coef(model)), ncol = 0)
  out_counter <- c()
  
  if (is.list(continuous_type)){
    
    v <- unlist(variable_list)
    v <- variable_type[match(v, names(variable_type))]
    v <- v[v == 'nnames']
    if (!all(names(v) %in% names(continuous_type))){
      stop('if continuous_type is a list, then it must have all numerical variables in it.')
    }
    if (!all(lengths(continuous_type) == 2)){
      stop('if continous_type is a list, it must have vectors of two elements.')
    }
  }
  
  continuous_f <- function(x, m, s, type){
    if (type == 'list'){
      s
    }else if (type == 'derivative'){
      x + s
    }else if (type == 'onesd'){
      m + s
    }else if (type %in% c('minmax', 'IQR')){
      s
    }else{NA}
  }
  for (cond_i in seq_len(ncond)){
    if (verbose & cond_i %% ceiling(ncond/ncond) == 0){
      message('|')
    }
    data <- raw_data
    
    if (!is.null(conditional)){
      cond_data_i <- conditional[cond_i,,drop=F]
      for (j in names(cond_data_i)){
        data[[j]] <- cond_data_i[[j]]
      }
    }
  
    data_out <- list()
    for (v in variable_list){
      
      type_v <- variable_type[v]
      
      packaged_data <- list()
      # Loop over each variable and create the "d0" and "d1" frames, respectively
      for (v_id in seq_len(length(v))){
        v_i <- v[v_id]
        type_i <- type_v[v_id]
        multiplier <- 1
        if (type_i == 'nnames'){
          
          if (is.list(continuous_type)){
            r_i <- continuous_type[[v_i]]
            step2 <- c('0' = r_i[1], '1' = r_i[2])
            step <- NULL
            ctype <- 'list' 
          }else if (continuous_type == 'derivative'){
            # Closely adapted from "margins" by Thomas Leeper
            step <- max(abs(data[[v_i]]), 1, na.rm=T) * epsilon
            multiplier <- 1/(2 * step)
            step2 <- c('0' = -step, '1' = step)
            step <- NULL
            ctype <- 'derivative'
          }else if (continuous_type == 'onesd'){
            step <- mean(data[[v_i]], na.rm=T)
            sd_i <- sd(data[[v_i]], na.rm=T)
            step2 <- c('0' = -sd_i, '1' = sd_i)
            ctype <- 'onesd'
          }else if (continuous_type == 'minmax'){
            r_i <- range(data[[v_i]], na.rm=T)
            names(r_i) <- NULL
            step2 <- c('0' = r_i[1], '1' = r_i[2])
            step <- NULL
            ctype <- 'minmax'
          }else if (continuous_type == 'IQR'){
            q_i <- quantile(data[[v_i]], c(0.25, 0.75), na.rm=T)
            names(q_i) <- NULL
            step2 <- c('0' = q_i[1], '1' = q_i[2])
            step <- NULL
            ctype <- 'IQR'
          }else{stop('invalid continuous_type')}
          if (v_id == 1){
            packaged_data <- list(list(
              data = list('d1' = data, 'd0' = data),
              weights = c(1, -1) * multiplier,
              raw = raw,
              type = continuous_type,
              name = v_i
            ))
            packaged_data[[1]]$data$d0[[v_i]] <- continuous_f(packaged_data[[1]]$data$d0[[v_i]], step, step2["0"], ctype)
            packaged_data[[1]]$data$d1[[v_i]] <- continuous_f(packaged_data[[1]]$data$d1[[v_i]], step, step2["1"], ctype)
            names(packaged_data) <- v_i
          }else{
            old_names <- names(packaged_data)
            packaged_data <- lapply(packaged_data, FUN=function(i){
              i$data$d0[[v_i]] <- continuous_f(i$data$d0[[v_i]], step, step2["0"], ctype)
              i$data$d1[[v_i]] <- continuous_f(i$data$d1[[v_i]], step, step2["1"], ctype)
              if (!any(i$type %in% 'derivative')){
                i$weights <- i$weights * multiplier
              }
              i$raw <- i$raw * raw
              i$type <- c(i$type, continuous_type)
              i$name <- c(i$name, v_i)
              return(i)
            })
            names(packaged_data) <- paste(names(packaged_data), ':', v_i, sep = '')
          }
          
        }else if (type_i == 'bnames'){
          if (v_id == 1){
            packaged_data <- list(list(
              data = list('d1' = data, 'd0' = data),
              weights = c(1, -1),
              raw = raw,
              type = 'binary',
              name = v_i
            ))
            packaged_data[[1]]$data$d0[[v_i]] <- 0
            packaged_data[[1]]$data$d1[[v_i]] <- 1
            names(packaged_data) <- v_i
          }else{
            old_names <- names(packaged_data)
            packaged_data <- lapply(packaged_data, FUN=function(i){
              i$data$d0[[v_i]] <- 0
              i$data$d1[[v_i]] <- 1
              i$weights <- i$weights
              i$raw <- i$raw * raw
              i$type <- c(i$type, 'binary')
              i$name <- c(i$name, v_i)
              return(i)
            })
            names(packaged_data) <- paste(names(packaged_data), ':', v_i, sep = '')
          }
        }else if (type_i == 'lnames'){
          if (v_id == 1){
            packaged_data <- list(list(
              data = list('d1' = data, 'd0' = data),
              weights = c(1, -1),
              raw = raw,
              type = 'logical',
              name = v_i
            ))
            packaged_data[[1]]$d0[[v_i]] <- FALSE
            packaged_data[[1]]$d1[[v_i]] <- TRUE
            
            names(packaged_data) <- v_i
          }else{
            stop('Invalid logical data.')
          }
        }else if (type_i == 'fnames'){
          levs <- levels(as.factor(data[[v_i]]))
          base <- levs[1L]
          levs <- levs[-1L]
          if (v_id == 1){
            packaged_data <- list()
            for (i in seq_along(levs)) {
              tmp_name <- paste0("factor(", v_i, ")", levs[i])
              packaged_data[[tmp_name]] <- list(
                data = list('d1' = data, 'd0' = data),
                weights = c(1, -1),
                raw = raw,
                type = 'factor',
                name = tmp_name
              )
              packaged_data[[tmp_name]]$data$d0[[v_i]] <- base
              packaged_data[[tmp_name]]$data$d1[[v_i]] <- levs[i]
              
            }
          }else{
            
            old_names <- names(packaged_data)
            added_data <- list()
            for (l in seq_along(levs)){
              temp_name <- paste0("factor(", v_i, ")", levs[l])
              add_l <- mapply(packaged_data, names(packaged_data), SIMPLIFY = FALSE, FUN=function(i, m){
                i$data$d0[[v_i]] <- base
                i$data$d1[[v_i]] <- levs[l]
                i$weights <- i$weights
                i$raw <- i$raw * raw
                i$type <- c(i$type, 'factor')
                i$name <- c(i$name, temp_name)
                return(i)
              })
              names(add_l) <- paste0(names(packaged_data), ':', temp_name)
              added_data <- c(added_data, add_l)
            }
            packaged_data <- added_data
          }
        }else{
          stop('Unknown type')
        }
      }
      gc()
      
      fit_mfx <- lapply(packaged_data, FUN=function(data_i){
        fit_i <- weighted_mfx(model = model, 
            data_list = data_i$data, 
            weights = list('AME' = data_i$weights),
            vcov = vcov, 
            individual = individual, raw = data_i$raw)
        if (data_i$raw){
          fit_i$aggregate[["...id"]] <- c('effect', 'raw_1', 'raw_0')
        }
        fit_i$name <- paste(data_i$name, collapse = ':')
        fit_i$type <- paste(data_i$type, collapse = ':')
        return(fit_i)
      })
      data_out <- c(data_out, fit_mfx)
    }
    
    out_mfx_i <- do.call('rbind', lapply(data_out, FUN=function(i){
      i$aggregate$variable <- i$name
      i$aggregate$name <- NULL
      i$aggregate$type <- i$type
      return(i$aggregate)
    }))
    
    if (individual){
      out_mfx_i_ind <- do.call('rbind', lapply(data_out, FUN=function(i){
        i$individual$variable <- i$name
        i$individual$name <- NULL
        i$individual$type <- i$type
        return(i$individual)
      }))
      rownames(out_mfx_i_ind) <- NULL
    }
    
    if ('...id' %in% colnames(out_mfx_i)){
      id <- out_mfx_i[['...id']]
      if (individual){
        id_individual <- out_mfx_i_ind[['...id']]
      }
    }else{
      id <- NULL
      id_individual <- NULL
    }
    out_mfx_i <- out_mfx_i[,c('variable', 'type', 'est', 'se')]
    out_mfx_i[['...id']] <- id
    rownames(out_mfx_i) <- NULL
    
    if (individual){
      out_mfx_i_ind <- out_mfx_i_ind[,c('obs', 'variable', 'type', 'est', 'se')]
      out_mfx_i_ind[['...id']] <- id
    }
    
    if (!is.null(conditional)){
      for (j in names(cond_data_i)){
        out_mfx_i[[j]] <- cond_data_i[[j]]
        if (individual){
          out_mfx_i_ind[[j]] <- cond_data_i[[j]]
        }
      }
    }
    
    out_counter <- c(out_counter, rep(cond_i, nrow(out_mfx_i)))
    out_mfx <- rbind(out_mfx, out_mfx_i)
    out_jacobian_i <- do.call('cbind', lapply(data_out, FUN=function(i){i$jacobian}))
    out_jacobian <- cbind(out_jacobian, out_jacobian_i)
    
    if (individual){
      out_mfx_individual <- rbind(out_mfx_individual, out_mfx_i_ind)
    }
    
  }
  cat('\n')
  if (ncol(out_jacobian) != nrow(out_mfx)){
    stop('Unusual alignment error between jacobian and marginal effects.')
  }
  out_mfx$t <- out_mfx$est/out_mfx$se
  out_mfx$p.value <- 2 * pt(-abs(out_mfx$t), df = N_eff)

  out <- list(marginal_effects = out_mfx, 
              jacobian = out_jacobian,
              counter = out_counter)
  if (individual){
    out_mfx_individual$t <- out_mfx_individual$est/out_mfx_individual$se
    out_mfx_individual$p.value <- 2 * pt(-abs(out_mfx_individual$t), df = N_eff)
    out$individual <- out_mfx_individual
  }
  out$N_eff <- N_eff
  out$N <- N
  class(out) <- 'gKRLS_mfx'
  return(out)
}

kernel_interactions <- function(model, 
  variables, QOI = c('AMIE', 'ACE', 'AIE', 'AME'), ...){
  
  QOI <- match.arg(QOI, several.ok = TRUE)
  args <- list(...)
  if (isTRUE(args$raw)){
    stop('raw=T not permitted for interactions.')
  }
  if (!is.null(args$conditional)){
    if (any(names(args$conditional) %in% c('variable', 'est', 'se', 'QOI', '...counter'))){
      stop('conditional may not contain "variable", "est", "se", or "QOI" as column names.')
    }
  }
  if (any(lengths(variables) != 2)){
    stop('variables must be a list of two-length vectors.')
  }
  if (is.null(args$vcov)){
    vcov_mdl <- vcov(model)
  }else{
    vcov_mdl <- args$vcov
  }
  
  all_inter <- data.frame()
  for (v in variables){
    
    fmt_v <- c(list(v), as.list(v))
    
    fit_var <- do.call('calculate_effects', 
      c(list(model = model, variables = fmt_v), 
        args))
    id <- seq_len(nrow(fit_var$marginal_effects))
    split_var <- strsplit(fit_var$marginal_effects$variable, split = ':')
    split_var <- unique(split_var[which(lengths(split_var) == 2)])
    
    out_inter <- lapply(split_var, FUN=function(i){
      out <- data.frame()
      v_i <- paste(i, collapse = ':')
      id_two <- which(fit_var$marginal_effects$variable == v_i)
      id_one <- which(fit_var$marginal_effects$variable %in% i)
      id_two <- split(id_two, fit_var$counter[id_two])
      id_one <- split(id_one, fit_var$counter[id_one])  
      
      if ('AME' %in% QOI){
        out <- rbind(out,
           data.frame(
             fit_var$marginal_effects[unlist(id_one), c('variable', 'est', 'se')],
             QOI = 'AME',
             `...counter` = rep(seq_len(length(id_one)), lengths(id_one))
           )
        )
      }
      if ('ACE' %in% QOI){
        out <- rbind(out,
          data.frame(
            fit_var$marginal_effects[unlist(id_two), c('variable', 'est', 'se')],
            QOI = 'ACE',
            `...counter` = seq_len(length(unlist(id_two)))
          )
        )
      }
      if ('AMIE' %in% QOI){
        id_matrix <- do.call('rbind', mapply(id_two, id_one, seq_len(length(id_one)), SIMPLIFY = FALSE, FUN=function(x,y, z){
          rbind(c(x, z, 1), cbind(y, z, -1))
        }))
        id_matrix <- sparseMatrix(i = id_matrix[,1], j = id_matrix[,2], x = id_matrix[,3],
                                  dims = c(nrow(fit_var$marginal_effects), length(id_one)))
        point_estimate <- as.vector(fit_var$marginal_effects$est %*% id_matrix)
        se_estimate <- apply(fit_var$jacobian %*% id_matrix, 
                             MARGIN = 2, FUN = function(i){as.numeric(sqrt(t(i) %*% vcov_mdl %*% i))})
        out <- rbind(out,
            data.frame(variable = v_i, QOI = 'AMIE', est = point_estimate, se = se_estimate, `...counter` = seq_len(length(point_estimate)), stringsAsFactors = F)
        )
      }
      return(out)
    })
    out_inter <- do.call('rbind', out_inter)
    all_inter <- rbind(all_inter, out_inter)
  }
  rownames(all_inter) <- NULL
  all_inter <- unique(all_inter)
  all_inter$t <- all_inter$est/all_inter$se
  if ('...id' %in% names(all_inter)){
    id <- all_inter[['...id']]
  }else{
    id <- NULL
  }
  all_inter <- all_inter[,c('QOI', 'variable', 'est', 'se', 't', '...counter')]
  all_inter[["...id"]] <- id
  if (!is.null(args$conditional)){
    conditional <- args$conditional
    for (v in colnames(conditional)){
      all_inter[[v]] <- conditional[[v]][all_inter[['...counter']]]
    }
  }
  all_inter[["...counter"]] <- NULL
  out <- list(marginal_effects = all_inter, jacobian = NA)
  class(out) <- 'gKRLS_mfx'
  return(out)
}

# Adapted from Leeper's "margins"
gam_terms <- function(model, variables = NULL){
  
  classes <- attributes(terms(model))$dataClasses[-1]
  if (is.null(classes)){stop('terms not found from gam')}
  classes <- classes[!names(classes) %in% "(weights)"]
  classes[classes == "character"] <- "factor"
  # List of types, again following scheme in Leeper's "margins"
  vars <- list(nnames = unique(names(classes)[!classes %in% c("factor", "ordered", "logical")]), 
               lnames = unique(names(classes)[classes == "logical"]), 
               fnames = unique(names(classes)[classes %in% c("factor", "ordered")]))
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

#' @rdname calculate_effects
#' @param x Object fit with \code{calculate_effects} or \code{legacy_marginal_effect}
#' @method print gKRLS_mfx
#' @importFrom stats quantile pt
#' @export
print.gKRLS_mfx <- function(x, ...){
  
  if (!is.null(x$mfx_type)){

    summary_pointwise <- apply(x$ME_pointwise, MARGIN = 2, FUN=function(i){quantile(i, c(0.25, 0.5, 0.75))})
    cat(paste0('Distribution of Pointwise Marginal Effects: N = ', nrow(x$ME_pointwise), '\n'))
    print(summary_pointwise)
    
    out <- data.frame(est = x$AME_pointwise, se = sqrt(x$AME_pointwise_var))
    
    out$t.stat <- out$est/out$se
    out$p.value <- 2 * pt(-abs(out$t.stat), df = x$N_eff)
    cat('\nSummary of Average Marginal Effects\n')
    print(out)
  }else{
    cat('\nSummary of Average Marginal Effects\n\n')
    print(x$marginal_effects)
  }
}

#' @rdname calculate_effects
#' @method summary gKRLS_mfx
#' @param ... Additional arguments (unused).
#' @export
summary.gKRLS_mfx <- function(object, ...){object$marginal_effects}

#' @importFrom stats coef vcov na.pass predict
weighted_mfx <- function(model, data_list, vcov,
                         weights, raw = FALSE, individual = FALSE){
  
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
    stop('weights may not be null..')
  }else{
    if (any(lengths(weights) != length(data_list))){
      stop('"weights" must be list of vectors, each with the same length as "data_list".')
    }
  }
  weights <- do.call('cbind', weights)
  if (raw){
    add_cols <- diag(length(data_list))
    colnames(add_cols) <- paste0('raw_', 1:length(data_list))
    weights <- cbind(weights, add_cols)
  }
  
  # Get the predictions for each covariate profile provided
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
  # Get the jacobian/gradient for each weighted average
  jacobian_net <- sapply(raw_predictions, FUN=function(i){i$jacobian}) %*% weights
  
  out_se <- sqrt(apply(jacobian_net, MARGIN = 2, FUN=function(i){as.vector(t(i) %*% vcov %*% i)}))
  out_est <- as.numeric(sapply(raw_predictions, FUN=function(i){i$expectation}) %*% weights)
  
  out_aggregate <- data.frame(name = colnames(weights), est = out_est, se = out_se)  
  
  if (individual){
    checksum_ind <- length(unique(sapply(raw_predictions, FUN=function(i){i$nrow_valid})))
    if (checksum_ind != 1){
      stop('individual=TRUE requires same pattern of missing data across all elements of data_list.')
    }
    jacob_net_i <- Reduce('+', mapply(sapply(raw_predictions, FUN=function(i){i$jacobian_i}), weights, FUN=function(i,j){i * j}))
    # out_se_i <- apply(jacob_net_i, MARGIN = 1, FUN=function(i){sqrt(as.numeric(t(i) %*% vcov %*% i))})
    out_se_i <- sqrt(rowSums( (jacob_net_i %*% vcov) * jacob_net_i ))
    names(out_se_i) <- NULL
    out_est <- as.vector(sapply(raw_predictions, FUN=function(i){i$expectation_i}) %*% weights)
    names(out_se) <- NULL
    out_individual <- data.frame(est = out_est, se = out_se_i, obs = 1:length(out_est))
  }else{
    out_individual <- NULL
  }
  
  return(list(
    aggregate = out_aggregate,
    individual = out_individual,
    jacobian = jacobian_net
  ))
}