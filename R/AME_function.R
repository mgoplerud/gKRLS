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
    browser()
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
    individual = out_individual,
    jacobian = jacobian_net
  ))
}

#' @export
calculate_effects <- function(model, data = NULL, 
  variables = NULL, vcov = NULL, raw = FALSE,
  conditional = NULL, epsilon = 1e-7, verbose = FALSE,
  continuous_type = c('IQR', 'minmax', 'derivative', 'onesd')){
  
  if (!is.list(continuous_type)){
    continuous_type <- match.arg(continuous_type)
  }
  individual <- FALSE

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
            f <- function(x, m, s){s}
          }else if (continuous_type == 'derivative'){
            # Closely adapted from "margins" by Thomas Leeper
            step <- max(abs(data[[v_i]]), 1, na.rm=T) * epsilon
            multiplier <- 1/(2 * step)
            step2 <- c('0' = -step, '1' = step)
            step <- NULL
            f <- function(x, s1, s2){
              x + s2
            }
          }else if (continuous_type == 'onesd'){
            step <- mean(data[[v_i]], na.rm=T)
            sd_i <- sd(data[[v_i]], na.rm=T)
            step2 <- c('0' = -sd_i, '1' = sd_i)
            f <- function(x, m, s){
              m + s
            }
          }else if (continuous_type == 'minmax'){
            r_i <- range(data[[v_i]], na.rm=T)
            names(r_i) <- NULL
            step2 <- c('0' = r_i[1], '1' = r_i[2])
            step <- NULL
            f <- function(x, m, s){s}
          }else if (continuous_type == 'IQR'){
            q_i <- quantile(data[[v_i]], c(0.25, 0.75), na.rm=T)
            names(q_i) <- NULL
            step2 <- c('0' = q_i[1], '1' = q_i[2])
            step <- NULL
            f <- function(x, m, s){s}
          }else{stop('invalid continuous_type')}
          if (v_id == 1){
            packaged_data <- list(list(
              data = list('d1' = data, 'd0' = data),
              weights = c(1, -1) * multiplier,
              raw = raw,
              type = continuous_type,
              name = v_i
            ))
            packaged_data[[1]]$data$d0[[v_i]] <- f(packaged_data[[1]]$data$d0[[v_i]], step, step2["0"])
            packaged_data[[1]]$data$d1[[v_i]] <- f(packaged_data[[1]]$data$d1[[v_i]], step, step2["1"])
            names(packaged_data) <- v_i
          }else{
            old_names <- names(packaged_data)
            packaged_data <- lapply(packaged_data, FUN=function(i){
              i$data$d0[[v_i]] <- f(i$data$d0[[v_i]], step, step2["0"])
              i$data$d1[[v_i]] <- f(i$data$d1[[v_i]], step, step2["1"])
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
            browser()
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
    
    if ('...id' %in% colnames(out_mfx_i)){
      id <- out_mfx_i[['...id']]
    }else{
      id <- NULL
    }
    out_mfx_i <- out_mfx_i[,c('variable', 'type', 'est', 'se')]
    out_mfx_i[['...id']] <- id
    rownames(out_mfx_i) <- NULL  
    if (!is.null(conditional)){
      for (j in names(cond_data_i)){
        out_mfx_i[[j]] <- cond_data_i[[j]]
      }
    }
    out_counter <- c(out_counter, rep(cond_i, nrow(out_mfx_i)))
    out_mfx <- rbind(out_mfx, out_mfx_i)
    out_jacobian_i <- do.call('cbind', lapply(data_out, FUN=function(i){i$jacobian}))
    out_jacobian <- cbind(out_jacobian, out_jacobian_i)
  }
  cat('\n')
  if (ncol(out_jacobian) != nrow(out_mfx)){
    stop('Unusual alignment error between jacobian and marginal effects.')
  }
  out_mfx$t <- out_mfx$est/out_mfx$se
  
  out <- list(marginal_effects = out_mfx, 
              jacobian = out_jacobian,
              counter = out_counter)
  class(out) <- 'gKRLS_mfx'
  return(out)
}

#' @export
print.gKRLS_mfx <- function(x){
  if (!is.null(x$mfx_type)){
    print(x$AME_pointwise)
  }else{
    print(x$marginal_effects)
  }
}

#' @export
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