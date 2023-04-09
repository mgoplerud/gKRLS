#' Marginal Effects
#'
#' This function calculates the marginal effects after estimating a model with
#' \code{gam} or \code{bam}. For continuous predictors, a numerical
#' approximation of the partial derivative is available following Leeper (2016).
#'
#' @name calculate_effects
#' @param model A model estimated from \code{mgcv}.
#' @param object A model estimated from \code{mgcv}.
#' @param data A new data frame that used to calculate the marginal effect, or
#'   set to \code{NULL}, which the data used to estimate the model will be used.
#'   The default is \code{NULL}.
#' @param variables A character vector that specifies the variables for which to
#'   calculate effects. The default, \code{NULL}, calculates marginal effects
#'   for all variables.
#' @param vcov A matrix that specifies the covariance matrix on the parameters.
#'   The default, \code{NULL}, uses the standard covariance matrix from
#'   \code{mgcv}. This can be used to specify clustered or robust matrices
#'   using, e.g., the \code{sandwich} package.
#' @param raw Argument used for internal functions only. Default is \code{FALSE}.
#' @param individual A value of \code{TRUE} calculates individual effects (i.e.
#'   an effect for each observation in the provided data). The default is
#'   \code{FALSE}.
#' @param conditional This is an analogue of Stata's \code{at()} option and the
#'   \code{at} argument in the \code{margins}  package. For a marginal effect on
#'   some variable \code{"a"}, this can specify the values for other covariates,
#'   e.g. \code{"b"}, to be held at. Examples are provided below. This should be
#'   either \code{NULL} (default) or a data.frame where one marginal effect (per
#'   variable) is provided for each row of \code{conditional}.
#' @param epsilon A numerical value to define the step when calculating
#'   numerical derivatives. See Leeper (2016) for details.
#' @param verbose A logical value indicates whether to report progress when
#'   calculating the marginal effects.
#' @param continuous_type A character string indicating the type of marginal
#'   effects to estimate when the variable is continuous (i.e. not binary,
#'   logical, factor, or character). Options are \code{"IQR"} (compares the
#'   variable at its 25\% and 75\% percentile), \code{"minmax"} (compares the
#'   variable at its minimum and maximum), \code{"derivative"} (numerically
#'   approximates the derivative at each observed value),
#'   \code{"second_derivative"} (numerically approximates the second derivative
#'   at each observed value), \code{"onesd"} (compares one standard deviation
#'   below and one standard deviation above). It may also accepted a \bold{named
#'   list} where each named element corresponds to a numeric variable and has a
#'   two-length vector as each element. The two values are then compared.
#'   
#'   A special option (\code{"predict"}) produces predictions (e.g.,
#'   \code{predict(model, ..., type = "response")}) at each observed value and
#'   then averages them together. This in conjunction with \code{conditional}
#'   provides a way of calculating, for example, predicted probability curves
#'   using an "observed value" approach (e.g., Hanmer and Kalkan 2013). The
#'   examples below provide an illustration.
#'
#' @return  \code{calculate_effects} returns a data.frame of class
#'   \code{"gKRLS_mfx"} that contains a data.frame with the estimated average
#'   marginal effects. \code{"type"} reports the type of marginal effect
#'   calculated. For families with multiple predicted outcomes (e.g.,
#'   multinomial), the column \code{"response"} numbers the different outcomes
#'   in the same order as \code{predict.gam(object)} for the specified family.
#'   Many (but not all) extended and generalized families from \code{mgcv} are
#'   included.
#'   
#'   To estimated predicted values at different covariate strata (i.e., create
#'   an "observed value" predicted probability curve for a logistic regression),
#'   please use the \code{conditional} argument while setting
#'   \code{continuous_type = "predict"}. Please see the examples for an
#'   illustration.
#'   
#'   This data.frame has attributes, listed below, that are used in other
#'   functions. The function \code{get_individual_effects} extracts the
#'   individual-level effects that are estimated if \code{individual=TRUE}.
#'   
#'   \itemize{
#'   \item{"jacobian": } This reports the corresponding Jacobian used to
#'   calculate the standard error (via the delta method) for the estimate. There
#'   is one row for each row in "marginal_effects". This can be used to, for
#'   example, calculate a standard error on the difference between two estimated
#'   marginal effects.
#'   \item{"counter": } A placeholder for the number of marginal effects calculated.
#'   \item{"N_eff": The number of observations (in the estimation data) minus
#'   the effective degrees of freedom. This is used when calculating p-values as
#'   the degrees of freedom for the t-distribution.}
#'   \item{"N": The number of observations.}
#'   }
#'  
#'  \code{calculate_interactions} provides some simple functions for calculate
#'  interaction effects between variables. The default quantities it can produce
#'  are listed below. Egami and Imai (2019) provide a detailed exposition of
#'  these quantities. All marginalization is done using an "observed value"
#'  approach, i.e. over the estimation data or a custom dataset provided to
#'  \code{data}.
#'  \itemize{
#'    \item{"AME" or Average Marginal Effect: } This is the standard quantity
#'    reported from \code{calculate_effects}.
#'    \item{"ACE" or Average Combination Effect: } This is the effect of
#'    changing two variables simultaneously on the outcome.
#'    \item{"AMIE" or Average Marginal Interaction Effect: } This is ACE minus
#'    each corresponding AME.
#'    \item{"AIE" or Average Interaction Effect:} This has a "conditional
#'    effect" interpretation and reports the difference in average effect of one
#'    variable ("A") between two different levels of a second variable ("B").
#'  }
#'  
#' @references 
#' 
#' Egami, Naoki and Kosuke Imai. 2019. "Causal Interaction in Factorial
#' Experiments: Application to Conjoint Analysis." \emph{Journal of the American
#' Statistical Association}. 114(526):529-540.
#' 
#' Hanmer, Michael J. and Kerem Ozan Kalkan. 2013. "Behind the Curve: Clarifying
#' the Best Approach to Calculating Predicted Probabilities and Marginal Effects
#' from Limited Dependent Variable Models." \emph{American Journal of Political
#' Science} 57(1): 263-277.
#' 
#' Leeper, Thomas J. 2016. "Interpreting Regression Results using Average
#' Marginal Effects with R's \code{margins}." Working paper available at
#' \url{https://s3.us-east-2.amazonaws.com/tjl-sharing/assets/AverageMarginalEffects.pdf}.
#' 
#' @examples
#' n <- 50
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' state <- sample(letters[1:5], n, replace = TRUE)
#' y <- 0.3 * x1 + 0.4 * x2 + 0.5 * x3 + rnorm(n)
#' data <- data.frame(y, x1, x2, x3, state)
#'
#' # Make character variables into factors for mgcv
#' data$state <- factor(data$state)
#'
#' # A gKRLS model
#' gkrls_est <- mgcv::gam(y ~ state + s(x1, x2, x3, bs = "gKRLS"), data = data)
#'
#' # calculate marginal effect using derivative
#' calculate_effects(gkrls_est, variables = "x1", continuous_type = "derivative")
#'
#' # calculate marginal effect by specifying conditional variables
#' calculate_effects(gkrls_est,
#'   variables = "x1",
#'   conditional = data.frame(x2 = c(0.6, 0.8), x3 = 0.3)
#' )
#'
#' # calculate interaction effects between two variables
#' # use the default setting ("IQR") for the baseline and
#' # comparison categories for each variable
#' calculate_interactions(gkrls_est, 
#'    variables = list(c("x1", "x2")),
#'    QOI = c('AIE', 'AMIE')
#' )
#' 
#' # calculate marginal effect by specifying a factor conditional variable
#' # estimate the individual marginal effects
#' out <- calculate_effects(gkrls_est,
#'   variables = "x1", individual = TRUE,
#'   conditional = data.frame(state = c("a", "b", "c")), continuous_type = "derivative"
#' )
#' 
#' # Extract the individual marginal effects: 
#' # shorthand for attr(fit_main, 'individual')
#' get_individual_effects(out)
#' 
#' # calculated the average expected value across a grid of "x1" using an observed value approach
#' calculate_effects(gkrls_est, conditional = data.frame(x1 = c(0, 0.2, 0.4, 0.6)),
#'   continuous_type = 'predict'
#' )
#' @importFrom stats model.frame sd
#' @export
calculate_effects <- function(model, data = NULL,
    variables = NULL, vcov = NULL, raw = FALSE, individual = FALSE,
    conditional = NULL, epsilon = 1e-7, verbose = FALSE,
    continuous_type = c("IQR", "minmax", "derivative", 
      "onesd", "predict", "second_derivative")) {
  if (!is.list(continuous_type)) {
    continuous_type <- match.arg(continuous_type)
  }

  simple_family <- !inherits(model$family, c('general.family', 'extended.family'))
  N_eff <- length(model$y) - sum(model$edf)
  N <- length(model$y)

  if (!all(model$prior.weights == 1)){
    warning('calculate_effects ignores "weights" argument when calculating effects.')
  }
  if (!is.null(model$offset)){
    flat_offset <- unlist(model$offset)
    if (!all(flat_offset == 0)){
      stop('"offset" not set up for calculate_effects.')
    }
  }
  if (is.null(data)) {
    raw_data <- model.frame(model)
  } else {
    raw_data <- data
    rm(data)
  }

  if (identical(continuous_type, 'predict')){
    if (raw){
      raw <- FALSE
      message("raw = FALSE if continuous_type = 'predict'")
    }
  }
  if (is.list(continuous_type)){
    fmt_type <- 'custom'
  }else{
    fmt_type <- continuous_type
  }

  if (is.null(vcov)) {
    vcov <- stats::vcov(model)
  }

  all_variables <- gam_terms(model = model, variables = NULL)
  variable_list <- unlist(all_variables)
  variable_type <- rep(names(all_variables), lengths(all_variables))
  names(variable_type) <- unlist(all_variables)

  # Address the conditional variables
  if (!is.null(conditional)) {
    names_cond <- names(conditional)
    # Check the variables in conditional are valid
    if ((!all(names_cond %in% variable_list)) & (length(setdiff(names_cond, variable_list)) != 0)) {
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
  } else {
    ncond <- 1
  }

  any_binary <- sapply(names(variable_type), FUN=function(i){
    j <- raw_data[[i]]
    all(j %in% c(0,1)) & !all(is.logical(j))
  })
  variable_type[which(any_binary)] <- "bnames"

  if (identical(continuous_type, 'predict')){
    variable_list <- list('...placeholder')
  }else if (is.null(variables)) {
    variable_list <- as.list(variable_list)
  } else if (is.list(variables)) {
    check_valid <- all(sapply(variables, FUN = function(i) {
      all(i %in% variable_list)
    }))
    if (!check_valid) {
      stop('All elements of list of "variables" must be in model.')
    }
    check_valid <- all(sapply(variables, anyDuplicated) == 0)
    if (!check_valid) {
      stop('All elements of a list of "variables" must be distinct.')
    }
    variable_list <- variables
  } else if (is.vector(variables)) {
    if (!all(variables %in% variable_list)) {
      stop('All elements of "variables" must be in model.')
    }
    if (anyDuplicated(variables)) {
      stop('All elements of "variables" must be unique.')
    }
    variable_list <- as.list(variables)
  } else {
    stop('"variables" must be NULL, list or a vector.')
  }
  # "conditional" variables cannot be used in "variables"
  
  check_overlap <- intersect(unlist(variable_list), colnames(conditional))
  if (length(check_overlap) > 0){
    stop('"conditional" should not contain any variables in "variable" argument.')
  }
  
  out_mfx <- data.frame()
  if (individual) {
    out_mfx_individual <- data.frame()
  }else{
    out_mfx_individual <- NULL
  }
  
  if (simple_family){
    out_jacobian <- matrix(nrow = length(coef(model)), ncol = 0)
  }else{
    out_jacobian <- list()
  }
  out_counter <- c()

  if (is.list(continuous_type)) {
    v <- unlist(variable_list)
    v <- variable_type[match(v, names(variable_type))]
    cont_v <- v[v == "nnames"]
    # is a continous variable missing from "continuous_type"
    name_not_in <- !all(names(cont_v) %in% names(continuous_type))
    # invalid inclusion of non-continous variable
    invalid_inclusion <- any(names(v[v != 'nnames']) %in% names(continuous_type))
    if (name_not_in) {
      stop("if continuous_type is a list, then it must have all numerical variables listed in it.")
    }
    if (invalid_inclusion) {
      stop("if continuous_type is a list, then it must not include any non-continuous variable (including binary or logical).")
    }
    if (!all(lengths(continuous_type) == 2)) {
      stop("if continous_type is a list, it must have vectors of two elements.")
    }
    continuous_type <- lapply(continuous_type, FUN=function(i){
      names(i) <- NULL
      return(i)
    })
  }

  continuous_f <- function(x, m, s, type) {
    if (type == "list") {
      s
    } else if (type == "derivative") {
      x + s
    } else if (type == "onesd") {
      m + s
    } else if (type %in% c("minmax", "IQR")) {
      s
    } else {
      NA
    }
  }
  for (cond_i in seq_len(ncond)) {
    if (verbose & cond_i %% ceiling(ncond / ncond) == 0) {
      message("|")
    }
    data <- raw_data

    if (!is.null(conditional)) {
      cond_data_i <- conditional[cond_i, , drop = F]
      for (j in names(cond_data_i)) {
        data[[j]] <- cond_data_i[[j]]
      }
    }

    data_out <- list()
    for (v in variable_list) {
      
      type_v <- variable_type[v]
      packaged_data <- list()
      # if "AIE" then use a *different* type
      if (isTRUE(attr(v, 'aie'))){
        aie_estimate <- TRUE
        if (identical(continuous_type, 'predict') | identical(continuous_type, 'derivative') | identical(continuous_type, 'second_derivative')){
          stop('aie_estimate cannot be true with predict or derivative continuous type.')          
        }
        if (length(v) != 2){stop("aie_estimate requires length(v) == 2")}
      }else{
        aie_estimate <- FALSE
      }
      # Loop over each variable and create the "d0" and "d1" frames, respectively
      for (v_id in seq_len(length(v))) {
        v_i <- v[v_id]
        type_i <- type_v[v_id]
        multiplier <- 1
        if (type_i %in% c("nnames", NA)) {
          if (is.list(continuous_type)) {
            r_i <- continuous_type[[v_i]]
            step2 <- c("0" = r_i[1], "1" = r_i[2])
            step <- NULL
            ctype <- "list"
          } else if (continuous_type == "second_derivative") {
            step <- sqrt(max(abs(data[[v_i]]), 1, na.rm = T) * epsilon)
            # multiplier <- 1 / (2 * step)
            multiplier <- 1/step^2
            step2 <- c("-1" = -step, "0" = 0, "1" = step)
            step <- NULL
            ctype <- "derivative"
          } else if (continuous_type == "derivative") {
            # Closely adapted from "margins" by Thomas Leeper
            step <- max(abs(data[[v_i]]), 1, na.rm = T) * epsilon
            multiplier <- 1 / (2 * step)
            step2 <- c("0" = -step, "1" = step)
            step <- NULL
            ctype <- "derivative"
          } else if (continuous_type == "onesd") {
            step <- mean(data[[v_i]], na.rm = T)
            sd_i <- sd(data[[v_i]], na.rm = T)
            step2 <- c("0" = -sd_i, "1" = sd_i)
            ctype <- "onesd"
          } else if (continuous_type == "minmax") {
            r_i <- range(data[[v_i]], na.rm = T)
            names(r_i) <- NULL
            step2 <- c("0" = r_i[1], "1" = r_i[2])
            step <- NULL
            ctype <- "minmax"
          } else if (continuous_type == "IQR") {
            q_i <- quantile(data[[v_i]], c(0.25, 0.75), na.rm = T)
            names(q_i) <- NULL
            step2 <- c("0" = q_i[1], "1" = q_i[2])
            step <- NULL
            ctype <- "IQR"
          } else if (continuous_type == "predict") {
            ctype <- "predict"
            step <- NULL
            step2 <- NULL
          } else {
            stop("invalid continuous_type")
          }
          
          if (v_id == 1) {
            if (identical(continuous_type, "predict")){
              packaged_data <- list(list(
                data = list("d0" = data),
                weights = 1,
                raw = raw,
                type = fmt_type,
                name = v_i
              ))
            }else if (identical(continuous_type, 'second_derivative')){
              packaged_data <- list(list(
                data = list("d2" = data, "d1" = data, "d0" = data),
                weights = c(1, -2, 1) * multiplier,
                raw = raw,
                type = fmt_type,
                name = v_i
              ))
              packaged_data[[1]]$data$d0[[v_i]] <- continuous_f(packaged_data[[1]]$data$d0[[v_i]], step, step2["-1"], ctype)
              packaged_data[[1]]$data$d1[[v_i]] <- continuous_f(packaged_data[[1]]$data$d1[[v_i]], step, step2["0"], ctype)
              packaged_data[[1]]$data$d2[[v_i]] <- continuous_f(packaged_data[[1]]$data$d2[[v_i]], step, step2["1"], ctype)
            }else{
              packaged_data <- list(list(
                data = list("d1" = data, "d0" = data),
                weights = c(1, -1) * multiplier,
                raw = raw,
                type = fmt_type,
                name = v_i
              ))
              packaged_data[[1]]$data$d0[[v_i]] <- continuous_f(packaged_data[[1]]$data$d0[[v_i]], step, step2["0"], ctype)
              packaged_data[[1]]$data$d1[[v_i]] <- continuous_f(packaged_data[[1]]$data$d1[[v_i]], step, step2["1"], ctype)
            }
            names(packaged_data) <- v_i
          } else {
            
            old_names <- names(packaged_data)
            packaged_data <- lapply(packaged_data, FUN = function(i) {
              
              if (aie_estimate){
                # Get all four relevant terms if "aie_estimate" = TRUE
                
                d00 <- d01 <- i$data$d0
                d11 <- d10 <- i$data$d1
                
                d00[[v_i]] <- continuous_f(d00[[v_i]], step, step2["0"], ctype)
                d01[[v_i]] <- continuous_f(d01[[v_i]], step, step2["1"], ctype)
                d10[[v_i]] <- continuous_f(d10[[v_i]], step, step2["0"], ctype)
                d11[[v_i]] <- continuous_f(d11[[v_i]], step, step2["1"], ctype)
                
                i$data <- list(d00 = d00, d01 = d01, d10 = d10, d11 = d11)
                stopifnot(!any(i$type %in% "derivative"))
                i$weights <- c(1, -1, -1, 1)
                stopifnot(!any(i$raw == TRUE))
                i$type <- c(i$type, fmt_type)
                i$name <- c(paste0('%%AIE%%', i$name), v_i)
              }else{
                i$data$d0[[v_i]] <- continuous_f(i$data$d0[[v_i]], step, step2["0"], ctype)
                i$data$d1[[v_i]] <- continuous_f(i$data$d1[[v_i]], step, step2["1"], ctype)
                if (!any(i$type %in% "derivative")) {
                  i$weights <- i$weights * multiplier
                }
                i$raw <- i$raw * raw
                i$type <- c(i$type, fmt_type)
                i$name <- c(i$name, v_i)
              }
              return(i)
            })
            names(packaged_data) <- paste(names(packaged_data), ":", v_i, sep = "")
          }
        } else if (type_i == "bnames") {
          if (v_id == 1) {
            packaged_data <- list(list(
              data = list("d1" = data, "d0" = data),
              weights = c(1, -1),
              raw = raw,
              type = "binary",
              name = v_i
            ))
            packaged_data[[1]]$data$d0[[v_i]] <- 0
            packaged_data[[1]]$data$d1[[v_i]] <- 1
            names(packaged_data) <- v_i
          } else {
            
            old_names <- names(packaged_data)
            packaged_data <- lapply(packaged_data, FUN = function(i) {
              
              if (aie_estimate){
                # Get all four relevant terms if "aie_estimate" = TRUE
                
                d00 <- d01 <- i$data$d0
                d11 <- d10 <- i$data$d1
                
                d00[[v_i]] <- 0
                d01[[v_i]] <- 1
                d10[[v_i]] <- 0
                d11[[v_i]] <- 1
                
                i$data <- list(d00 = d00, d01 = d01, d10 = d10, d11 = d11)
                i$weights <- c(1, -1, -1, 1)
                stopifnot(!any(i$raw == TRUE))
                i$type <- c(i$type, "binary")
                i$name <- c(paste0('%%AIE%%', i$name), v_i)
                
              }else{
                i$data$d0[[v_i]] <- 0
                i$data$d1[[v_i]] <- 1
                i$weights <- i$weights
                i$raw <- i$raw * raw
                i$type <- c(i$type, "binary")
                i$name <- c(i$name, v_i)
                
              }
              return(i)
            })
            names(packaged_data) <- paste(names(packaged_data), ":", v_i, sep = "")
          }
        } else if (type_i == "lnames") {
          if (v_id == 1) {
            packaged_data <- list(list(
              data = list("d1" = data, "d0" = data),
              weights = c(1, -1),
              raw = raw,
              type = "logical",
              name = v_i
            ))
            packaged_data[[1]]$data$d0[[v_i]] <- FALSE
            packaged_data[[1]]$data$d1[[v_i]] <- TRUE
            names(packaged_data) <- v_i
          } else {
            
            old_names <- names(packaged_data)
            packaged_data <- lapply(packaged_data, FUN = function(i) {
              
              if (aie_estimate){
                
                d00 <- d01 <- i$data$d0
                d11 <- d10 <- i$data$d1
                
                d00[[v_i]] <- FALSE
                d01[[v_i]] <- TRUE
                d10[[v_i]] <- FALSE
                d11[[v_i]] <- TRUE
                
                i$data <- list(d00 = d00, d01 = d01, d10 = d10, d11 = d11)
                i$weights <- c(1, -1, -1, 1)
                stopifnot(!any(i$raw == TRUE))
                i$type <- c(i$type, "logical")
                i$name <- c(paste0('%%AIE%%', i$name), v_i)
                
              }else{
                
                i$data$d0[[v_i]] <- FALSE
                i$data$d1[[v_i]] <- TRUE
                i$weights <- i$weights
                i$raw <- i$raw * raw
                i$type <- c(i$type, "logical")
                i$name <- c(i$name, v_i)
                
              }
              return(i)
            })
            names(packaged_data) <- paste(names(packaged_data), ":", v_i, sep = "")
          }
        } else if (type_i == "fnames") {
          levs <- levels(as.factor(raw_data[[v_i]]))
          base <- levs[1L]
          levs <- levs[-1L]
          if (v_id == 1) {
            packaged_data <- list()
            for (i in seq_along(levs)) {
              tmp_name <- paste0("factor(", v_i, ")", levs[i])
              packaged_data[[tmp_name]] <- list(
                data = list("d1" = data, "d0" = data),
                weights = c(1, -1),
                raw = raw,
                type = "factor",
                name = tmp_name
              )
              packaged_data[[tmp_name]]$data$d0[[v_i]] <- base
              packaged_data[[tmp_name]]$data$d1[[v_i]] <- levs[i]
            }
          } else {
            old_names <- names(packaged_data)
            added_data <- list()
            for (l in seq_along(levs)) {
              temp_name <- paste0("factor(", v_i, ")", levs[l])
              add_l <- mapply(packaged_data, names(packaged_data), SIMPLIFY = FALSE, FUN = function(i, m) {
                if (aie_estimate){
                  d00 <- d01 <- i$data$d0
                  d11 <- d10 <- i$data$d1
                  d00[[v_i]] <- base
                  d01[[v_i]] <- levs[l]
                  d10[[v_i]] <- base
                  d11[[v_i]] <- levs[l]
                  
                  i$data <- list(d00 = d00, d01 = d01, d10 = d10, d11 = d11)
                  i$weights <- c(1, -1, -1, 1)
                  stopifnot(!any(i$raw == TRUE))
                  i$type <- c(i$type, 'factor')
                  i$name <- c(paste0('%%AIE%%', temp_name), v_i)
                }else{
                  i$data$d0[[v_i]] <- base
                  i$data$d1[[v_i]] <- levs[l]
                  i$weights <- i$weights
                  i$raw <- i$raw * raw
                  i$type <- c(i$type, "factor")
                  i$name <- c(i$name, temp_name)
                }
                return(i)
              })
              names(add_l) <- paste0(names(packaged_data), ":", temp_name)
              added_data <- c(added_data, add_l)
            }
            packaged_data <- added_data
          }
        } else {
          stop("Unknown type")
        }
      }
      
      gc()

      fit_mfx <- lapply(packaged_data, FUN = function(data_i) {
        fit_i <- weighted_mfx(
          model = model,
          data_list = data_i$data,
          weights = list("AME" = data_i$weights),
          vcov = vcov,
          individual = individual, raw = data_i$raw
        )
        if (data_i$raw) {
          fit_i$aggregate[["...id"]] <- c("effect", "raw_1", "raw_0")
        }
        fit_i$name <- paste(data_i$name, collapse = ":")
        fit_i$type <- paste(data_i$type, collapse = ":")
        return(fit_i)
      })
      data_out <- c(data_out, fit_mfx)
    }

    out_mfx_i <- do.call("rbind", lapply(data_out, FUN = function(i) {
      i$aggregate$variable <- i$name
      i$aggregate$name <- NULL
      i$aggregate$type <- i$type
      return(i$aggregate)
    }))

    if (individual) {
      out_mfx_i_ind <- do.call("rbind", lapply(data_out, FUN = function(i) {
        i$individual$variable <- i$name
        i$individual$name <- NULL
        i$individual$type <- i$type
        return(i$individual)
      }))
      rownames(out_mfx_i_ind) <- NULL
    }

    if ("...id" %in% colnames(out_mfx_i)) {
      id <- out_mfx_i[["...id"]]
      if (individual) {
        id_individual <- out_mfx_i_ind[["...id"]]
      }
    } else {
      id <- NULL
      id_individual <- NULL
    }
    select_col <- c("variable", "type", "est", "se")
    if ('response' %in% names(out_mfx_i)){
      select_col <- c(select_col, 'response')
    }
    out_mfx_i <- out_mfx_i[ , select_col]
    out_mfx_i[["...id"]] <- id
    rownames(out_mfx_i) <- NULL

    if (individual) {
      select_col <- c("obs", "variable", "type", "est", "se")
      if ('response' %in% names(out_mfx_i_ind)){
        select_col <- c(select_col, 'response')
      }
      out_mfx_i_ind <- out_mfx_i_ind[, select_col]
      out_mfx_i_ind[["...id"]] <- id
    }

    if (!is.null(conditional)) {
      for (j in names(cond_data_i)) {
        out_mfx_i[[j]] <- cond_data_i[[j]]
        if (individual) {
          out_mfx_i_ind[[j]] <- cond_data_i[[j]]
        }
      }
    }

    out_counter <- c(out_counter, rep(cond_i, nrow(out_mfx_i)))
    out_mfx <- rbind(out_mfx, out_mfx_i)
    
    if (simple_family){
      out_jacobian_i <- do.call("cbind", lapply(data_out, FUN = function(i) {
        i$jacobian
      }))
      out_jacobian <- cbind(out_jacobian, out_jacobian_i)
    }else{

      out_jacobian_i <- lapply(data_out, FUN = function(i) {
        i$jacobian
      })
      out_jacobian <- c(out_jacobian, list(out_jacobian_i))
      
    }

    if (individual) {
      out_mfx_individual <- rbind(out_mfx_individual, out_mfx_i_ind)
    }

  }
  
  if (simple_family){
    if (ncol(out_jacobian) != nrow(out_mfx)) {
      stop("Unusual alignment error between jacobian and marginal effects.")
    }
  }else{
    # out_jacobian: One for each "conditional"
    # one for each variable
    checksum_jacobian <- sum(sapply(out_jacobian, FUN=function(i){
      sum(sapply(i, FUN=function(j){sapply(j, ncol)}))
    }))
    if (checksum_jacobian != nrow(out_mfx)){
      stop("Unusual alignment error between jacobian and marginal effects.")
    }
  }
  
  out_mfx$t <- out_mfx$est / out_mfx$se
  out_mfx$p.value <- 2 * pt(-abs(out_mfx$t), df = N_eff)

  if (individual) {
    out_mfx_individual$t <- out_mfx_individual$est / out_mfx_individual$se
    out_mfx_individual$p.value <- 2 * pt(-abs(out_mfx_individual$t), df = N_eff)
  }

  if (identical(continuous_type, 'predict')){
    if (individual){
      out_mfx_individual <- out_mfx_individual[, !(names(out_mfx_individual) %in% c('variable'))]
    }
    out_mfx <- out_mfx[, !(names(out_mfx) %in% c('variable'))]
  }

  out <- out_mfx
  attr(out, 'jacobian') <- out_jacobian
  attr(out, 'counter') <- out_counter
  attr(out, 'individual') <- out_mfx_individual
  attr(out, 'N_eff') <- N_eff
  attr(out, 'N') <- N
  class(out) <- c("gKRLS_mfx", "data.frame")
  return(out)
}

#' @rdname calculate_effects
#' @param QOI Quantity of interest to calculate for
#'   \code{calculate_interactions}. Options include AME (average marginal
#'   effect), ACE (average combination effect), AIE (average interaction effect)
#'   and AMIE (average marginal interaction effect); see "Details" for more
#'   information.
#' @param ... Used for \code{calculate_interactions} to pass arguments to
#'   \code{calculate_effects}.
#' @export
calculate_interactions <- function(model,
      variables, QOI = c("AMIE", "ACE", "AME", "AIE"), ...) {

  QOI <- match.arg(QOI, several.ok = TRUE)
  args <- list(...)

  if (isTRUE(args$individual)){
    stop("individual=TRUE not set up for calculate_interactions.")
  }
  
  if (!is.list(args$continuous_type)){
    if (isTRUE(args$continuous_type == 'predict')){
      stop('continous_type="predict" not set up for calculate_interactions.')
    }  
    if (isTRUE(grepl(args$continuous_type, pattern='derivative'))){
      stop('continous_type="derivative" or "second_derivative" not set up for calculate_interactions.')
    } 
  }
  
  if (isTRUE(args$raw)) {
    stop("raw=T not permitted for interactions.")
  }
  
  if (!is.null(args$conditional)) {
    if (any(names(args$conditional) %in% c("variable", "est", "se", "QOI", "...counter"))) {
      stop('conditional may not contain "variable", "est", "se", or "QOI" as column names.')
    }
  }
  
  if (any(lengths(variables) != 2)) {
    stop("variables must be a list of two-length vectors.")
  }
  
  if (is.null(args$vcov)) {
    vcov_mdl <- vcov(model)
  } else {
    vcov_mdl <- args$vcov
  }

  all_inter <- data.frame()
  for (v in variables) {
    
    fmt_v <- list()
    if (any(c('AME', 'ACE', 'AMIE') %in% QOI)){
      fmt_v <- c(fmt_v, as.list(v))
      fmt_v <- c(fmt_v, list(v))
    }
    if ('AIE' %in% QOI){
      aie_list <- v
      attr(aie_list, 'aie') <- TRUE
      fmt_v <- c(fmt_v, list(aie_list))
    }

    fit_var <- do.call(
      "calculate_effects",
      c(
        list(model = model, variables = fmt_v),
        args
      )
    )
    
    if (!('response' %in% names(fit_var))){
      fit_var$response <- 1
    }
    fit_var$counter <- paste(attr(fit_var, 'counter'), '@', fit_var$response)
    
    id <- seq_len(nrow(fit_var))
    split_var <- strsplit(fit_var$variable, split = ":")
    split_var <- unique(split_var[which(lengths(split_var) == 2)])

    out_inter <- lapply(split_var, FUN = function(i) {
      
      out <- data.frame()
      v_i <- paste(i, collapse = ":")
      is_aie <- grepl(v_i, pattern='^%%AIE%%')
      
      id_inter <- which(fit_var$variable == v_i)
      i <- gsub(i, pattern='^%%AIE%%', replacement='')
      
      id_marg <- which(fit_var$variable %in% i)
      id_marg <- split(id_marg, fit_var$counter[id_marg])
      id_inter <- split(id_inter, fit_var$counter[id_inter])
      fmt_names <- do.call('rbind', strsplit(names(id_inter), split=' @ '))
      colnames(fmt_names) <- c('counter', 'response')
      stopifnot(all(lengths(id_inter) == 1))
      stopifnot(all(lengths(id_marg) == 2))
      if (any(c('AME', 'AMIE', 'ACE') %in% QOI)){
        stopifnot(identical(names(id_inter), names(id_marg)))
      }
      
      if ( ("AME" %in% QOI) & !is_aie ) {
        out <- rbind(
          out,
          data.frame(
            fit_var[unlist(id_marg), c("variable", "est", "se", "type")],
            QOI = "AME",
            `...counter` = rep(seq_len(length(id_marg)), lengths(id_marg))
          )
        )
      }
      if ( ("ACE" %in% QOI) & !is_aie ) {
        out <- rbind(
          out,
          data.frame(
            fit_var[unlist(id_inter), c("variable", "est", "se", "type")],
            QOI = "ACE", 
            `...counter` = seq_len(length(unlist(id_inter)))
          )
        )
      }
      
      if ( ("AIE" %in% QOI) & is_aie){
        
        extract_AIE <- data.frame(
          fit_var[unlist(id_inter), c("variable", "est", "se", "type")],
          QOI = "AIE", 
          `...counter` = seq_len(length(unlist(id_inter)))
        )
        extract_AIE$variable <- gsub(extract_AIE$variable, pattern='^%%AIE%%', replacement = '')
        out <- rbind(out, extract_AIE)
        
      }
      
      if ( ("AMIE" %in% QOI) & !is_aie ) {
        
        id_matrix <- do.call("rbind", mapply(id_inter, id_marg, seq_len(length(id_marg)), SIMPLIFY = FALSE, FUN = function(x, y, z) {
          rbind(c(x, z, 1), cbind(y, z, -1))
        }))
        id_matrix <- sparseMatrix(
          i = id_matrix[, 1], j = id_matrix[, 2], x = id_matrix[, 3],
          dims = c(nrow(fit_var), length(id_marg))
        )
        rownames(id_matrix) <- paste(fit_var$variable, '@', fit_var$response)
        
        if (all(fit_var$response == 1)){
          point_estimate <- as.vector(fit_var$est %*% id_matrix)
          se_estimate <- apply(attr(fit_var, 'jacobian') %*% id_matrix,
             MARGIN = 2, FUN = function(i) {
               as.numeric(sqrt(t(i) %*% vcov_mdl %*% i))
             }
          )
        }else{
          point_estimate <- as.vector(fit_var$est %*% id_matrix)
          # Loop over the jacobian for each conditional group
          se_estimate <- do.call('c', lapply(attr(fit_var, 'jacobian'), FUN=function(ji){
            # Loop over the jacobian for each *response*
            se_loop <- mapply(ji[[i[1]]], ji[[i[2]]], ji[[v_i]], SIMPLIFY = TRUE, FUN=function(marg_1, marg_2, inter){
               jacobian_ji <- cbind(inter, marg_1, marg_2) %*% c(1, -1, -1)
               out <- sqrt(as.numeric(t(jacobian_ji) %*% vcov_mdl %*% jacobian_ji))
               return(out)
            })
            return(se_loop)
          }))
        }
        extract_AMIE <- data.frame(variable = v_i, QOI = "AMIE",
         type = fit_var$type[unlist(id_inter)], 
         est = point_estimate, se = se_estimate,
         `...counter` = seq_len(length(point_estimate)),
         stringsAsFactors = F)
        
        out <- rbind(out, extract_AMIE)
        
      }
      out$response <- as.numeric(fmt_names[,'response'])[out$`...counter`]
      out$`...counter` <- as.numeric(fmt_names[,'counter'])[out$`...counter`]
      
      return(out)
    })
    out_inter <- do.call("rbind", out_inter)
    all_inter <- rbind(all_inter, out_inter)
  }

  rownames(all_inter) <- NULL
  all_inter <- unique(all_inter)
  all_inter$t <- all_inter$est / all_inter$se
  if ("...id" %in% names(all_inter)) {
    id <- all_inter[["...id"]]
  } else {
    id <- NULL
  }

  all_inter <- all_inter[, c("QOI", "type", "variable", "est", "se", "t", "response", "...counter")]
  all_inter[["...id"]] <- id
  if (!is.null(args$conditional)) {
    conditional <- args$conditional
    for (v in colnames(conditional)) {
      all_inter[[v]] <- conditional[[v]][all_inter[["...counter"]]]
    }
  }
  all_inter[["...counter"]] <- NULL
  out <- all_inter
  class(out) <- c("gKRLS_mfx", "data.frame")
  return(out)
}

# Adapted from Leeper's "margins"
gam_terms <- function(model, variables = NULL) {
  classes <- attributes(terms(model))$dataClasses[-1]
  if (is.null(classes)) {
    stop("terms not found from gam")
  }
  classes <- classes[!names(classes) %in% "(weights)"]
  classes[classes == "character"] <- "factor"
  # List of types, again following scheme in Leeper's "margins"
  vars <- list(
    nnames = unique(names(classes)[!classes %in% c("factor", "ordered", "logical")]),
    lnames = unique(names(classes)[classes == "logical"]),
    fnames = unique(names(classes)[classes %in% c("factor", "ordered")])
  )
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
#' @export
get_individual_effects <- function(x){
  return(attr(x, 'individual'))
}

#' @rdname calculate_effects
#' @param x Object fit with \code{calculate_effects} or \code{legacy_marginal_effect}
#' @method print gKRLS_mfx
#' @importFrom stats quantile pt
#' @export
print.gKRLS_mfx <- function(x, ...) {
  if (!is.null(x$mfx_type)) {
    summary_pointwise <- apply(x$ME_pointwise, MARGIN = 2, FUN = function(i) {
      quantile(i, c(0.25, 0.5, 0.75))
    })
    cat(paste0("Distribution of Pointwise Marginal Effects: N = ", nrow(x$ME_pointwise), "\n"))
    print(summary_pointwise)

    out <- data.frame(est = x$AME_pointwise, se = sqrt(x$AME_pointwise_var))

    out$t.stat <- out$est / out$se
    out$p.value <- 2 * pt(-abs(out$t.stat), df = x$N_eff)
    cat("\nSummary of Average Marginal Effects\n")
    print(out)
  } else {
    print.data.frame(x)
  }
}


#' @rdname calculate_effects
#' @method summary gKRLS_mfx
#' @param ... Additional arguments (unused).
#' @export
summary.gKRLS_mfx <- function(object, ...) {
  if (!is.null(object$type)){
    out <- data.frame(est = object$AME_pointwise, se = sqrt(object$AME_pointwise_var))
    out$t.stat <- out$est / out$se
    out$p.value <- 2 * pt(-abs(out$t.stat), df = object$N_eff)
    out <- cbind(data.frame(variable = rownames(out), stringsAsFactors = F), out)
    rownames(out) <- NULL
  }else{
    out <- object
  }
  return(out)
}

#' @importFrom stats coef vcov na.pass predict
weighted_mfx <- function(model, data_list, vcov,
                         weights, raw = FALSE, individual = FALSE) {
  model_coef <- coef(model)

  simple_family <- !inherits(model$family, c('general.family', 'extended.family'))
  
  if (missing(vcov)) {
    vcov <- vcov(model)
  } else if (identical(vcov, "none")) {
    vcov <- NULL
  } else {
    if (nrow(vcov) != length(coef(model))) {
      stop("If vcov is provided manually, it must be the same size as the coefficient vector.")
    }
  }

  if (missing(weights)) {
    stop("weights may not be null..")
  } else {
    if (any(lengths(weights) != length(data_list))) {
      stop('"weights" must be list of vectors, each with the same length as "data_list".')
    }
  }
  weights <- do.call("cbind", weights)
  if (raw) {
    add_cols <- diag(length(data_list))
    colnames(add_cols) <- paste0("raw_", 1:length(data_list))
    weights <- cbind(weights, add_cols)
  }

  # Get the predictions for each covariate profile provided
  raw_predictions <- lapply(data_list, FUN = function(data_i) {
    # Get the design
    matrix_i <- predict(model, newdata = data_i, 
      na.action = na.pass, type = "lpmatrix")
    if (simple_family){
      lp_i <- as.vector(matrix_i %*% model_coef)
      e_i <- model$family$linkinv(lp_i)
      ex <- mean(e_i, na.rm = T)
      
      if (individual) {
        jacob_i <- Diagonal(x = model$family$mu.eta(lp_i)) %*% matrix_i
        jacob <- colMeans(jacob_i, na.rm = T)
        nrow_valid <- sum(!is.na(lp_i))
      } else {
        se_i <- NULL
        e_i <- NULL
        jacob_i <- NULL
        nrow_valid <- NULL
        jacob <- colMeans(Diagonal(x = model$family$mu.eta(lp_i)) %*% matrix_i, na.rm = T)
      }
      
      out <- list(
        expectation = ex,
        expectation_i = e_i,
        jacobian_i = jacob_i,
        jacobian = jacob,
        nrow_valid = nrow_valid
      )
    }else{
      out <- predict_extended(object = model, 
        X = matrix_i, individual = individual)
    }
    return(out)
  })
  
  if (!simple_family){
    noutcome <- length(raw_predictions[[1]]$jacobian)
    complex_outcome <- isTRUE(raw_predictions[[1]]$complex_extended)
    
    jacobian_net <- lapply(1:noutcome, FUN=function(d){
      sapply(raw_predictions, FUN = function(i) {
        i$jacobian[[d]]
      }) %*% weights
    })
    
    out_est <- sapply(1:noutcome, FUN=function(d){
      sapply(raw_predictions, FUN=function(i){
        i$expectation[[d]]
      }) %*% weights
    })
    
    if (!complex_outcome){
      # If "simple extended family", then do each outcome
      # with its own linear predictor
      lpi <- lapply(1:noutcome, FUN=function(i){
        li <- sapply(raw_predictions, FUN=function(j){j$lpi[[i]]})
        range_li <- max(abs(sweep(li, MARGIN = 1, STATS = rowMeans(li), FUN = '-')))
        if (range_li != 0){stop('...')}
        return(li[,1])
      })
    
      out_se <- sapply(1:noutcome, FUN=function(d){
        sqrt(rowSums( (t(jacobian_net[[d]]) %*% vcov[lpi[[d]], lpi[[d]]]) * t(jacobian_net[[d]]) ))
      })
      
    }else{
      # If "complex" family, e.g., multinomial, do this one.
      out_se <- sapply(1:noutcome, FUN=function(d){
        sqrt(rowSums( (t(jacobian_net[[d]]) %*% vcov) * t(jacobian_net[[d]]) ))
      })
    }
    
    if (ncol(weights) == 1){
      out_est <- t(matrix(out_est))
      out_se <- t(matrix(out_se))
    }
    
    out_aggregate <- do.call('rbind', lapply(1:noutcome, FUN=function(d){
      out_aggregate <- data.frame(
        name = colnames(weights), 
        est = out_est[,d], 
        se = out_se[,d])
      out_aggregate$response <- d
      return(out_aggregate)
    }))
    rownames(out_aggregate) <- NULL
  }else{
    # Get the jacobian/gradient for each weighted average
    jacobian_net <- sapply(raw_predictions, FUN = function(i) {
      i$jacobian
    }) %*% weights
    
    out_se <- sqrt(rowSums( (t(jacobian_net) %*% vcov) * t(jacobian_net) ))
    # out_se <- sqrt(apply(jacobian_net, MARGIN = 2, FUN = function(i) {
    #   as.vector(t(i) %*% vcov %*% i)
    # }))
    out_est <- as.numeric(sapply(raw_predictions, FUN = function(i) {
      i$expectation
    }) %*% weights)
    
    out_aggregate <- data.frame(name = colnames(weights), est = out_est, se = out_se)
    
  }

  if (individual) {
    
    checksum_ind <- length(unique(sapply(raw_predictions, FUN = function(i) {
      i$nrow_valid
    })))
    if (checksum_ind != 1) {
      stop("individual=TRUE requires same pattern of missing data across all elements of data_list.")
    }
    
    if (simple_family){
      
      extract_ei <- sapply(raw_predictions, FUN = function(i) {
        i$expectation_i
      })
      extract_jacob_i <- sapply(raw_predictions, FUN = function(i) {
        i$jacobian_i
      })
      
      out_individual <- do.call('rbind', lapply(1:ncol(weights), FUN=function(w){
        ji <- Reduce('+', mapply(extract_jacob_i, weights[,w], FUN = function(i, j) {
          i * j
        }))
        est <- as.vector(extract_ei %*% weights[,w])
        se <- sqrt(rowSums( (ji %*% vcov) * ji))
        out <- data.frame(est = est, se = se, obs = 1:length(est), variable = w)
        return(out)
      }))
      out_individual$variable <- colnames(weights)[out_individual$variable]
    } else {

      out_individual <- lapply(1:noutcome, FUN=function(d){
        
        extract_ei <- sapply(raw_predictions, FUN = function(i) {
          i$expectation_i[,d]
        })
        extract_jacob_i <- sapply(raw_predictions, FUN = function(i) {
          i$jacobian_i[[d]]
        })
        
        if (!complex_outcome){
          lpi_d <- lpi[[d]]
          vcov_d <- vcov[lpi_d, lpi_d]
        }else{
          vcov_d <- vcov
        }
        
        out_individual <- do.call('rbind', lapply(1:ncol(weights), FUN=function(w){
          ji <- Reduce('+', mapply(extract_jacob_i, weights[,w], FUN = function(i, j) {
            i * j
          }))
          est <- as.vector(extract_ei %*% weights[,w])
          se <- sqrt(rowSums( (ji %*% vcov_d) * ji))
          out <- data.frame(est = est, se = se, obs = 1:length(est), variable = w)
          return(out)
        }))
        out_individual$variable <- colnames(weights)[out_individual$variable]
        out_individual$response <- d
        return(out_individual)
      })
      out_individual <- do.call('rbind', out_individual)
    }
  } else {
    out_individual <- NULL
  }

  return(list(
    aggregate = out_aggregate,
    individual = out_individual,
    jacobian = jacobian_net
  ))
}

predict_extended <- function(object, X, individual){
  
  stopifnot(all(rowMeans(is.na(X)) %in% c(0,1)))
  coef_object <- coef(object)
  lpi <- attr(X, 'lpi')
  if (isTRUE(attr(lpi, 'overlap'))){
    stop('Not set up for "lpi" with overlap.')
  }
  family_object <- object$family
  if (family_object$family == 'multinom' & is.null(lpi)){
    if (family_object$nlp == 1){
      lpi <- list(1:ncol(X))
    }
  }
  # Set up a flag for complex extended families (e.g., multinomial)
  complex_extended <- FALSE
  if (!is.null(family_object$predict)){
    # [1] Does it have a predict function, if so, then use that
    if (is.null(lpi)){
      
      # Use modifications of functions from "mgcv" to
      # calculate the jacobian needed for the delta method
      lp_i <- as.vector(X %*% coef_object)
      if (grepl(family_object$family, pattern='^Ordered Categorical\\(')){
        pred_obj <- ocat_jacobian(lp_i = as.vector(lp_i), 
          family_object = family_object)
      }else if (grepl(family_object$family, pattern='^Zero inflated Poisson\\(')){
        pred_obj <- zip_jacobian(lp_i = as.vector(lp_i), 
         family_object = family_object)      
      }else{
        stop('This extended family not set up for calculate_effect.')
      }
      
      pred_obj <- lapply(pred_obj, FUN=function(i){
        if (!is.matrix(i)){
          i <- matrix(i, ncol = 1)
        }
        return(i)
      })

      # # Incorrect as doesn't account for possibility of pos/neg signs      
      # old_jacob_i <- exp(sweep(log(matrix(pred_jacobian$se.fit)), MARGIN = 1,
      #                      STATS = log(abs(lp_i)), FUN = '-'
      # ))
      
      e_i <- pred_obj$fit
      jacob_i <- pred_obj$jacobian
      nlp <- ncol(e_i)
    }else{
      
      complex_extended <- TRUE
      
      lp_i <- sapply(lpi, FUN=function(lpi_d){
        as.vector(X[, lpi_d, drop = F] %*% coef_object[lpi_d])
      })
      attr(lp_i, 'lpi') <- as.list(1:length(lpi))
      
      if (grepl(family_object$family, pattern='^multinom$')){
        
        pred_obj <- multinom_jacobian(lp_i = lp_i, 
          family_object = family_object)
        
      }else{
        stop('This extended family not set up for calculate_effect.')
      }
      e_i <- pred_obj$fit
      jacob_i <- pred_obj$jacobian
      nlp <- ncol(e_i)
      
    }
  }else{
    # [2] If not, then use the linkinv for all of the relevant "link" components
    if ('linfo' %in% names(family_object)){
      list.mu.eta <- lapply(family_object$linfo, FUN=function(i){i$mu.eta})
      list.linkinv <- lapply(family_object$linfo, FUN=function(i){i$linkinv})
      nlp <- length(list.mu.eta)
      stopifnot(length(list.linkinv) == nlp)
    }else{
      list.mu.eta <- list(family_object$mu.eta)
      list.linkinv <- list(family_object$linkinv)
      nlp <- 1
      lpi <- list(1:ncol(X))
    }
    
    jacob_i <- sapply(1:nlp, FUN=function(d){
      lpi_d <- lpi[[d]]
      list.mu.eta[[d]](as.vector(X[, lpi_d] %*% coef_object[lpi_d]))
    })
    e_i <- sapply(1:nlp, FUN=function(d){
      lpi_d <- lpi[[d]]
      list.linkinv[[d]](as.vector(X[, lpi_d] %*% coef_object[lpi_d]))
    })
  }
  
  ex <- colMeans(e_i, na.rm=T)
  nrow_valid <- sum(!is.na(e_i[,1]))
  
  coef_names <- names(coef(object))
  if (individual){
    
    if (!is.list(jacob_i)){
      jacob_i <- lapply(1:nlp, FUN=function(d){
        if (is.null(lpi)){
          lpi_d <- 1:ncol(X)
        }else{
          lpi_d <- lpi[[d]]
        }
        return(Diagonal(x = jacob_i[,d]) %*% X[, lpi_d])
      })
      jacob <- lapply(jacob_i, colMeans, na.rm=T)
    }else{
      jacob_i <- lapply(jacob_i, FUN=function(ji){
        out <- do.call('cbind', sapply(1:ncol(ji), FUN=function(d){
          lpi_d <- lpi[[d]]
          return(Diagonal(x = ji[,d]) %*% X[, lpi_d])
        }))
        # out <- out[, unlist(lpi)]
        if (!isTRUE(identical(colnames(out), coef_names))){
          stop("Name misalignment when creating jacobian")
        }
        return(out)
      })
      jacob <- lapply(jacob_i, colMeans, na.rm=T)
    }
    
  }else{
    
    if (!is.list(jacob_i)){
      jacob <- lapply(1:nlp, FUN=function(d){
        if (is.null(lpi)){
          lpi_d <- 1:ncol(X)
        }else{
          lpi_d <- lpi[[d]]
        }
        ji <- colMeans(Diagonal(x = jacob_i[,d]) %*% X[, lpi_d, drop = F], na.rm=T)
        return(ji)
      })
    }else{
      jacob <- lapply(jacob_i, FUN=function(ji){
        out <- do.call('cbind', sapply(1:ncol(ji), FUN=function(d){
          lpi_d <- lpi[[d]]
          return(Diagonal(x = ji[,d]) %*% X[, lpi_d, drop = F])
        }))
        out <- out[, unlist(lpi), drop = F]
        if (!isTRUE(identical(colnames(out), coef_names))){
          stop("Name misalignment when creating jacobian")
        }
        out <- colMeans(out, na.rm=T)
        return(out)
      })
    }
    jacob_i <- NULL
    e_i <- NULL
  }
  
  if (is.null(lpi)){
    lpi <- lapply(1:nlp, FUN=function(i){1:ncol(X)})
  }
  return(
    list(
      complex_extended = complex_extended,
      expectation = ex,
      expectation_i = e_i,
      jacobian_i = jacob_i,
      jacobian = jacob,
      nrow_valid = nrow_valid,
      lpi = lpi
    )
  )
}