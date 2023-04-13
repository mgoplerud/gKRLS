#' Estimating Robust/Clustered Standard Errors with \code{mgcv}
#' 
#' This extracts the score of the log-likelihood for each observation for a
#' model estimated using \code{gam} or \code{bam}. It is necessary to allow
#' \code{mgcv} to work correctly with functions from \code{sandwich} to product
#' cluster or robust standard errors.
#' 
#' @keywords internal
#' @param x A model estimated using \code{gam} or \code{bam}.
#' @param correct_df The default, \code{TRUE}, adjusts \code{sandwich} to use the
#'   effective degrees of freedom for \code{k} instead of the number of columns
#'   of the design. \code{FALSE} uses \code{k}.
#' @param override_check The default, \code{FALSE}, allows this function to be
#'   used inside of \code{sandwich}. If only the matrix of scores is required,
#'   set to \code{TRUE}.
#' @param ... Not used for \code{estfun.gam}.
#'
#' @import sandwich
#' @importFrom stats family residuals weights model.matrix
#' @export
estfun.gam <- function(x, correct_df = TRUE, override_check = FALSE, ...){
  
  options <- list(...)
  
  if (length(options) > 0){
    stop('... not used for estfun.gam.')
  }
  
  family_x <- family(x)
  
  if (!override_check){
    # Get the function from which estfun.gam is called
    sandwich_func <- gsub(deparse(sys.calls()[[sys.nframe()-2]]), 
                          pattern='^(meat[A-Z]+|NROW)\\(.*\\)', 
                          replacement = '\\1',
                          perl = TRUE)
    # Get the environment from which the function is called
    # use to extract options
    e <- parent.frame(n = 1) 
  }
  
  if (!inherits(family_x, what = c('general.family', 'extended.family'))){
    # Extract slightly differently from gam to ensure that it works for
    # different families and link functions. See Wood (2017, Ch. 3) for discussion
    out <- as.vector( residuals(x, 'pearson') * sqrt(weights(x, 'working')) ) * model.matrix(x)
  }else{
    stop('Robust SE from sandwich not set up for general.family or extended.family')
  }

  if (correct_df & !override_check){
    
    # Get the effective degrees of freedom
    edf <- NROW(out) - sum(x$edf)
    raw_df <- NROW(out) - NCOL(out)
    
    robust_with_adjust <- c('meatHAC')
    robust_with_type <- c('meatHC', 'meatCL')
    robust_with_neither <- c('meatPC', 'meatPL')
    invalid_robust <- c('sandwich', 'vcovBS', 'vcovOPG')
    if (sandwich_func == 'NROW'){
      adjust <- 1
    }else if (sandwich_func %in% robust_with_adjust){
      
    }else if (sandwich_func %in% robust_with_type){
      type_sandwich <- e$type
      if (sandwich_func == 'meatCL'){
        if (is.null(type_sandwich)){
          adjust <- 1
        }else if (type_sandwich %in% 'HC1'){
          # adjust <- nrow(out)/(nrow(out) - edf)
        }else{stop('Only type HC1 or HC0 permitted with mgcv.')}
      }else if (sandwich_func == 'meatHC'){
        if (!is.null(e$omega) | is.null(type_sandwich)){
          adjust <- 1
        }else if (type_sandwich %in% c("HC", "HC0")){
          adjust <- 1
        }else if (type_sandwich %in% c('const', 'HC1')){
          adjust <- raw_df/edf
        }else{stop('HC2, HC3, HC4, HC5 are not set up for mgcv.')}
      }else{stop('meatCL and meatHC are only options.')}
    }else if (sandwich_func %in% robust_with_neither){
      
    }else if (sandwich_func %in% invalid_robust){
      stop(paste0(sandwich_func, ' not yet implemented for mgcv.'))
    }else{
      stop(paste0(sandwich_func, ' not recognized as one from sandwich::;', 
        ' set override_check = TRUE to use in other settings.'))
    }
  }else{
    adjust <- 1
  }
  adjust <- sqrt(adjust)
  # Do degree of freedom correction for penalized GAMs
  out <- adjust * out
  return(out)
}