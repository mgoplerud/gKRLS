#' @export
SL.mgcv <- function(Y, X, newX, formula, family, id, obsWeights, bam = FALSE, ...) {
  if(!require('mgcv')) {stop("SL.mgcv requires the mgcv package, but it isn't available")} 
  
  if (is.character(formula)){
    formula <- as.formula(formula)
  }
  # https://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object
  getResponseFromFormula = function(formula) {
    if (attr(terms(as.formula(formula))    , which = 'response'))
      all.vars(formula)[1]
    else
      NULL
  }
  rformula <- getResponseFromFormula(formula)
  if (!is.null(rformula)){
    if (rformula %in% names(X)){
      warning(paste0('Outcome "', rformula, '" seems to be in "X". This is likely ill-advised'))
    }
  }
  if ('...Y' %in% names(X)){
    stop('SL.glmer cannot accept a column in "X" called "...Y". Please rename.')
  }
  X[['...Y']] <- Y
  formula <- update.formula(formula, '`...Y` ~ .')
  
  if (bam == FALSE){
    fit.mgcv <- mgcv::gam(formula, data = X, weights = obsWeights, family = family, ...)
  }else{
    fit.mgcv <- mgcv::bam(formula, data = X, weights = obsWeights, family = family, ...)
  }
  pred <- predict(fit.mgcv, newdata = newX, allow.new.levels = TRUE, type = 'response')
  fit <- list(object = fit.mgcv)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.glmer")
  return(out)
}

#' @export
predict.SL.mgcv <- function(object, newdata, allow_missing_levels = TRUE, ...){
  if(!require('mgcv')) {stop("SL.mgcv requires the mgcv package, but it isn't available")} 
  
  pred <- predict(object$object, newdata = newdata, allow.new.levels = TRUE, type = 'response', ...)
  return(pred)
  
}
