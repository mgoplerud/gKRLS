#' Machine Learning with gKRLS
#' 
#' This provides a number of functions to integrate machine learning with gKRLS
#' (and more generally `mgcv`). Integrations into `SuperLearner` and `DoubleML`
#' (via `mlr3`) are provided below.
#' 
#' SuperLearner integration is provided by `SL.mgcv` and the corresponding
#' predict method. `mgcv::bam` can be enabled by using `bam = TRUE`. Note that a
#' formula must be explicitly provided as a *character*, e.g. " ~ X1 + X2".
#' 
#' `DoubleML` integration is provided in two ways. First, one could load
#' `mlr3extralearners` to access `regr.gam` and `classif.gam`. Second, this
#' package provides `mgcv::bam` integration directly with a slight adaption of
#' the `mlr3extralearner` implementation. These can be either manually added to
#' the list of `mlr3` learners by calling `add_bam_to_mlr3()` or direct usage.
#' Examples are provided below. For `classif.bam` and `regr.bam`, the formula
#' argument is mandatory.
#' 
#' @rdname ml_gKRLS
#' @importFrom stats as.formula terms update.formula
#' @param Y Specify the outcome variable.
#' @param X All independent variables include variables in and outside the kernel. 
#' @param newX A new dataset uses for prediction. If no data provided, the original 
#' data will be used.
#' @param formula A gKRLS style formula. See details in the help(gKRLS) and help(gam).
#' @param family A string variable indicate the distribution and link function to use. 
#' The default is gaussian distribution. 
#' @param obsWeights The weights for observations.
#' @param bam A logical variable indicates whether a gKRLS model is applying to a 
#' very large dataset. The default is False.
#' @param ... Additional arguments to gam/bam.
#' 
#' @examples 
#' library(gKRLS)
#' library(DoubleML)
#' library(glue)
#' library(tidyverse)
#' 
#' n <- 5000
#' treatment <- sample(c(0,1), n, replace = T)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' state <- sample(letters, n, replace = T)
#' y = 0.3*treatment + 0.4*x2 +0.5*x3 + rnorm(n)
#' data <- data.frame(y,treatment, x2, x3, state)
#' data$state <- as.character(data$state)
#' 
#' # Double machine learning intergration
#' 
#' data_dml <- cbind(data,  model.matrix(~ state, data)[, -1]) %>%
#'  select(-state) %>% as.data.frame()
#'  
#' gkrls_formul <- as.formula(glue("~ {{paste(grep(colnames(data_dml), pattern='state', value = T), collapse = '+')}} + 
#'                                 s( x2,x3, bs = 'gKRLS')", .open = '{{', .close = '}}'))
#'                                 
#' ml_g <- LearnerRegrBam$new()
#' ml_g$param_set$values$formula <- update.formula(gkrls_formul, as.character(glue('y ~ .')))
#' ml_m <- LearnerRegrBam$new()
#' ml_m$param_set$values$formula <- update.formula(gkrls_formul, as.character(glue('treatment ~ .')))
#' 
#' ml_g$param_set$values$method <- 'REML'
#' ml_m$param_set$values$method <- 'REML'
#' 
#' data_DML <- double_ml_data_from_data_frame(
#'  df = data_dml,
#'  y_col = "y",
#'  d_cols = 'treatment',
#'  x_cols = setdiff(names(data_dml), c("y", 'treatment'))
#' ) 
#' 
#' ## Fit Partial Linear Regression
#' dml_plr <- DoubleMLPLR$new(data_DML, ml_g, ml_m)
#' dml_plr$fit()
#' dml_plr
#' 
#' ## Fit ATE
#' ### change treatment formula to classification class.
#' 
#' ml_m1 <- LearnerClassifBam$new()
#' ml_m1$param_set$values$formula <- update.formula(gkrls_formul, as.character(glue('treatment ~ .')))
#' 
#' dml_ATE <- DoubleMLIRM$new(data_DML, ml_g, ml_m1, 
#'                           score = 'ATE')
#' dml_ATE$fit()
#' dml_ATE
#' 
#' # Super learner intergration
#' 
#' library(SuperLearner) # need to pre-install and load this package 
#' 
#' data$state <- as.factor(data$state)
#' 
#' gkrls_sl <- function(...){SL.mgcv(..., bam = TRUE, formula = "~ state + s(treatment, 
#'                                                          x2, x3, bs = 'gKRLS')")}
#'                                                          
#' fit_sl <- SuperLearner(Y = data$y, 
#'                       X = data, family = 'gaussian',
#'                       SL.library = 'gkrls_sl',
#'                       verbose = T)
#'                       
#' pred <- predict(fit_sl, newdata = data)
#' @export
SL.mgcv <- function(Y, X, newX, formula, family, obsWeights, bam = FALSE, ...) {
  if(!requireNamespace('mgcv', quietly = TRUE)) {stop("SL.mgcv requires the mgcv package, but it isn't available")} 
  
  if (is.character(formula)){
    formula <- as.formula(formula)
  }
  # Modified to deal with "y ~ ."
  # https://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object
  getResponseFromFormula = function(formula, X) {
    if (attr(terms(as.formula(formula), data = X)    , which = 'response'))
      all.vars(formula)[1]
    else
      NULL
  }
  rformula <- getResponseFromFormula(formula, X)
  if (!is.null(rformula)){
    if (rformula %in% names(X)){
      warning(paste0('Outcome "', rformula, '" seems to be in "X". This is likely ill-advised'))
    }
  }
  if ('...Y' %in% names(X)){
    stop('SL.mgcv cannot accept a column in "X" called "...Y". Please rename.')
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
  class(out$fit) <- c("SL.mgcv")
  return(out)
}

#' @rdname ml_gKRLS
#' @param object A gKRLS model. 
#' @param newdata A new dataset uses for prediction. If no data provided, the original 
#' data will be used.
#' @param allow_missing_levels A logical variable indicates whether missing levels are 
#' allowed for prediction. The default is True. 
#' @export
predict.SL.mgcv <- function(object, newdata, allow_missing_levels = TRUE, ...){
  if(!requireNamespace('mgcv', quietly = TRUE)) {stop("SL.mgcv requires the mgcv package, but it isn't available")} 
  
  pred <- predict(object$object, newdata = newdata, allow.new.levels = TRUE, type = 'response', ...)
  return(pred)
  
}

#' @rdname ml_gKRLS
#' @export
add_bam_to_mlr3 <- function(){
  message('Adding "classif.bam" and "regr.bam" to DictionaryLearner.')
  mlr3::mlr_learners$add('classif.bam', LearnerClassifBam)
  mlr3::mlr_learners$add('regr.bam', LearnerClassifBam)
}

#' @importFrom R6 R6Class
test_func <- function(x){x}

#' mlr3 integrations with `bam`
#' 
#' This contains the integration of `mgcv::bam` into `mlr3` without requiring
#' explicit loading of `mlr3extralearners.`
#' @rdname mlr3_gKRLS
#' @importFrom mlr3 LearnerRegr
#' @importFrom mlr3misc invoke
#' @importFrom R6 R6Class
#' @export
LearnerRegrBam <- R6Class("LearnerRegrBam", inherit = LearnerRegr,
                         
 public = list(
   #' @description
   #' Creates a new instance of this [R6][R6::R6Class] class.
   initialize = function() {
      
     if (!requireNamespace('paradox', quietly = TRUE)){
       stop('paradox must be installed.')
     }     
     ps = paradox::ps(
       family = paradox::p_fct(default = "gaussian", levels = c("gaussian", "poisson"), tags = "train"),
       formula = paradox::p_uty(tags = "train"),
       offset = paradox::p_uty(default = NULL, tags = "train"),
       method = paradox::p_fct(
         levels = c("fREML", "GCV.Cp", "GACV.Cp", "REML", "P-REML", "ML", "P-ML"),
         default = "fREML",
         tags = "train"
       ),
       scale = paradox::p_dbl(default = 0, tags = "train"),
       select = paradox::p_lgl(default = FALSE, tags = "train"),
       knots = paradox::p_uty(default = NULL, tags = "train"),
       sp = paradox::p_uty(default = NULL, tags = "train"),
       min.sp = paradox::p_uty(default = NULL, tags = "train"),
       gamma = paradox::p_dbl(default = 1, lower = 1, tags = "train"),
       paraPen = paradox::p_uty(default = NULL, tags = "train"),
       G = paradox::p_uty(default = NULL, tags = "train"),
       drop.unused.levels = paradox::p_lgl(default = TRUE, tags = "train"),
       drop.intercept = paradox::p_lgl(default = FALSE, tags = "train"),
       discrete = paradox::p_lgl(default = FALSE, tags = "train"),
       nthreads = paradox::p_int(default = 1L, lower = 1L, tags = c("train", "threads")),
       
       # gam.control arguments
       irls.reg = paradox::p_dbl(default = 0.0, lower = 0, tags = "train"),
       epsilon = paradox::p_dbl(default = 1e-07, lower = 0, tags = "train"),
       maxit = paradox::p_int(default = 200L, tags = "train"),
       trace = paradox::p_lgl(default = FALSE, tags = "train"),
       mgcv.tol = paradox::p_dbl(default = 1e-07, lower = 0, tags = "train"),
       mgcv.half = paradox::p_int(default = 15L, lower = 0L, tags = "train"),
       rank.tol = paradox::p_dbl(default = .Machine$double.eps^0.5, lower = 0, tags = "train"),
       nlm = paradox::p_uty(default = list(), tags = "train"),
       optim = paradox::p_uty(default = list(), tags = "train"),
       newton = paradox::p_uty(default = list(), tags = "train"),
       outerPIsteps = paradox::p_int(default = 0L, lower = 0L, tags = "train"),
       idLinksBases = paradox::p_lgl(default = TRUE, tags = "train"),
       scalePenalty = paradox::p_lgl(default = TRUE, tags = "train"),
       efs.lspmax = paradox::p_int(default = 15L, lower = 0L, tags = "train"),
       efs.tol = paradox::p_dbl(default = .1, lower = 0, tags = "train"),
       scale.est = paradox::p_fct(levels = c("fletcher", "pearson", "deviance"),
                         default = "fletcher", tags = "train"),
       edge.correct = paradox::p_lgl(default = FALSE, tags = "train"),
       # Predict arguments
       block.size = paradox::p_int(default = 50000L, tags = "predict"),
       unconditional = paradox::p_lgl(default = FALSE, tags = "predict")
     )
     
     super$initialize(
       id = "regr.bam",
       packages = "mgcv",
       feature_types = c("logical", "integer", "numeric", "factor", "ordered"),
       predict_types = c("response", "se"),
       param_set = ps,
       properties = c("missings", "weights"),
       man = "regr.bam"
     )
   }
 ),
 
 private = list(
   
   .train = function(task) {
     pars = self$param_set$get_values(tags = "train")
     control_pars = pars[names(pars) %in% formalArgs(mgcv::gam.control)]
     pars = pars[!(names(pars) %in% formalArgs(mgcv::gam.control))]
     
     # set column names to ensure consistency in fit and predict
     self$state$feature_names = task$feature_names
     
     data = task$data(cols = c(task$feature_names, task$target_names))
     if ("weights" %in% task$properties) {
       pars = insert_named(pars, list(weights = task$weights$weight))
     }
     
     if (is.null(pars$formula)) {
       stop('must provide formula to "bam"')
     }else{
       if (length(task$target_names) != 1){
         stop('regr.bam only accepts single target')
       }
       pars$formula <- stats::as.formula(pars$formula)
       # Modified to deal with "y ~ ."
       # https://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object
       getResponseFromFormula = function(formula, X) {
         if (attr(terms(as.formula(formula), data = X)    , which = 'response'))
           all.vars(formula)[1]
         else
           NULL
       }
       routcome <- getResponseFromFormula(pars$formula, data)
       if (!is.null(routcome)){
         if (routcome != task$target_names){
           stop('outcome from task does not align with outcome in formula. Remove outcome from formula.')
         }          
       }
       pars$formula <- update.formula(pars$formula, paste0(task$target_names, ' ~ . '))
       
     }
     
     if (length(control_pars)) {
       control_obj = mlr3misc::invoke(mgcv::gam.control, .args = control_pars)
       pars = pars[!names(pars) %in% names(control_pars)]
     } else {
       control_obj = mgcv::gam.control()
     }
     
     mlr3misc::invoke(mgcv::bam, data = data,
                      .args = pars, control = control_obj
     )
   },
   
   .predict = function(task) {
     # get parameters with tag "predict"
     
     pars = self$param_set$get_values(tags = "predict")
     
     # get newdata and ensure same ordering in train and predict
     newdata = task$data(cols = self$state$feature_names)
     
     include_se = (self$predict_type == "se")
     
     preds = mlr3misc::invoke(
       predict,
       self$model,
       newdata = newdata,
       type = "response",
       newdata.guaranteed = TRUE,
       se.fit = include_se,
       .args = pars
     )
     
     if (include_se) {
       list(response = preds$fit, se = preds$se)
     } else {
       list(response = preds)
     }
   }
 )
)
#' 
#' @rdname mlr3_gKRLS
#' @importFrom mlr3 LearnerClassif
#' @export
LearnerClassifBam = R6Class("LearnerClassifBam",
  inherit = LearnerClassif,
  
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      
      if (!requireNamespace('paradox', quietly = TRUE)){
        stop('paradox must be installed.')
      }     
      
      ps = paradox::ps(
        family = paradox::p_fct(default = "binomial", levels = c("binomial", "multinom", "ocat"), tags = "train"),
        formula = paradox::p_uty(tags = "train"),
        offset = paradox::p_uty(default = NULL, tags = "train"),
        method = paradox::p_fct(
          levels = c("fREML", "GCV.Cp", "GACV.Cp", "REML", "P-REML", "ML", "P-ML"),
          default = "fREML",
          tags = "train"
        ),
        scale = paradox::p_dbl(default = 0, tags = "train"),
        select = paradox::p_lgl(default = FALSE, tags = "train"),
        knots = paradox::p_uty(default = NULL, tags = "train"),
        sp = paradox::p_uty(default = NULL, tags = "train"),
        min.sp = paradox::p_uty(default = NULL, tags = "train"),
        gamma = paradox::p_dbl(default = 1, lower = 1, tags = "train"),
        paraPen = paradox::p_uty(default = NULL, tags = "train"),
        G = paradox::p_uty(default = NULL, tags = "train"),
        drop.unused.levels = paradox::p_lgl(default = TRUE, tags = "train"),
        drop.intercept = paradox::p_lgl(default = FALSE, tags = "train"),
        discrete = paradox::p_lgl(default = FALSE, tags = "train"),
        nthreads = paradox::p_int(default = 1L, lower = 1L, tags = c("train", "threads")),
        
        # gam.control arguments
        irls.reg = paradox::p_dbl(default = 0.0, lower = 0, tags = "train"),
        epsilon = paradox::p_dbl(default = 1e-07, lower = 0, tags = "train"),
        maxit = paradox::p_int(default = 200L, tags = "train"),
        trace = paradox::p_lgl(default = FALSE, tags = "train"),
        mgcv.tol = paradox::p_dbl(default = 1e-07, lower = 0, tags = "train"),
        mgcv.half = paradox::p_int(default = 15L, lower = 0L, tags = "train"),
        rank.tol = paradox::p_dbl(default = .Machine$double.eps^0.5, lower = 0, tags = "train"),
        nlm = paradox::p_uty(default = list(), tags = "train"),
        optim = paradox::p_uty(default = list(), tags = "train"),
        newton = paradox::p_uty(default = list(), tags = "train"),
        outerPIsteps = paradox::p_int(default = 0L, lower = 0L, tags = "train"),
        idLinksBases = paradox::p_lgl(default = TRUE, tags = "train"),
        scalePenalty = paradox::p_lgl(default = TRUE, tags = "train"),
        efs.lspmax = paradox::p_int(default = 15L, lower = 0L, tags = "train"),
        efs.tol = paradox::p_dbl(default = .1, lower = 0, tags = "train"),
        scale.est = paradox::p_fct(levels = c("fletcher", "pearson", "deviance"),
                          default = "fletcher", tags = "train"),
        edge.correct = paradox::p_lgl(default = FALSE, tags = "train"),
        # Predict arguments
        block.size = paradox::p_int(default = 50000L, tags = "predict"),
        unconditional = paradox::p_lgl(default = FALSE, tags = "predict")
      )
      
      super$initialize(
        id = "classif.bam",
        packages = "mgcv",
        feature_types = c("logical", "integer", "numeric", "factor", "ordered"),
        predict_types = c("response", "prob"),
        param_set = ps,
        properties = c("missings", "twoclass", "weights"),
        man = "mlr_learners_classif.bam"
      )
    }
    
    
  ),
  
  private = list(
    
    .train = function(task) {
      
      pars = self$param_set$get_values(tags = "train")
      control_pars = pars[names(pars) %in% formalArgs(mgcv::gam.control)]
      pars = pars[!(names(pars) %in% formalArgs(mgcv::gam.control))]
      
      # set column names to ensure consistency in fit and predict
      self$state$feature_names = task$feature_names
      
      data = task$data(cols = c(task$feature_names, task$target_names))
      if ("weights" %in% task$properties) {
        pars = insert_named(pars, list(weights = task$weights$weight))
      }
      
      if (is.null(pars$family)){
        pars$family <- 'binomial'
      }
      
      if (is.null(pars$formula)) {
        stop('must provide formula to "bam"')
      }else{
        if (length(task$target_names) != 1){
          stop('classif.bam only accepts single target')
        }
        pars$formula <- stats::as.formula(pars$formula)
        
        # Modified to deal with "y ~ ."
        # https://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object
        getResponseFromFormula = function(formula, X) {
          if (attr(terms(as.formula(formula), data = X)    , which = 'response'))
            all.vars(formula)[1]
          else
            NULL
        }
        routcome <- getResponseFromFormula(pars$formula, data)
        if (!is.null(routcome)){
          if (routcome != task$target_names){
            stop('outcome from task does not align with outcome in formula. Remove outcome from formula.')
          }          
        }
        pars$formula <- update.formula(pars$formula, paste0(task$target_names, ' ~ . '))
      }
      
      if (length(control_pars)) {
        control_obj = mlr3misc::invoke(mgcv::gam.control, .args = control_pars)
        pars = pars[!names(pars) %in% names(control_pars)]
      } else {
        control_obj = mgcv::gam.control()
      }
      
      mlr3misc::invoke(mgcv::bam, data = data,
                       .args = pars, control = control_obj
      )
      
    },
    
    .predict = function(task) {
      
      pars = self$param_set$get_values(tags = "predict")
      lvls = task$class_names
      newdata = task$data(cols = self$state$feature_names)
      
      prob = mlr3misc::invoke(
        predict,
        self$model,
        newdata = newdata,
        type = "response",
        newdata.guaranteed = TRUE,
        .args = pars
      )
      
      if (!("multiclass" %in% task$properties)) {
        prob = cbind(as.matrix(1 - prob), as.matrix(prob))
      }
      
      colnames(prob) = lvls
      
      if (self$predict_type == "response") {
        i = max.col(prob, ties.method = "random")
        response = factor(colnames(prob)[i], levels = lvls)
        list(response = response)
      } else if (self$predict_type == "prob") {
        list(prob = prob)
      }
    }
  )
)