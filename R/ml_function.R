#' Machine Learning with gKRLS
#'
#' @description This provides a number of functions to use \code{gKRLS} (and
#'   \code{mgcv} more generally) as part of machine learning algorithms.
#'   Integration into \code{SuperLearner} and \code{DoubleML} (and \code{mlr3})
#'   is described below.
#'
#' @details 
#'   \bold{Ensembles:} \code{SuperLearner} integration is provided by
#'   \code{SL.mgcv} and the corresponding predict method. \code{mgcv::bam} can be
#'   enabled by using \code{bam = TRUE}. A formula \bold{without an outcome}
#'   must be explicitly provided.
#'
#'   Please note that it is often useful to load \code{SuperLearner} before
#'   \code{gKRLS} or \code{mgcv} to avoid functions including \code{gam} and
#'   \code{s} being masked from other packages.
#'
#'   \bold{Double Machine Learning}: \code{DoubleML} integration is provided in
#'   two ways. First, one could load \code{mlr3extralearners} to access
#'   \code{regr.gam} and \code{classif.gam}.
#'
#'   Second, this package provides \code{mgcv::bam} integration directly with a
#'   slight adaption of the \code{mlr3extralearner} implementation (see
#'   \code{?LearnerClassifBam} for more details). These can be either manually
#'   added to the list of \code{mlr3} learners by calling
#'   \code{add_bam_to_mlr3()} or direct usage. Examples are provided below. For
#'   \code{classif.bam} and \code{regr.bam}, the formula argument is mandatory.
#'
#' @name ml_gKRLS
#' @importFrom stats as.formula terms update.formula
#' @param Y This is not usually directly specified in \code{SL.mgcv}, see the
#'   examples below and documentation in \code{SuperLearner} for more details.
#' @param X This is not usually directly specified in \code{SL.mgcv}, see the
#'   examples below and documentation in \code{SuperLearner} for more details.
#' @param newX This is not usually directly specified in \code{SL.mgcv}, see the
#'   examples below and documentation in \code{SuperLearner} for more details.
#' @param formula A formula used for \code{gam} or \code{bam} from
#'   \code{mgcv}. This must be specified, see the examples.
#' @param family This is not usually directly specified in \code{SL.mgcv}, see
#'   the examples below and documentation in \code{SuperLearner} for more
#'   details.
#' @param obsWeights This is not usually directly specified in \code{SL.mgcv},
#'   see the examples below and documentation in \code{SuperLearner} for more
#'   details.
#' @param bam A logical value for whether \code{mgcv::bam} should be used
#'   instead of \code{mgcv::gam}. Default is \code{FALSE}. For large datasets,
#'   this can dramatically improve estimation time. Wood et al. (2015) and
#'   \code{mgcv} provide details on \code{bam}.
#' @param ... Additional arguments to \code{mgcv::gam} and \code{mgcv::bam}.
#'
#' @references 
#' Wood, Simon N., Yannig Goude, and Simon Shaw. 2015. "Generalized Additive
#' Models for Large Data Sets." \emph{Journal of the Royal Statistical Society:
#' Series C (Applied Statistics)} 64(1):139-155.
#' 
#' @returns All three of the returned functions are usually called for use in
#'   other functions, i.e. creating objects for use in \code{SuperLearner} or
#'   adding \code{bam} models to \code{mlr3}.
#' 
#' @examples
#' set.seed(789)
#' N <- 100
#' x1 <- rnorm(N)
#' x2 <- rbinom(N, size = 1, prob = .2)
#' y <- x1^3 - 0.5 * x2 + rnorm(N, 0, 1)
#' y <- y * 10
#' X <- cbind(x1, x2, x1 + x2 * 3)
#' X <- cbind(X, "x3" = rexp(nrow(X)))
#'
#' if (requireNamespace("SuperLearner", quietly = TRUE)) {
#' # Estimate Ensemble with SuperLearner
#'   require(SuperLearner)
#'   sl_m <- function(...) { SL.mgcv(formula = ~ x1 + x2 + x3, ...) }
#'   fit_SL <- SuperLearner::SuperLearner(
#'     Y = y, X = data.frame(X),
#'     SL.library = "sl_m"
#'   )
#'   pred <- predict(fit_SL, newdata = data.frame(X))
#' }
#' # Estimate Double/Debiased Machine Learning
#' if (requireNamespace("DoubleML", quietly = TRUE)) {
#'   require(DoubleML)
#'   # Load the models; for testing *ONLY* have multiplier of 2
#'   double_bam_1 <- LearnerRegrBam$new()
#'   double_bam_1$param_set$values$formula <- ~ s(x1, x3, bs = "gKRLS", 
#'     xt = gKRLS(sketch_multiplier = NULL, sketch_size_raw = 2))
#'   double_bam_2 <- LearnerClassifBam$new()
#'   double_bam_2$param_set$values$formula <- ~ s(x1, x3, bs = "gKRLS", 
#'     xt = gKRLS(sketch_multiplier = NULL, sketch_size_raw = 2))
#'
#'   # Create data
#'   dml_data <- DoubleMLData$new(
#'     data = data.frame(X, y),
#'     x_cols = c("x1", "x3"), y_col = "y",
#'     d_cols = "x2"
#'   )
#'   # Estimate effects treatment (works for other DoubleML methods)
#'   dml_est <- DoubleMLIRM$new(
#'     data = dml_data,
#'     n_folds = 2,
#'     ml_g = double_bam_1,
#'     ml_m = double_bam_2
#'   )$fit()
#' }
#' @export
SL.mgcv <- function(Y, X, newX, formula, family, obsWeights, bam = FALSE, ...) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("SL.mgcv requires the mgcv package, but it isn't available")
  }

  if (is.character(formula)) {
    formula <- as.formula(formula)
  }
  
  # Modified to deal with "y ~ ."
  # https://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object
  getResponseFromFormula <- function(formula, X) {
    if (attr(terms(as.formula(formula), data = X), which = "response")) {
      all.vars(formula)[1]
    } else {
      NULL
    }
  }
  rformula <- getResponseFromFormula(formula, X)
  if (!is.null(rformula)) {
    if (rformula %in% names(X)) {
      warning(paste0('Outcome "', rformula, '" seems to be in "X". This is likely ill-advised'))
    }
  }
  if ("...Y" %in% names(X)) {
    stop('SL.mgcv cannot accept a column in "X" called "...Y". Please rename.')
  }
  X[["...Y"]] <- Y
  formula <- update.formula(formula, "`...Y` ~ .")

  if ("...obsWeights" %in% names(X)) {
    stop('...obsWeights cannot be a column in "X".')
  }
  X[["...obsWeights"]] <- obsWeights

  options <- list(...)

  if (any(c("formula", "data", "weights", "family") %in% names(options))) {
    stop('... must not contain "formula", "data", "weights", or "family".')
  }
  options$formula <- formula
  options$data <- X
  options$weights <- X[["...obsWeights"]]
  options$family <- family

  if (bam == FALSE) {
    gam <- mgcv::gam
    fit.mgcv <- do.call("gam", options)
  } else {
    options$id <- NULL
    bam <- mgcv::bam
    fit.mgcv <- do.call("bam", options)
  }
  pred <- predict(fit.mgcv, newdata = newX, allow.new.levels = TRUE, type = "response")
  fit <- list(object = fit.mgcv)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.mgcv")
  return(out)
}

#' @rdname ml_gKRLS
#' @param object This is not usually directly specified in \code{SL.mgcv}, see
#'   the examples below and documentation in \code{SuperLearner} for more
#'   details.
#' @param newdata This is not usually directly specified in \code{SL.mgcv}, see
#'   the examples below and documentation in \code{SuperLearner} for more
#'   details.
#' @param allow_missing_levels A logical variable that indicates whether missing
#'   levels in factors are allowed for prediction. The default is \code{TRUE}.
#'   
#' @export
predict.SL.mgcv <- function(object, newdata, allow_missing_levels = TRUE, ...) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("SL.mgcv requires the mgcv package, but it isn't available")
  }

  pred <- predict(object$object, newdata = newdata, allow.new.levels = TRUE, type = "response", ...)
  return(pred)
}

#' @rdname ml_gKRLS
#' @export
add_bam_to_mlr3 <- function() {
  message('Adding "classif.bam" and "regr.bam" to DictionaryLearner.')
  mlr3::mlr_learners$add("classif.bam", LearnerClassifBam)
  mlr3::mlr_learners$add("regr.bam", LearnerClassifBam)
}

#' mlr3 integration
#' 
#' @description This documents \code{LearnerRegrBam} and
#'   \code{LearnerClassifBam} that allow for \code{mgcv::bam} to be used in
#'   \code{mlr3} without explicitly loading \code{mlr3extralearners}. See
#'   \link{ml_gKRLS} for examples of how to use this and \code{mlr3} for
#'   discussion of the "Learner" objects.
#'   
#' @name mlr3_gKRLS
#' @importFrom mlr3 LearnerRegr
#' @importFrom R6 R6Class
#' @keywords internal
#' @references 
#' Wood, Simon N and Goude, Yannig and Simon Shaw. 2015. "Generalized Additive
#' Models for Large Data Sets." \emph{Journal of the Royal Statistical Society:
#' Series C (Applied Statistics)} 64(1):139-155.
#' @export
LearnerRegrBam <- R6Class("LearnerRegrBam",
  inherit = LearnerRegr,

  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      if (!requireNamespace("paradox", quietly = TRUE)) {
        stop("paradox must be installed.")
      }
      if (!requireNamespace("mlr3misc", quietly = TRUE)) {
        stop("mlr3misc must be installed.")
      }
      ps <- paradox::ps(
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
        scale.est = paradox::p_fct(
          levels = c("fletcher", "pearson", "deviance"),
          default = "fletcher", tags = "train"
        ),
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
      pars <- self$param_set$get_values(tags = "train")
      control_pars <- pars[names(pars) %in% formalArgs(mgcv::gam.control)]
      pars <- pars[!(names(pars) %in% formalArgs(mgcv::gam.control))]

      # set column names to ensure consistency in fit and predict
      self$state$feature_names <- task$feature_names

      data <- task$data(cols = c(task$feature_names, task$target_names))
      if ("weights" %in% task$properties) {
        pars <- insert_named(pars, list(weights = task$weights$weight))
      }

      if (is.null(pars$formula)) {
        stop('must provide formula to "bam"')
      } else {
        if (length(task$target_names) != 1) {
          stop("regr.bam only accepts single target")
        }
        pars$formula <- stats::as.formula(pars$formula)
        # Modified to deal with "y ~ ."
        # https://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object
        getResponseFromFormula <- function(formula, X) {
          if (attr(terms(as.formula(formula), data = X), which = "response")) {
            all.vars(formula)[1]
          } else {
            NULL
          }
        }
        routcome <- getResponseFromFormula(pars$formula, data)
        if (!is.null(routcome)) {
          stop("outcome must not be provided in formula to LearnerRegrBam")
          # if (routcome != task$target_names){
          #   stop('outcome from task does not align with outcome in formula. Remove outcome from formula.')
          # }
        }
        pars$formula <- update.formula(pars$formula, paste0(task$target_names, " ~ . "))
      }

      if (length(control_pars)) {
        control_obj <- mlr3misc::invoke(mgcv::gam.control, .args = control_pars)
        pars <- pars[!names(pars) %in% names(control_pars)]
      } else {
        control_obj <- mgcv::gam.control()
      }

      mlr3misc::invoke(mgcv::bam,
        data = data,
        .args = pars, control = control_obj
      )
    },

    .predict = function(task) {
      # get parameters with tag "predict"

      pars <- self$param_set$get_values(tags = "predict")

      # get newdata and ensure same ordering in train and predict
      newdata <- task$data(cols = self$state$feature_names)

      include_se <- (self$predict_type == "se")

      preds <- mlr3misc::invoke(
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

#' @rdname mlr3_gKRLS
#' @importFrom mlr3 LearnerClassif
#' @export
LearnerClassifBam <- R6Class("LearnerClassifBam",
  inherit = LearnerClassif,

  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    initialize = function() {
      if (!requireNamespace("paradox", quietly = TRUE)) {
        stop("paradox must be installed.")
      }

      ps <- paradox::ps(
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
        scale.est = paradox::p_fct(
          levels = c("fletcher", "pearson", "deviance"),
          default = "fletcher", tags = "train"
        ),
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
      pars <- self$param_set$get_values(tags = "train")
      control_pars <- pars[names(pars) %in% formalArgs(mgcv::gam.control)]
      pars <- pars[!(names(pars) %in% formalArgs(mgcv::gam.control))]

      # set column names to ensure consistency in fit and predict
      self$state$feature_names <- task$feature_names

      data <- task$data(cols = c(task$feature_names, task$target_names))
      if ("weights" %in% task$properties) {
        pars <- insert_named(pars, list(weights = task$weights$weight))
      }

      if (is.null(pars$family)) {
        pars$family <- "binomial"
      }

      if (is.null(pars$formula)) {
        stop('must provide formula to "bam"')
      } else {
        if (length(task$target_names) != 1) {
          stop("classif.bam only accepts single target")
        }
        pars$formula <- stats::as.formula(pars$formula)

        # Modified to deal with "y ~ ."
        # https://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object
        getResponseFromFormula <- function(formula, X) {
          if (attr(terms(as.formula(formula), data = X), which = "response")) {
            all.vars(formula)[1]
          } else {
            NULL
          }
        }
        routcome <- getResponseFromFormula(pars$formula, data)
        if (!is.null(routcome)) {
          stop("outcome must not be provided in formula to LearnerClassifBam")
          # if (routcome != task$target_names){
          #   stop('outcome from task does not align with outcome in formula. Remove outcome from formula.')
          # }
        }
        pars$formula <- update.formula(pars$formula, paste0(task$target_names, " ~ . "))
      }

      if (length(control_pars)) {
        control_obj <- mlr3misc::invoke(mgcv::gam.control, .args = control_pars)
        pars <- pars[!names(pars) %in% names(control_pars)]
      } else {
        control_obj <- mgcv::gam.control()
      }

      mlr3misc::invoke(mgcv::bam,
        data = data,
        .args = pars, control = control_obj
      )
    },

    .predict = function(task) {
      pars <- self$param_set$get_values(tags = "predict")
      lvls <- task$class_names
      newdata <- task$data(cols = self$state$feature_names)

      prob <- mlr3misc::invoke(
        predict,
        self$model,
        newdata = newdata,
        type = "response",
        newdata.guaranteed = TRUE,
        .args = pars
      )

      if (!("multiclass" %in% task$properties)) {
        prob <- cbind(as.matrix(1 - prob), as.matrix(prob))
      }

      colnames(prob) <- lvls

      if (self$predict_type == "response") {
        i <- max.col(prob, ties.method = "random")
        response <- factor(colnames(prob)[i], levels = lvls)
        list(response = response)
      } else if (self$predict_type == "prob") {
        list(prob = prob)
      }
    }
  )
)
