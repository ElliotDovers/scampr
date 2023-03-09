#' vcov for objects of class 'scampr'
#'
#' @param object a scampr model object
#' @param ... NA
#' @param getJointPrecision a logical indicating whether to calculate the full joint precision matrix, including random effects (TRUE) or just the fixed effects (FALSE).
#'
#' @return A matrix of the estimated covariances between the parameter estimates in the linear or non-linear predictor of the model. This should have row and column names corresponding to the parameter names given by the coef method.
#' @export
#'
#' @importFrom TMB MakeADFun sdreport
#'
#' @examples
#' # Get the gorilla nesting data
#' dat <- gorillas
#'
#' # Standardise the elevation covariate
#' dat$elev.std <- scale(dat$elevation)
#'
#' # Fit a scampr model to the point pattern
#' m <- scampr(pres ~ elev.std, data = dat, model.type = "ipp")
#'
#' vcov(m)
vcov.scampr <- function(object, ..., getJointPrecision = T) {

  # get the TMB data inputs so we can recalculate the Hessian
  inputs <- do.call("get.TMB.data.input", object$tmb.call.list)

  # adjust the start pars to be the fitted par values
  inputs$tmb.pars <- scampr:::update.starting.parameters(object, inputs$tmb.pars, target.approx.type = object$approx.type)

  # create the objective function at the newly updated parameters
  obj <- make.objective.function(inputs)

  # get the sdreport that produces the full precision matrix
  tmp <- TMB::sdreport(obj, getJointPrecision = getJointPrecision)
  if (is.null(tmp$jointPrecision)) {
    vcov_mat <- tmp$cov.fixed
  } else {
    Hess <- as.matrix(tmp$jointPrecision)
    vcov_mat <- solve(Hess)
  }
  # re-name the rows and columns
  if (object$approx.type == "laplace") {
    rownames(vcov_mat)[rownames(vcov_mat) != "random"] <- names(object$coefficients)
    colnames(vcov_mat)[colnames(vcov_mat) != "random"] <- names(object$coefficients)
  } else {
    rownames(vcov_mat) <- names(object$coefficients)
    colnames(vcov_mat) <- names(object$coefficients)
  }

  return(vcov_mat)
}
