#'  Akaike's Information Criteria for objects of class 'scampr'
#'
#' @param object a scampr model
#' @param k a numeric describing the penalty per parameter - defaults to k = 2 i.e. classical AIC.
#' @param ... NA
#'
#' @return a numeric value with the corresponding AIC (or BIC, or ..., depending on k).
#' @export
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
#' AIC(m)
AIC.scampr <- function(object, ..., k = 2) {
  # get an identifier for the model type
  mod.id <- NULL
  if (is.na(object$approx.type)) {
    mod.id <- "ipp"
  } else {
    mod.id <- object$approx.type
  }
  # Need to get the random effect coefficients from Laplace models
  add.coef <- switch(mod.id,
         ipp = 0,
         variational = 0,
         laplace = sum(object$basis.per.res)
         )
  aic <- -2*logLik.scampr(object) + k*length(object$coefficients) + k*add.coef
  return(aic)
}
