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
#' m <- ipp(pres ~ elev.std, data = dat)
#'
#' AIC(m)
AIC.scampr <- function(object, ..., k = 2) {
  aic <- -2*logLik.scampr(object) + k*length(object$coefficients)
  return(aic)
}
