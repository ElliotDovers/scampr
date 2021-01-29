#' Marginal log-likelihood for objects of class 'scampr'
#'
#' @param object a scampr model
#' @param ... NA
#'
#' @return numeric describing the approximate marginal log-likelihood of the scampr model object
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
#' m <- ippm(pres ~ elev.std, data = dat)
#'
#' logLik(m)
logLik.scampr <- function(object, ...) {
  return(-object$value)
}
