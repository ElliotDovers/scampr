#' Marginal log-likelihood for objects of class 'scampr'
#'
#' @param object a scampr model
#'
#' @return numeric describing the approximate marginal log-likelihood of the scampr model object
#' @export
#'
#' @examples
logLik.scampr <- function(object) {
  return(-object$value)
}
