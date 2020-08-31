#' Marginal log-likelihood for objects of class 'scampr'
#'
#' @param object
#'
#' @return numeric describing the approximate marginal log-likelihood of the scampr model object
#' @export
#'
#' @examples
logLik.scampr <- function(object) {
  return(-object$value)
}
