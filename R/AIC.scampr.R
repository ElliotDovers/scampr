#' AIC for objects of class 'scampr'
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
AIC.scampr <- function(object) {
  # Approximate marginal log-likelihood
  tmp.loglik <- logLik(object)
  # Akaike's Information Criteria
  tmp.aic <- -2*(tmp.loglik - length(unlist(object$starting.pars)))
  return(tmp.aic)
}
