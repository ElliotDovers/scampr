#'  Akaike's Information Criteria for objects of class 'scampr'
#'
#' @param object a scampr model
#'
#' @return
#' @export
#'
AIC.scampr <- function(object, k = 2) {
  aic <- -2*logLik(object) + k*length(object$coefficients)
  return(aic)
}
