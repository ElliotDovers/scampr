#'  Akaike's Information Criteria for objects of class 'scampr'
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
AIC.scampr <- function(object, ...) {
  aic <- -2*(logLik(object) - length(object$coefficients))
  return(aic)
}
