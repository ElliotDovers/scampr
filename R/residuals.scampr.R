#' Extract Model Residuals for objects of class 'scampr'
#'
#' @param object
#' @param type As in residuals.ppm: A string indicating the type of residuals or weights to be used. Current options are "raw" for the raw residuals, "inverse" for the inverse-lambda residuals, and "pearson" for the Pearson residuals.
#'
#' @return
#' @export
#'
#' @examples
residuals.scampr <- function(object, type = "raw") {
  z <- object$pt.quad.id
  lambda <- exp(object$fitted.values)
  w <- object$data[ , object$quad.weights.name]
  switch(type,
         raw = z - (w * lambda),
         inverse = (z/lambda) - w,
         pearson = (z/sqrt(lambda)) - (w * sqrt(lambda)),
  )
}
