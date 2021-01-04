#' Extract Model Residuals for objects of class 'scampr'
#'
#' @param object a scampr model object
#' @param type As in residuals.ppm: A string indicating the type of residuals or weights to be used. Current options are "raw" for the raw residuals, "inverse" for the inverse-lambda residuals, and "pearson" for the Pearson residuals.
#' @param data.type a character string indicating which data type to produce residuals from.
#'
#' @return
#' @export
#'
#' @examples
residuals.scampr <- function(object, type = c("raw", "inverse", "pearson"), data.type = c("po", "pa")) {
  type <- match.arg(type)
  data.type <- match.arg(data.type)
  if (data.type == "po") {
    if (object$data.model.type == "pa") {
      stop("Residuals from model and 'data.type' are not compatible")
    }
    z <- object$pt.quad.id
    lambda <- exp(object$fitted.values)
    w <- object$data[ , object$quad.weights.name]
    resids <- switch(type,
                     raw = z - (w * lambda),
                     inverse = (z/lambda) - w,
                     pearson = (z/sqrt(lambda)) - (w * sqrt(lambda)),
    )
  } else if (data.type == "pa") {
    if (object$data.model.type == "po") {
      stop("Residuals from model and 'data.type' are not compatible")
    } else if (object$data.model.type == "popa") {
      pres_prob <- 1 - exp(-exp(attr(object$fitted.values, "abundance")))
      tmp <- scampr:::get.data(object)
      resp.name <- attr(tmp, "response")[2]
      Y <- as.numeric(attr(object$data, "pa")[ , resp.name] > 0)
    } else {
      pres_prob <- 1 - exp(-exp(object$fitted.values))
      tmp <- scampr:::get.data(object)
      Y <- as.numeric(tmp[ , attr(tmp, "response")] > 0)
    }
    resids <- switch(type,
                     raw = Y - pres_prob,
                     inverse = (Y/pres_prob),
                     pearson = (Y - pres_prob)/sqrt(pres_prob * (1 - pres_prob)),
    )
  } else {
    stop("Unknown 'data.type'. Must be one of 'po' or 'pa'")
  }
  return(resids)
}
