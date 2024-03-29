#' Extract Model Residuals for objects of class 'scampr' ### WIP ###
#'
#' @param object a scampr model object
#' @param type As in residuals.ppm: A string indicating the type of residuals or weights to be used. Current options are "raw" for the raw residuals, "inverse" for the inverse-lambda residuals, and "pearson" for the Pearson residuals.
#' @param which.data Optionally, a character string indicating which data to calculate residuals on. Only relevant to integrated data models, one of "PO" or "PA".
#' @param ... NA
#'
#' @return a numeric vector, of length of the model data, containing the extracted residuals.
#' @export
#'
#' @examples
#' # Get the flora data for one of the species
#' dat_po <- flora$po$sp1
#'
#' # obtain a sample of 10,000 quadrature points for the point process model
#' set.seed(1)
#' quad.pts <- flora$quad[sample(1:nrow(flora$quad), 10000, replace = F), ]
#' set.seed(NULL)
#'
#' # Attach the quadrature points to the presence-only data
#' dat_po <- rbind.data.frame(dat_po, quad.pts)
#'
#' # Point Process Model
#' m <- scampr(pres ~ MNT + D.Main, dat_po, include.sre = F)
#'
#' # Get the residuals
#' res <- residuals(m)
residuals.scampr <- function(object, ..., type = c("raw", "inverse", "pearson"), which.data = c("PO", "PA")) {
  type <- match.arg(type)
  which.data <- match.arg(which.data)
  if (which.data == "PO") {
    if (object$model.type == "PA") {
      stop("Residuals from model and 'which.data' are not compatible")
    }
    z <- object$pt.quad.id
    lambda <- exp(object$fitted.values)
    w <- object$data[ , object$quad.weights.name]
    resids <- switch(type,
                     raw = z - (w * lambda),
                     inverse = (z/lambda) - w,
                     pearson = (z/sqrt(lambda)) - (w * sqrt(lambda)),
    )
  } else if (which.data == "PA") {
    if (object$model.type == "PO") {
      stop("Residuals from model and 'which.data' are not compatible")
    } else if (object$model.type == "popa") {
      pres_prob <- 1 - exp(-exp(attr(object$fitted.values, "abundance")))
      resp.name <- all.vars(attr(object$formula, "PA")[[2L]])
      Y <- as.numeric(attr(object$data, "PA")[ , resp.name] > 0)
    } else {
      pres_prob <- 1 - exp(-exp(object$fitted.values))
      Y <- as.numeric(object$data[ , all.vars(object$formula[[2L]])] > 0)
    }
    resids <- switch(type,
                     raw = Y - pres_prob,
                     inverse = (Y/pres_prob),
                     pearson = (Y - pres_prob)/sqrt(pres_prob * (1 - pres_prob)),
    )
  } else {
    stop("Unknown 'which.data'. Must be one of 'po' or 'pa'")
  }
  return(resids)
}
