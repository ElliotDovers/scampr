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
#' dat_po <- eucalypt[["po"]]
#' dat_pa <- eucalypt[["pa"]]
#'
#' # Set up a simple 2D grid of basis functions to fit a LGCP model to the data
#' bfs <- simple_basis(nodes.on.long.edge = 9, data = dat_po)
#'
#' # Fit an IPP model to the point pattern
#' m.ipp <- ippm(pres ~ TMP_MIN, data = dat_po)
#'
#' # Fit a combined data model
#' m.popa <- popa(pres ~ TMP_MIN + D_MAIN_RDS, Y ~ TMP_MIN, po.data = dat_po, pa.data = dat_pa, model.type = "ipp")
#'
#' # Fit presence/absence model
#' m.pa <- pa(Y ~ TMP_MIN, pa.data = dat_pa, model.type = "ipp")
#'
#' # Fit a LGCP model to the point pattern
#' m.lgcp_va1 <- po(pres ~ TMP_MIN + D_MAIN_RDS, data = dat_po, model.type = "variational", simple.basis = bfs)
#' # Or
#' m.lgcp_va2 <- lgcpm(pres ~ TMP_MIN + D_MAIN_RDS, data = dat_po, approx.with = "variational", simple.basis = bfs)
#'
#' residuals(m.ipp, test_po)
#' residuals(m.popa, test_po)
#' residuals(m.popa, test_pa, data.type = "pa")
#' residuals(m.pa, data.type = "pa")
#' residuals(m.lgcp_va1)
#' residuals(m.lgcp_va2, type = "pearson")
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
