#' Extract Model Residuals for objects of class 'scampr'
#'
#' @param object a scampr model object
#' @param type As in residuals.ppm: A string indicating the type of residuals or weights to be used. Current options are "raw" for the raw residuals, "inverse" for the inverse-lambda residuals, and "pearson" for the Pearson residuals.
#' @param data.type a character string indicating which data type to produce residuals from.
#' @param ... NA
#'
#' @return a numeric vector, of length of the model data, containing the extracted residuals.
#' @export
#'
#' @examples
#' # Get the flora data for one of the species
#' dat_po <- flora$po$sp1
#' dat_pa <- flora$pa
#'
#' # Attach the quadrature to the PO data
#' dat_po <- rbind.data.frame(dat_po, flora$quad)
#'
#' # Set up a simple 2D grid of basis functions to fit a LGCP model to the data
#' bfs <- simple_basis(nodes.on.long.edge = 9, data = dat_po)
#'
#' # Fit an IPP model to the point pattern
#' m.ipp <- scampr(pres ~ MNT + D.Main, dat_po, model.type = "ipp")
#'
#' # Fit a combined data model
#' m.comb <- scampr(pres ~ MNT + D.Main, dat_po, sp1 ~ MNT,
#' dat_pa, model.type = "ipp")
#'
#' # Fit presence/absence model
#' m.pa <- scampr(pa.formula = sp1 ~ MNT, pa.data = dat_pa, model.type = "ipp")
#'
#' \dontrun{
#' # Fit a LGCP model to the point pattern
#' m.lgcp_va <- scampr(pres ~ MNT + D.Main, dat_po,
#' model.type = "variational", simple.basis = bfs)
#' # Or
#' m.lgcp_lp <- scampr(pres ~ MNT + D.Main, dat_po,
#' approx.with = "laplace", simple.basis = bfs)
#' }
#' residuals(m.ipp)
#' residuals(m.comb)
#' residuals(m.comb, data.type = "pa")
#' residuals(m.pa, data.type = "pa")
#' \dontrun{
#' residuals(m.lgcp_va)
#' residuals(m.lgcp_lp, type = "pearson")
#' }
residuals.scampr <- function(object, ..., type = c("raw", "inverse", "pearson"), data.type = c("po", "pa")) {
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
      resp.name <- all.vars(attr(object$formula, "pa")[[2L]])
      Y <- as.numeric(attr(object$data, "pa")[ , resp.name] > 0)
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
    stop("Unknown 'data.type'. Must be one of 'po' or 'pa'")
  }
  return(resids)
}
