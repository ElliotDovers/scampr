#' getCall function for objects of class 'scampr'
#'
#' @description get the call from a scampr model
#'
#' @param x a scampr model object
#' @param ... Additional arguments to the call, or arguments with changed values
#'
#' @return the updated call
#' @exportS3Method stats::getCall scampr
#'
#' @examples
#' # Get the flora PO data for one of the species
#' dat_po <- flora$po$sp1
#'
#' # Attach the quadrature to the PO data
#' dat_po <- rbind.data.frame(dat_po, flora$quad)
#'
#' # Fit an IPP model to the point pattern
#' m.ipp <- scampr(pres ~ MNT + D.Main, dat_po, model.type = "ipp")
#'
#' getCall(m.ipp)
getCall.scampr <- function(x, ...) {

  ## checks ####################################################################
  if (!is(x, "scampr")) {
    stop("x must be a fitted scampr model")
  }
  return(x$call)
}
