#' Internal scampr function that extracts data list of presence and quadrature from a fitted model
#'
#' @param object
#'
#' @return
#'
#' @examples
get.data<- function(object) {
  dat <- object$data
  pres.id <- object$pt.quad.id
  pres <- dat[pres.id == 1, ]
  quad <- dat[pres.id == 0, ]
  res <- list(pres = pres, quad = quad)
  attr(res, "response") <- all.vars(object$formula[[2]])
  attr(res, "predictors") <- all.vars(object$formula[[3]])
  attr(res, "coords") <- object$coord.names
  attr(res, "quad.wts") <- object$quad.weights.name
  return(res)
}
