#' Internal scampr function that extracts data information from a fitted model
#'
#' @param object
#'
#' @return
#' @noRd
#'
#' @importFrom stats as.formula
#'
#' @examples
#' # Get the gorilla nesting data
#' dat <- gorillas
#'
#' # Standardise the elevation covariate
#' dat$elev.std <- scale(dat$elevation)
#'
#' # Fit an IPP model to the point pattern
#' m.ipp <- ippm(pres ~ elev.std, data = dat)
#'
#' # Extract data
#' data.info <- get.data(m.ipp)
#' str(data.info)
get.data<- function(object) {
  if (object$data.model.type != "pa") {
    dat <- object$data
    pres.id <- object$pt.quad.id
    pres <- dat[pres.id == 1, ]
    quad <- dat[pres.id == 0, ]
    res <- list(pres = pres, quad = quad)
    if (grepl("|&|", paste(object$formula, collapse = ","), fixed = T)) {
      po.form <- stats::as.formula(strsplit(object$formula, " |&| ", fixed = T)[[1L]][1L])
      pa.form <- stats::as.formula(strsplit(object$formula, " |&| ", fixed = T)[[1L]][2L])
      attr(res, "response") <- c(all.vars(po.form[[2L]]), all.vars(pa.form[[2]]))
      attr(res, "predictors") <- c(all.vars(po.form[[3L]]), all.vars(po.form[[3L]]))
    } else {
      attr(res, "response") <- all.vars(object$formula[[2L]])
      attr(res, "predictors") <- all.vars(object$formula[[3L]])
    }
    attr(res, "coords") <- object$coord.names
    attr(res, "quad.wts") <- object$quad.weights.name
  } else {
    res <- object$data
    attr(res, "response") <- all.vars(object$formula[[2L]])
    attr(res, "predictors") <- all.vars(object$formula[[3L]])
    attr(res, "coords") <- object$coord.names
    attr(res, "quad.wts") <- object$quad.weights.name
  }
  return(res)
}
