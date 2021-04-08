#' Marginal log-likelihood for objects of class 'scampr'
#'
#' @param object a scampr model
#' @param ... Optionally, additional scampr model objects
#'
#' @return numeric describing the approximate marginal log-likelihood of the scampr model object
#' @export
#'
#' @examples
#' # Get the gorilla nesting data
#' dat <- gorillas
#'
#' # Standardise the elevation covariate
#' dat$elev.std <- scale(dat$elevation)
#'
#' # Fit a scampr model to the point pattern
#' m <- scampr(pres ~ elev.std, data = dat, model.type = "ipp")
#'
#' logLik(m)
logLik.scampr <- function(object, ...) {
  if (class(object) != "scampr") {
    stop("provided model is not of class 'scampr'")
  }
  extra.args <- list(...)
  extra.arg.names <- as.list(substitute(list(...)))[-1L]
  if (length(extra.args) == 0) {
    return.obj <- -object$value
  } else {
    lls <- NULL
    mod.names <- NULL
    for (i in 1:length(extra.args)) {
      mod.names[i] <- as.character(extra.arg.names[i])
      if (class(extra.args[[i]]) == "scampr") {
        lls[i] <- -extra.args[[i]]$value
      } else {
        lls[i] <- NA
      }
      return.obj <- cbind.data.frame(model = c(deparse(substitute(object)), mod.names), logLik = c(-object$value, lls))
    }
  }
  return(return.obj)
}
