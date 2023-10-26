#'  Akaike's Information Criteria for objects of class 'scampr'. Note that this currently calculates mAIC
#'
#' @param object a scampr model
#' @param ... Optionally, additional scampr model objects
#' @param k a numeric describing the penalty per parameter - defaults to k = 2 i.e. classical AIC.
#'
#' @return a numeric value with the corresponding AIC (or BIC, or ..., depending on k).
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
#' m <- scampr(pres ~ elev.std, data = dat, model.type = "PO")
#'
#' AIC(m)
AIC.scampr <- function(object, ..., k = 2) {
  extra.args <- list(...)
  extra.arg.names <- as.list(substitute(list(...)))[-1L]
  if (length(extra.args) == 0) {
    return.obj <- get.single.model.aic(object, k = k)
  } else {
    aics <- NULL
    mod.names <- NULL
    for (i in 1:length(extra.args)) {
      mod.names[i] <- as.character(extra.arg.names[i])
      if (is(extra.args[[i]], "scampr")) {
        aics[i] <- get.single.model.aic(extra.args[[i]], k = k)
      } else {
        warning(paste0("argument ", mod.names[i], " is not a scampr model. AIC for this object will appear as NA"))
        aics[i] <- NA
      }
      return.obj <- cbind.data.frame(model = c(deparse(substitute(object)), mod.names), AIC = c(get.single.model.aic(object, k = k), aics))
    }
  }
  return(return.obj)
}
