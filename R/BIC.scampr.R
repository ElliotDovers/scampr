#'  Wrapper for Akaike's Information Criteria for objects of class 'scampr' with k = log(N)
#'
#' @param object a scampr model
#' @param ... Optionally, additional scampr model objects
#'
#' @return a numeric value with the corresponding BIC
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
#' BIC(m)
BIC.scampr <- function(object, ...) {

  extra.args <- list(...)
  extra.arg.names <- as.list(substitute(list(...)))[-1L]
  if (length(extra.args) == 0) {
    if (object$model.type == "IDM") {
      N <- nrow(object$data) + nrow(attr(object$data, "PA"))
    } else {
      N <- nrow(object$data)
    }
    return.obj <- get.single.model.aic(object, k = log(N))
  } else {
    bics <- NULL
    mod.names <- NULL
    for (i in 1:length(extra.args)) {
      mod.names[i] <- as.character(extra.arg.names[i])
      if (is(extra.args[[i]], "scampr")) {
        if (extra.args[[i]]$model.type == "IDM") {
          N <- nrow(extra.args[[i]]$data) + nrow(attr(extra.args[[i]]$data, "PA"))
        } else {
          N <- nrow(extra.args[[i]]$data)
        }
        bics[i] <- get.single.model.aic(extra.args[[i]], k = log(N))
      } else {
        warning(paste0("argument ", mod.names[i], " is not a scampr model. AIC for this object will appear as NA"))
        bics[i] <- NA
      }
      if (object$model.type == "IDM") {
        N <- nrow(object$data) + nrow(attr(object$data, "PA"))
      } else {
        N <- nrow(object$data)
      }
      return.obj <- cbind.data.frame(model = c(deparse(substitute(object)), mod.names), BIC = c(get.single.model.aic(object, k = log(N)), bics))
    }
  }
  return(return.obj)
}
