predict.scampr <- function(object, newdata, type = c("link", "response"), dens = c("posterior", "prior")) {
  type <- match.arg(type)
  match.arg(dens)
  if (dens == "posterior") {
    if (missing(newdata)) {
      pred <- switch(link = object$fitted.values,
                     response = exp(object$fitted.values))
    } else {
      X <- get.desgin.matrix(object$formula, newdata)
      Z <-
      pred <- switch(link = X %*% object$fixed.effects[ , 1L])
    }
  } else {

  }

}
