#' Estimates the effective correlation range of fields from a scampr model using the estimated covariance function via \code{spatstat.core::spatcov}
#'
#' @param object a scampr model object
#' @param domain.data a data frame describing describing the domain via point locations on a regular grid, as well as any predictors used in the provided model.
#' @param cutoff.correlation a numeric describing the correlation at which the effective range is reached. Default is 0.1
#' @param field a character string, one of 'latent', 'mean' or 'bias'. Describes the field (component) for which the range is to be estimated. 'mean' refers to the fitted intensity.
#' @param plotting a logical indicating whether or not to plot the field and estimated spatial covariance function.
#'
#' @return a numeric describing the estimated range (with a numeric attribute describing the correlation at this range)
#' @export
#'
#' @examples
#' #' # Get the gorilla nesting data
#' dat <- gorillas
#'
#' # Standardise the elevation covariate
#' dat$elev.std <- scale(dat$elevation)
#'
#' # Fit a scampr model to the point pattern
#' m <- scampr(pres ~ elev.std, data = dat, model.type = "PO")
#'
#' # Use the gorillas quadrature as the domain grid
#' domain <- gorillas[gorillas$pres == 0, ]
#' estimate.correlation.range(m, domain, field = "mean", plotting = T)
estimate.correlation.range <- function(object, domain.data, cutoff.correlation = 0.1, field = c("latent", "mean", "bias"), plotting = FALSE) {
  field <- match.arg(field)

  # use the model data if domain.data is missing
  if (missing(domain.data)) {
    warning("'domain.data' is missing. Using model data however, estimated fields will be poor if this is not a representative cover of the domain.")
    domain.data <- object$data
  }
  # predict the field over the domain
  preds <- predict(object, newdata = domain.data, include.bias.accounting = TRUE)
  # set up image object
  fld.im <- switch(field,
                   latent = vec2im(attr(preds, "Zmu"), domain.data[,object$coord.names[1]], domain.data[,object$coord.names[2]]),
                   mean = vec2im(preds, domain.data[,object$coord.names[1]], domain.data[,object$coord.names[2]]),
                   bias = stop("Functionality for bias fields are not yet implemented")
  )
  cov.fn <- spatstat.core::spatcov(fld.im, correlation = TRUE) # exponentiation makes no difference
  first.drop.to.cutoff <- min(which(cov.fn$est <= cutoff.correlation))
  if (is.infinite(first.drop.to.cutoff)) {
    warning("Spatial Cov. fn. does not drop to 'cutoff.correlation' provided. First minimum correlation across reliable range provided (alternatively, try a higher cutoff).")
    first.drop.to.cutoff <- which.min(cov.fn$est[cov.fn$r >= attr(cov.fn, "alim")[1] & cov.fn$r <= attr(cov.fn, "alim")[2]])
  }
  if (plotting) {
    par(mar = c(0, 1, 1, 1))
    layout(mat = matrix(c(1, 2), nrow = 1, ncol = 2, byrow = TRUE), heights = 1, widths = c(0.5,0.5))
    plot(fld.im, main = paste0("Fitted ", field, " field"), ribbon = FALSE)
    par(mar = c(2, 1, 1, 1))
    plot(cov.fn, main = "Estimted Spatial Cov. fn.")
    abline(v = cov.fn$r[first.drop.to.cutoff], col = "blue")
    abline(h = cutoff.correlation, col = "red", lty = "dashed")
    par(mar = c(5.1, 4.1, 4.1, 2.1))
    layout(matrix(1, nrow = 1, ncol = 1))
  }
  ret.value <- cov.fn$r[first.drop.to.cutoff]
  attr(ret.value, "cor") <- cov.fn$est[first.drop.to.cutoff]
  return(ret.value)
}
