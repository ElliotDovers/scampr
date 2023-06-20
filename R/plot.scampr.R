#' Plot diagnositcs for objects of class 'scampr' ### WIP ###
#'
#' @description Currently plots residuals as in 'spatstat', as well as the fitted intensity, over the quadrature used to model the point pattern. In the case of presence/absence data fitted vs. residuals is plotted.
#'
#' @param x a scampr model object
#' @param ... additional plotting arguments
#' @param which if a subset of the plots is required, specify one of 'residuals' or 'fitted': plots Pearson residuals and fitted intensity respectively for PPM (for a PA data model: plots fitted vs. raw residuals and the survey data for 'residuals' and 'fitted' resp.).
#' @param add.points logical indicating whether to add the presence point locations to plots
#'
#' @return See spatstat::plot.ppm
#' @exportS3Method graphics::plot scampr
#'
#' @importFrom graphics plot.default points
#'
#' @examples
#' # Get the flora data for one of the species
#' dat_po <- flora$po$sp1
#' dat_pa <- flora$pa
#'
#' # Attach the quadrature to the PO data
#' dat_po <- rbind.data.frame(dat_po, flora$quad)
#'
#' # Fit an IPP model to the point pattern
#' m.ipp <- scampr(pres ~ MNT + D.Main, data = dat_po, model.type = "ipp")
#'
#' # Fit a combined data model
#' m.comb <- scampr(pres ~ MNT + D.Main, dat_po, sp1 ~ MNT,
#' dat_pa, model.type = "ipp")
#'
#' # Fit presence/absence model
#' m.pa <- scampr(pa.formula = sp1 ~ MNT, pa.data = dat_pa, model.type = "ipp")
#'
#' plot(m.ipp)
#' plot(m.comb)
#' plot(m.pa)
plot.scampr <- function(x, ..., which = c("residuals", "fitted"), add.points = F) {
  xtrargs <- list(...)
  which <- match.arg(which, several.ok = T)
  if (x$model.type == "PA") {
    if ("residuals" %in% which) {
      graphics::plot.default(1 - exp(-exp(x$fitted.values)), residuals.scampr(x, type = "raw", which.data = "PA"),
                             xlab = "fitted values", ylab = "raw residuals")
    }
    # readline(prompt="Hit <Return> to see next plot:")
    if ("fitted" %in% which) {
      image.scampr(x, ...)
    }
  } else {
    if ("residuals" %in% which) {
      if (!"residual.type" %in% names(xtrargs)) {
        xtrargs$residual.type <- "pearson"
      }
      xtrargs$x <- x
      xtrargs$z <- "residuals"
      do.call(image.scampr, xtrargs)
      # image.scampr(x, "residuals", residual.type = "pearson", ...)
      if (add.points) {
        graphics::points(x$data[x$pt.quad.id == 1, x$coord.names])
      }
    }

    # if (x$model.type == "IDM") {
    #   readline(prompt="Hit <Return> to see next plot:")
    #   resids <- residuals.scampr(x, type = "raw", which.data = "PA")
    #   fits <- 1 - exp(-exp(attr(x$fitted.values, "abundance")))
    #   graphics::plot.default(fits, resids, xlab = "fitted values", ylab = "raw residuals",
    #        main = "Presence/Absence Residuals")
    # }
    # readline(prompt="Hit <Return> to see next plot:")
    if ("fitted" %in% which) {
      image.scampr(x, "fitted", ...)
      if (add.points) {
        graphics::points(x$data[x$pt.quad.id == 1, x$coord.names])
      }
    }
  }
}
