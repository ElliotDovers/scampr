#' Plot diagnositcs for objects of class 'scampr'
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
#' # Get the Eucalypt data
#' dat_po <- eucalypt[["po"]]
#' dat_pa <- eucalypt[["pa"]]
#'
#' # Fit an IPP model to the point pattern
#' m.ipp <- ipp(pres ~ TMP_MIN + D_MAIN_RDS, data = dat_po)
#'
#' # Fit a combined data model
#' m.popa <- popa(pres ~ TMP_MIN + D_MAIN_RDS, Y ~ TMP_MIN,
#' po.data = dat_po, pa.data = dat_pa, model.type = "ipp")
#'
#' # Fit presence/absence model
#' m.pa <- pa(Y ~ TMP_MIN, pa.data = dat_pa, model.type = "ipp")
#'
#' plot(m.ipp)
#' plot(m.popa)
#' plot(m.pa)
plot.scampr <- function(x, ..., which = c("residuals", "fitted"), add.points = F) {
  which <- match.arg(which, several.ok = T)
  if (x$data.model.type == "pa") {
    if ("residuals" %in% which) {
      graphics::plot.default(1 - exp(-exp(x$fitted.values)), residuals.scampr(x, type = "raw", data.type = "pa"),
                             xlab = "fitted values", ylab = "raw residuals")
    }
    # readline(prompt="Hit <Return> to see next plot:")
    if ("fitted" %in% which) {
      image.scampr(x, ...)
    }
  } else {
    if ("residuals" %in% which) {
      image.scampr(x, "residuals", residual.type = "pearson", ...)
      if (add.points) {
        graphics::points(x$data[x$pt.quad.id == 1, x$coord.names])
      }
    }

    # if (x$data.model.type == "popa") {
    #   readline(prompt="Hit <Return> to see next plot:")
    #   resids <- residuals.scampr(x, type = "raw", data.type = "pa")
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
