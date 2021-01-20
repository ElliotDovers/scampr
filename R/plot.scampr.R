#' Plot function for objects of class 'scampr'
#'
#' @param x a scampr model object
#' @param add.points logical indicating whether to add the presence point locations to plots
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' # Get the Eucalypt data
#' dat_po <- eucalypt[["po"]]
#' dat_pa <- eucalypt[["pa"]]
#'
#' # Fit an IPP model to the point pattern
#' m.ipp <- ippm(pres ~ elev.std, data = dat)
#'
#' # Fit a combined data model
#' m.popa <- popa(pres ~ TMP_MIN + D_MAIN_RDS, Y ~ TMP_MIN, po.data = dat_po, pa.data = dat_pa, model.type = "ipp")
#'
#' # Fit presence/absence model
#' m.pa <- pa(Y ~ TMP_MIN, pa.data = dat_pa, model.type = "ipp")
#'
#' plot(m.ipp)
#' plot(m.popa)
#' plot(m.pa)
plot.scampr <- function(x, add.points = F) {
  if (x$data.model.type == "pa") {
    plot(1 - exp(-exp(x$fitted.values)), residuals(x, type = "raw", data.type = "pa"),
         xlab = "fitted values", ylab = "raw residuals")
    readline(prompt="Hit <Return> to see next plot:")
    image(x)
  } else {
    image(x, "residuals", residual.type = "pearson")
    if (add.points) {
      tmp <- scampr:::get.data(x)
      points(tmp$pres[ , x$coord.names])
    }
    if (x$data.model.type == "popa") {
      readline(prompt="Hit <Return> to see next plot:")
      resids <- residuals(x, type = "raw", data.type = "pa")
      fits <- 1 - exp(-exp(attr(x$fitted.values, "abundance")))
      plot(fits, resids, xlab = "fitted values", ylab = "raw residuals",
           main = "Presence/Absence Residuals")
    }
    readline(prompt="Hit <Return> to see next plot:")
    image(x, "fitted")
    if (add.points) {
      points(tmp$pres[ , x$coord.names])
    }
  }
}
