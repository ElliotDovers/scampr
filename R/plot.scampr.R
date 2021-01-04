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
