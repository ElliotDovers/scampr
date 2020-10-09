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
plot.scampr <- function(x, add.points = F, ...) {
  image(x, "residuals", residual.type = "pearson")
  if (add.points) {
    tmp <- scampr:::get.data(x)
    points(tmp$pres[ , x$coord.names])
  }
  readline(prompt="Hit <Return> to see next plot:")
  image(x, "fitted")
  if (add.points) {
    points(tmp$pres[ , x$coord.names])
  }
}
