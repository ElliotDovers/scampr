#' Plot basis functions for objects of class 'bf.df' ### WIP ###
#'
#' @description Plots the bi-square basis functions created by \code{scampr::simple_basis()}.
#'
#' @param x a basis function data frame ("bf.df") object
#' @param ... additional plotting arguments
#' @param add logical indicating whether to add the basis functions to the current plot
#'
#' @return NULL
#' @export
#'
#' @importFrom graphics symbols
#'
#' @examples
#' # Create the basis functions on the flora quadrature
#' bfs <-  simple_basis(9, flora$quad)
#'
#' # plot the quadrature points
#' plot(flora$quad[,c("x", "y")], asp = 1, pch = ".", col = "blue")
#'
#' # add the plot of the bi-square basis functions
#' plot_basis(bfs, add = T)
plot_basis <- function(x, ..., add = FALSE) {
  if (!is(x, "bf.df")) {
    stop("Please provide an object of class 'bf.df'")
  }
  for (res in unique(x$res)) {
    if (which(res == unique(x$res)) > 1) {
      graphics::symbols(x[x$res == res,1], x[x$res == res,2],
              circles = x$scale[x$res == res], inches = FALSE,
              add = TRUE, ...)
    } else {
      graphics::symbols(x[x$res == res,1], x[x$res == res,2],
              circles = x$scale[x$res == res], inches = FALSE,
              add = add, xlab = colnames(x)[1], ylab = colnames(x)[2], ...)
    }
  }
}
