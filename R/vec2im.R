#' Internal scampr function that converts a vector to spatstat image given spatial locations.
#'
#' @param vec A vector of field values of equal length to locations x.loc and y.loc
#' @param x.loc A vector of horizontal locations in 2D
#' @param y.loc A vector of vertical locations in 2D
#'
#' @return a spatstat image (im) object
#' @export
#'
#' @importFrom spatstat.geom im
#'
#' @examples
#' #' # Get the quadrature from the gorillas data
#' quad <- gorillas[gorillas$pres == 0, ]
#'
#' elevation <- scampr:::vec2im(vec = quad$elevation, x.loc = quad$x, y.loc = quad$y)
vec2im <- function(vec, x.loc, y.loc) {
  var.mat <- vec2mat(vec, x.loc, y.loc)
  var.im <- spatstat.geom::im(t(var.mat), xcol = sort(unique(x.loc)), yrow = sort(unique(y.loc)))
  return(var.im)
}
