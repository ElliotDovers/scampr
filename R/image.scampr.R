#' Plot the image of a field over a model quadrature
#'
#' @param x A scampr model object
#' @param z A vector of length equal to the number (and in order) of the quadrature in model x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
image.scampr <- function(x, z, ...) {
  if (length(z) != sum(x$pt.quad.id == 0)) {
    stop("z vector of field values must match the number of quadrature points in the model provided")
  }
  tmp <- get.data(x)
  quad <- tmp$quad
  xs <- sort(unique(quad[ , attr(tmp, "coords")[1]]))
  ys <- sort(unique(quad[ , attr(tmp, "coords")[2]]))
  zs <- vec2mat(z, quad[ , attr(tmp, "coords")[1]], quad[ , attr(tmp, "coords")[2]])
  # image(x = xs, y = ys, z = zs, col = topo.colors(100))
  fields::image.plot(x = xs, y = ys, z = zs, col = topo.colors(100), asp = 1)
}
