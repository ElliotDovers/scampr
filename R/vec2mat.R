#' Internal scampr function that converts a vector to matrix given spatial locations. For use in images
#'
#' @param vec A vector of field values of equal length to locations x.loc and y.loc
#' @param x.loc A vector of horizontal locations in 2D
#' @param y.loc A vector of vertical locations in 2D
#'
#' @return a matrix representing the locations of a regular spatial domain
#' @noRd
#'
#' @examples
#' # Get the quadrature from the gorillas data
#' quad <- gorillas[gorillas$pres == 0, ]
#'
#' elevation <- scampr:::vec2mat(vec = quad$elevation, x.loc = quad$x, y.loc = quad$y)
#' image(elevation, asp = 1)
vec2mat <- function(vec, x.loc, y.loc){
  ux <- sort(unique(x.loc))
  uy <- sort(unique(y.loc))
  nx <- length(ux)
  ny <- length(uy)
  row.ref <- match(x.loc, ux)
  col.ref <- match(y.loc, uy)
  Vec <- rep(NA, max(row.ref)*max(col.ref))
  vec.ref <- (col.ref - 1)*max(row.ref) + row.ref
  Vec[vec.ref] <- vec
  grid.mask <- matrix(Vec, max(row.ref), max(col.ref), dimnames = list(ux, uy))

  return(grid.mask)
}
