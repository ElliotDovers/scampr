#' Plot the image of a field over a model quadrature
#'
#' @param x A scampr model object
#' @param z Either a single characater string of the variable name (in the model data) or one of 'fitted', 'residuals'. Alternatively, a vector of numeric values to be plotted
#' @param residual.type an optional character string for residual type if z == 'residual'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
image.scampr <- function(x, z, residual.type, ...) {

  # checks
  if (length(z == 1)) {
    if (typeof(z) == "character") {
      if (!(z %in% c(colnames(x$data), 'fitted', 'residuals'))) {
        stop(paste0(z, " must be one of 'fitted', 'residuals' or one of the columns of the model data frame.\nAlternatively supply a vector of z values directly"))
      } else {
        if (z %in% colnames(x$data)) {
          z.name <- z
          z <- x$data[x$pt.quad.id == 0, z.name]
        } else if (z == "fitted") {
          z.name <- "Fitted log-Intensity"
          z <- x$fitted.values[x$pt.quad.id == 0]
        } else {
          if (missing(residual.type)) {
            z.name <- "Residuals (raw)"
            tmp.z <- residuals.scampr(x, type = "raw")
            z <- tmp.z[x$pt.quad.id ==0]
          } else {
            z.name <- paste0("Residuals (", residual.type, ")")
            tmp.z <- residuals.scampr(x, residual.type)
            z <- tmp.z[x$pt.quad.id ==0]
          }
        }
      }
    }
  } else {
    z.name <- quote(z)
  }
  if (length(z) != sum(x$pt.quad.id == 0)) {
    stop("z vector of field values must match the number (and order) of quadrature points in the model provided")
  }
  tmp <- get.data(x)
  quad <- tmp$quad
  xs <- sort(unique(quad[ , attr(tmp, "coords")[1]]))
  ys <- sort(unique(quad[ , attr(tmp, "coords")[2]]))
  zs <- vec2mat(z, quad[ , attr(tmp, "coords")[1]], quad[ , attr(tmp, "coords")[2]])
  # image(x = xs, y = ys, z = zs, col = topo.colors(100)) # could write my own version with legend to remove need for fields
  fields::image.plot(x = xs, y = ys, z = zs, col = topo.colors(100), asp = 1, xlab = attr(tmp, "coords")[1], ylab = attr(tmp, "coords")[2], main = z.name, bty = 'n')
}
