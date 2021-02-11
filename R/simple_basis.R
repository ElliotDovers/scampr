#' Create a single resolution, regular grid of basis function nodes in flat 2D space
#'
#' @description Uses data provided to place points across the range of the provided coordinates. The points are set out along the wider coordinate range according to the number provided. Points along the shorter ranging coordinate then adhere to this spacing. The points are then returned as a data frame of the same format as in 'FRK::auto_basis', the scale is provided based on one of two radius options. For use in bi-square basis functions with local support.
#'
#' @param nodes.on.long.edge integer describing the number of basis function nodes to place along the longest edge of the domain
#' @param data a data frame containing the two coordinates described by 'coord.names'
#' @param radius.type character string describing the type of radius length to use. One of 'diag' = diagonal dist. between nodes or 'limiting' = sqrt(Domain Area)/log(k).
#' @param coord.names vector of character strings describing the names of the coordinates in 'data'. Ordered horizontal axis to vertical
#'
#' @return a simple basis data frame of class 'bf.df'. Consisting of columns: horizontal axis location, vertical axis location, scale (radius), res (resolution id).
#' @export
#'
#' @examples
#' # Base the basis function nodes on the locations of presence records and quadrature
#' bfs1 <- simple_basis(nodes.on.long.edge = 9, data = eucalypt[["po"]])
#'
#' # Base the basis function nodes on the locations of survey sites from presence/absence data
#' bfs2 <- simple_basis(nodes.on.long.edge = 9, data = eucalypt[["pa"]])
simple_basis <- function(nodes.on.long.edge, data, radius.type = c("diag", "limiting"), coord.names = c("x", "y")) {
  if (!all(coord.names %in% colnames(data))) {
    stop("at least one of 'coord.names' not found in the data provided")
  }
  radius.type <- match.arg(radius.type)

  # Get the data ranges
  xrange <- range(data[ , coord.names[1]], na.rm = T)
  yrange <- range(data[ , coord.names[2]], na.rm = T)
  # Get the axis sizes
  dx <- diff(xrange)
  dy <- diff(yrange)
  # Get the largest axis identifier
  big.axis.id <- which.max(c(dx, dy))
  small.axis.id <- c(1, 2)[!c(1, 2) %in% big.axis.id]
  # Jitter the first node placement by 1/4 the distance between nodes (so basis fn doesn't land on an edge point)
  big.axis.jitter <- c(dx, dy)[big.axis.id] / (nodes.on.long.edge * 4)
  big.axis.nodes <- seq(c(xrange[1], yrange[1])[big.axis.id] - big.axis.jitter, c(xrange[2], yrange[2])[big.axis.id] + big.axis.jitter, length.out = nodes.on.long.edge)
  # Create small axis nodes using the same distance between nodes
  small.axis.nodes <- seq(c(xrange[1], yrange[1])[small.axis.id] - big.axis.jitter, c(xrange[2], yrange[2])[small.axis.id] + big.axis.jitter, by = diff(big.axis.nodes)[1])
  node.list <- list()
  node.list[[big.axis.id]] <- big.axis.nodes
  node.list[[small.axis.id]] <- small.axis.nodes
  bf.info <- as.data.frame(expand.grid(node.list[[1]], node.list[[2]]))
  colnames(bf.info) <- coord.names
  bf.info$scale <- switch(radius.type,
                          diag = diff(big.axis.nodes)[1] * sqrt(2),
                          limiting = sqrt(dx * dy) / log(nrow(bf.info))
  )
  bf.info$res <- 1
  class(bf.info) <- c(class(bf.info), "bf.df")
  return(bf.info)
}
