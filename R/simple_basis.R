#' Create a single resolution, regular grid of bisquare basis function nodes in flat 2D space
#'
#' @param nodes.on.long.edge integer describing the number of basis function nodes to place along the longest edge of the domain
#' @param data a data frame containing the two coordinates described by 'coord.names'
#' @param coord.names vector of character strings describing the names of the coordinates in 'data'. Ordered horizontal axis to vertical
#'
#' @return
#' @export
#'
#' @examples
simple_basis <- function(nodes.on.long.edge, data, coord.names = c("x", "y")) {
  if (!all(coord.names %in% colnames(data))) {
    stop("at least one of 'coord.names' not found in the data provided")
  }
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
  bf.info$scale <- diff(big.axis.nodes)[1]
  bf.info$res <- 1
  class(bf.info) <- c(class(bf.info), "bf.df")
  return(bf.info)
}
