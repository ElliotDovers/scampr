#' Pruning for regular grids of bi-square basis functions in flat 2D space
#'
#' @param pres.locs a data frame or matrix of 2 columns describing the horizontal and vertical coordinate locations for presence locations on which to base pruning.
#' @param simple.bfs a data frame of class \code{bf.df} created by \code{simple_basis()} that describes the unpruned basis functions.
#' @param min.points Optional. a postive integer describing the least number of points required within each basis function for it to be retained in the pruning. Default is 1, meaning all basis functions with no points within their radius are pruned.
#'
#' @return a simple basis data frame of class 'bf.df' pruned to the presence locations provided. Consisting of columns: horizontal axis location, vertical axis location, scale (radius), res (resolution id).
#' @export
#'
#' @examples
#' #' # Base the basis function nodes on the locations of presence records and quadrature
#' dat <- dat <- gorillas
#' bfs <- simple_basis(nodes.on.long.edge = 9, data = dat)
#' pruned.bfs <- prune_simple_basis(dat[dat$pres == 1, c("x", "y")], bfs)
prune_simple_basis <- function(pres.locs, simple.bfs, min.points = 1) {
  # checks
  if (!"bf.df" %in% class(simple.bfs)) {
    stop("'simple.bfs': Basis functions need to be of class created by simple_basis()")
  }
  if (min.points < 1) {
    stop("'min.points' must be a postive integer")
  }
  tmp.dists <- fields::rdist(pres.locs, simple.bfs[ , 1:2])
  for.pruning <- NULL
  for (i in 1:nrow(simple.bfs)) {
    for.pruning[i] <- sum(tmp.dists[ , i] < simple.bfs$scale[i]) < min.points
  }
  new.bfs <- simple.bfs[!for.pruning, ]
  simple.bfs$pruned <- for.pruning
  attr(new.bfs, "unpruned") <- simple.bfs
  return(new.bfs)
}
