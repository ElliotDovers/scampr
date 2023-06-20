#' Partition data into k spatially blocked folds for cross validation
#'
#' @description This function takes in a data frame containing coordinates and creates a spatially blocked set of k folds.
#'
#' @param data a data frame with columns \code{coord.names} to be split into spatial k folds.
#' @param rangeX a vector of length 2 denoting the range of the horizontal coordinate. If missing this will be obtained from \code{data} - note that this assumes the region of interest is completely represented by points in the data frame provided.
#' @param rangeY a vector of length 2 denoting the range of the vertical coordinate. If missing this will be obtained from \code{data} - note that this assumes the region of interest is completely represented by points in the data frame provided.
#' @param k an integer denoting the number of folds. Default is 5.
#' @param coord.names a vector of length two describing the column names of the horizontal and vertical coordinates in \code{data} in this order. Default is "x" and "y" respectively.
#'
#' @return a vector of length \code{nrow(data)} with integers denoting the fold into which each point falls.
#' @export
#'
#' @examples
#' #' # Get the gorilla nesting data
#' dat <- gorillas
#'
#' # create spatial CV folds
#' dat$fold <- make.spatial.folds(dat)
#'
#' # view the folds
#' table(dat$fold)
make.spatial.folds <- function(data, rangeX, rangeY, k = 5, coord.names = c("x", "y")) {
  # if ranges are missing use data provided - assumes the data spans the region of interest
  if (missing(rangeX)) {
    rangeX <- range(data[ , coord.names[1]])
  }
  if (missing(rangeY)) {
    rangeY <- range(data[ , coord.names[2]])
  }
  # calculate the minimums
  min_X <- min(rangeX)
  min_Y <- min(rangeY)
  # calculate the minimum span in either coordinate
  min_span <- min(diff(rangeX), diff(rangeY))
  # set up spatial breaks according to the number of folds (plus one as we are after gaps in the sequence)
  break_pts <- seq(0, min_span, length.out = k + 1)
  # initialise the folds in coordinates
  fold.X <- NULL
  fold.Y <- NULL
  # loop through break points and assign coordinates (individually) to folds
  for (i in 1:k) {
    fold.X[((data[ , coord.names[1]] - min_X) %% min_span < break_pts[i + 1]) & ((data[ , coord.names[1]] - min_X) %% min_span >= break_pts[i])] <- i
    fold.Y[((data[ , coord.names[2]] - min_Y) %% min_span < break_pts[i + 1]) & ((data[ , coord.names[2]] - min_Y) %% min_span >= break_pts[i])] <- i
  }
  # combine individual coordinate folds for assigned folds
  fold.id <- ((fold.X + fold.Y) %% k) + 1
  return(factor(fold.id))
}
