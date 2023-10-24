#' Internal scampr function that calculates the basis function evaluations at point locations provided. From either a fitted model or data.frame of class 'bf.df'
#'
#' @param object either a fitted model or data.frame of class 'bf.df' created by 'simple_basis'
#' @param point.locations either a matrix or data.frame describing the point locations in the same order as in 'object'
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to create sparse or dense matrix.
#'
#' @return a data.frame (sparse or dense depending on parameter bf.matrix.type)
#' @export
#'
#' @importFrom FRK eval_basis
#' @importFrom methods as
#' @importFrom sp spDists
#'
#' @examples
#' # Get the gorilla nesting data
#' dat <- gorillas
#'
#' # Standardise the elevation covariate
#' dat$elev.std <- scale(dat$elevation)
#'
#' # Set up a simple 2D grid of basis functions to fit a LGCP model to the data
#' bfs <- simple_basis(nodes.on.long.edge = 9, data = dat)
#'
#' # Fit a LGCP model using variational approximation
#' m.lgcp_va <- scampr(pres ~ elev.std, data = dat, model.type = "variational", simple.basis = bfs)
#'
#' # Get some new locations
#' new.sites <- dat[sample(1:nrow(dat), 10), c("x", "y")]
#'
#' # Calculate the basis function matrix at the new sites
#' scampr:::get.bf.matrix(m.lgcp_va, new.sites)
get.bf.matrix <- function(object, point.locations, bf.matrix.type = c("sparse", "dense")) {

  ## WANT TO ADD IN FUNCTIONALITY TO MAKE BF MATRIX ORTHOGONAL TO FIXED EFFECT MATRIX AS IN MGCV
  # Uses: PX = X(X′X)−1X′ and QX = I − PX

  bf.matrix.type <- match.arg(bf.matrix.type)

  # check if the provided object is a scampr model object and if so extract the basis functions
  if (is(object, "scampr")) {
    if (!is.null(object$basis.functions)) { # check that the model has basis functions
      object <- object$basis.functions
    } else {
      stop("No 'basis.functions' found in the scampr model provided to 'get.bf.matrix'")
    }
  }

  # alter approach based on whether basis functions are to be scampr's simple version or those of FRK
  if (is(object, "Basis")) {
    # evaluate the basis function matrix using FRK::eval_basis
    bf.mat <- FRK::eval_basis(basis = object, as.matrix(point.locations))
    # since these default to sparse we need to convert to dense if needed
    if (bf.matrix.type == "dense") {
      bf.mat <- methods::as(bf.mat, "matrix")
    }
    # set the basis information matrix
    bf.info <- object@df
    # add in a class for the basis information matrix
    class(bf.info) <- c(class(bf.info), "bf.df")
  } else if (is(object, "bf.df")) {
    # initialise the basis function matrix
    bf.mat <- NULL
    # loop through each resolution of basis functions
    for (res in unique(object$res)) {
      # obtain the radius
      radius <- object$scale[object$res == res][1]
      # calculate the distances from points to basis function nodes
      if (nrow(point.locations) == 0) { # allows the edge case where there are no points (could arise in cross-validation for example)
        dist.mat <- matrix(NA, nrow = 0, ncol = nrow(as.matrix(object[,1:2][object$res == res, ])))
      } else {
        dist.mat <- sp::spDists(as.matrix(point.locations), as.matrix(object[,1:2][object$res == res, ]), longlat = attr(object, "longlat"))
      }
      # if (is.null(attr(object, "longlat"))) {
      #   dist.mat <- fields::rdist(point.locations, object[,1:2][object$res == res, ])
      # } else {
      #   dist.mat <- sp::spDists(as.matrix(point.locations), as.matrix(object[,1:2][object$res == res, ]), longlat = attr(object, "longlat"))
      # }
      # We want to get these matrices into sparse format ASAP to save memory and computation times

      # record instances of any point exactly on the nodes of the basis functions
      pts.on.nodes <- methods::as(dist.mat == 0, "sparseMatrix")
      # set distances beyond the radius to zero
      dist.mat[dist.mat > radius] <- 0
      # if sparse we can save computation time here
      if (bf.matrix.type == "sparse") {
        dist.mat <- methods::as(dist.mat, "sparseMatrix")
      }
      # calculate the bi-square basis function matrix
      Z <- ((dist.mat != 0 | pts.on.nodes) - (dist.mat / radius)^2)^2
      # add resolution to matrix via columns
      bf.mat <- cbind(bf.mat, Z)
    }
    # set the basis information matrix
    bf.info <- object
  } else {
    stop("Basis functions provided are not of the correct type. See documentation for details")
  }
  # ensure the matrix is sparse if required
  if (bf.matrix.type == "sparse") {
    bf.mat <- methods::as(bf.mat, "sparseMatrix")
  }
  # add in the basis information matrix to the return object
  attr(bf.mat, "bf.df") <- bf.info

  return(bf.mat)
}
