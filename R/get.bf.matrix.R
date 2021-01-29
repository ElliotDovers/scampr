#' Internal scampr function that calculates the basis function evaluations at point locations provided. From either a fitted model or data.frame of class 'bf.df'
#'
#' @param object either a fitted model or data.frame of class 'bf.df' created by 'simple_basis'
#' @param point.locations either a matrix or data.frame describing the point locations in the same order as in 'object'
#'
#' @return
#' @noRd
#'
#' @importFrom FRK eval_basis
#' @importFrom fields rdist
#' @importFrom methods as
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
#' m.lgcp_va <- po(pres ~ elev.std, data = dat, model.type = "variational", simple.basis = bfs)
#'
#' # Get some new locations
#' new.sites <- dat[sample(1:nrow(dat), 10), c("x", "y")]
#'
#' # Calculate the basis function matrix at the new sites
#' scampr:::get.bf.matrix(m.lgcp_va, new.sites)
get.bf.matrix <- function(object, point.locations) {
  # Check the object is of one of the two correct types
  if(!"bf.df" %in% class(object)) {
    if (class(object)[1] != "scampr") {
      stop("'object' must be a model of class 'scampr' or a basis function data.frame of class 'bf.df'")
    }
  }
  # If the object is a scampr model created without an FRK auto_basis use the basis function info data frame
  object.switch <- T
  if (class(object)[1] == "scampr") {
    if (is.na(object$approx.type)) {
      stop("Cannot get a basis function matrix from an IPP model")
    } else {
      if (!is.null(object$FRK.basis.functions)) {
        bf.mat <- FRK::eval_basis(object$FRK.basis.functions, as.matrix(point.locations))
        attr(bf.mat, "bf.df") <- object$basis.fn.info
        object.switch <- F
      } else {
        tmp <- object$basis.fn.info
        object <- tmp
      }
    }
  }
  if (object.switch) {
    bf.mat <- NULL
    for (res in unique(object$res)) {
      radius <- object$scale[object$res == res][1]
      dist.mat <- fields::rdist(point.locations, object[,1:2][object$res == res, ])
      Z <- matrix(0, nrow = nrow(dist.mat), ncol = ncol(dist.mat))
      Z[dist.mat <= radius] <- (1 - (dist.mat[dist.mat <= radius] / radius)^2)^2
      bf.mat <- cbind(bf.mat, Z)
    }
    bf.mat <- methods::as(bf.mat, "sparseMatrix")
    attr(bf.mat, "bf.df") <- object
  }
  return(bf.mat)
}
