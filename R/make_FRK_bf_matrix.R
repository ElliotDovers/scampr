#' Create an FRK basis function matrix
#'
#' @param data
#' @param coord.names
#' @param manifold
#' @param regular
#' @param nres
#' @param prune
#' @param max_basis
#' @param subsamp
#' @param type
#' @param isea3h_lo
#' @param bndary
#' @param scale_aperture
#' @param verbose
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
make_FRK_bf_matrix <- function(data, coord.names = c("x", "y"), manifold = FRK::plane(), regular = 1, nres = 2, prune = 0,
                           max_basis = NULL, subsamp = 10000, type = c("bisquare", "Gaussian",
                                                                       "exp", "Matern32"), isea3h_lo = 2, bndary = NULL,
                           scale_aperture = ifelse(is(manifold, "sphere"), 1, 1.25), verbose = 0L,
                           ...) {

  # checks #
  if (!all(coord.names %in% colnames(data))) {
    stop(paste0("coord.names, ", coord.names, " not found in the data provided"))
  }

  # create a spatial pixels data frame as required by FRK::auto_basis
  sp.data <- data
  coordinates(sp.data) <- coord.names

  # Use the FRK auto_basis function
  bf.as.functions <- FRK::auto_basis(manifold = manifold, data = sp.data, regular = regular, nres = nres, prune = prune,
                         max_basis = max_basis, subsamp = subsamp, type = type, isea3h_lo = isea3h_lo, bndary = bndary,
                         scale_aperture = scale_aperture, verbose = verbose)

  # Evaluate the basis functions at the data points
  bf.list <- lapply(bf.as.functions@fn, function(x){x(data[,coord.names])})

  # Combine into a matrix and make sparse
  Z <- as(as.matrix(do.call(cbind, bf.list)), "sparseMatrix")

  # Add an attribute that is a vector of length # basis functions indicating which resolution each belongs
  attr(Z, "resolution_id") <- bf.as.functions@df$res

  return(Z)
}
