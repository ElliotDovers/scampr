#' Search algorithm for simple 2D basis functions configurations on scampr models
#'
#' @description This function takes in a scampr model and calculates likelihoods and AIC for the list of basis functions supplied. If none are supplied then the algorithm fits increasingly dense regular grids of basis functions (of the type created by scampr::simple_basis). The algorithm starts with an IPP (i.e. zero basis functions) and increases to <= 'max.basis.functions'. If 'po.fold.id' and/or 'pa.fold.id' are supplied then the function will perform a k-fold hold-one-out cross validation to calculate out-of-sample likelihoods (conditional on the latent field).
#'
#' @param data a data frame of columns 'coord.names' that contains at least the extremities of the domain of interest. Useful to ensure the same basis function configurations are created by 'simple_basis' if comparing to various searches.
#' @param max.basis.functions Optional. An integer describing a rough upper limit to the number of basis functions to search. Defaults to half the number of presences in the data sets.
#' @param coord.names vector of character strings describing the names of the coordinates in 'data'. Ordered horizontal axis to vertical
#' @param radius.type a character string describing the type of radius length to use. One of 'diag' = diagonal dist. between nodes or 'limiting' = sqrt(Domain Area)/log(k).
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to use sparse or dense matrix computations for the basis functions created.
#'
#' @return a list of basis function configurations up to the maximum supplied
#' @noRd
#'
#' @examples
#' #' # Get the gorilla nesting data
#' dat <- gorillas
#'
#' # Standardise the elevation covariate
#' dat$elev.std <- scale(dat$elevation)
#'
#' # Fit an IPP model to the point pattern
#' m.ipp <- scampr(pres ~ elev.std, data = dat, model.type = "ipp")
#'  \dontrun{
#' # Search through an increasingly dense regular grid of basis functions
#' res <- simple_basis_search(m.ipp)
#' }
candidate.bfs <- function(data, max.basis.functions, coord.names = c("x", "y"), radius.type = c("diag", "limiting"), bf.matrix.type = c("sparse", "dense")) {
  # checks
  radius.type <- match.arg(radius.type)
  bf.matrix.type <- match.arg(bf.matrix.type)
  if (!all(coord.names %in% colnames(data))) {
    stop("at least one of 'coord.names' not found in the data provided")
  }

  if (missing(max.basis.functions)) { # if a maximum number of basis functions is not provided
    # set the max number of basis functions to one fifth the number rows in the data provided
    max.basis.functions <- 0.2 * nrow(data)
  }
  # initialise objects
  tmp.nodes <- NULL
  tmp.k <- NULL
  basis.functions.list <- list()
  # first iteration
  tmp.nodes <- c(tmp.nodes, 2)
  # simple basis function configuration
  tmp.bfs <- simple_basis(tmp.nodes[1], data = data, radius.type = radius.type)
  basis.functions.list[[1]] <- tmp.bfs
  tmp.nodes <- c(tmp.nodes, 3)
  tmp.k <- c(tmp.k, nrow(tmp.bfs))
  counter <- 2
  while (nrow(tmp.bfs) < max.basis.functions) {
    # create new basis config.
    tmp.bfs <- simple_basis(tmp.nodes[counter], data = data, radius.type = radius.type)
    # store in basis function list
    basis.functions.list[[counter]] <- tmp.bfs
    # record the number of basis functions
    tmp.k <- c(tmp.k, nrow(tmp.bfs))

    # set the next number of nodes to use
    tmp.nodes <- c(tmp.nodes, tmp.nodes[[counter]] + 1)
    # increase the counter (loop control variable)
    counter <- counter + 1
  }
  return(basis.functions.list)
}

