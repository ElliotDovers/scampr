#' Search algorithm for simple 2D basis functions configurations on scampr PO models
#'
#' @description This function takes in a scampr model and calculates likelihoods and AIC for the list of basis functions supplied. If none are supplied then the algorithm fits increasingly dense regular grids of basis functions (of the type created by scampr::simple_basis). The algorithm starts with an IPP (i.e. zero basis functions) and increases to <= 'max.basis.functions'. If 'po.fold.id' and/or 'pa.fold.id' are supplied then the function will perform a k-fold hold-one-out cross validation to calculate out-of-sample likelihoods (conditional on the latent field).
#'
#' @param object a scampr model: object of class 'scampr' that provides the framework for the search algorithm. Recommended that an IPP model of the appropriate type is used.
#' @param search.rate an integer determining the rate of increasingly dense basis function configurations trialled. Default is \code{search.rate = 1}, however this can be increased to reduce computation time (at the expense of how fine-scale the search will be).
#' @param return.model a logical indicating whether to return the model with the lowest BIC found through the search. Default is \code{FALSE}, meaning the full search results are returned.
#' @param max.basis.functions Optional. An integer describing a rough upper limit to the number of basis functions to search. Defaults to half the number of presences in the data sets.
#' @param radius.type a character string describing the type of radius length to use. One of 'diag' = diagonal dist. between nodes or 'limiting' = sqrt(Domain Area)/log(k).
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to use sparse or dense matrix computations for the basis functions created.
#' @param domain.data Optional. A data frame of columns 'coord.names' that contains at least the extremities of the domain of interest. Useful to ensure the same basis function configurations are created by 'simple_basis' if comparing to various searches.
#' @param start.nodes an integer determining the effective number of basis functions to start the search from (\code{k = start.nodes^2} on a square domain). Default is \code{start.nodes = 4}, however this can be increased so that the search is started from a denser basis function configuration (and will likely increase computation time).
#' @param which.approx a character string indicating the type of approximation to be used to marginalise over the spatial random effects. May be one of 'laplace' or 'variational' - the latter should result in faster search times.
#'
#' @return a data.frame with columns including- 'nodes': number used in scampr::simple_basis to create basis configuration. 'k': the number of basis functions. 'radius': the radius of the basis function configuration. 'll': the fitting marginal log-likelihood. 'BIC': the corresponding Bayesian Info. Crit. 'cpu': the computation time for the model fits. 'convergence': indicator for whether the model converged properly (0 = convergence).#' @export
#' @export
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
basis.search.po <- function(object, search.rate = 1, return.model = FALSE, max.basis.functions, radius.type = c("diag", "limiting"), bf.matrix.type = c("sparse", "dense"), domain.data, start.nodes = 4, which.approx = c("variational", "laplace")) {

  if (object$model.type == "PA") {
    stop("Model provided must be of type 'PO' or 'IDM'")
  }

  # checks not covered by model fitting
  radius.type <- match.arg(radius.type)
  bf.matrix.type <- match.arg(bf.matrix.type)
  which.approx <- match.arg(which.approx)

  # Use provided model data as domain.data if missing and check coords are present
  if (missing(domain.data)) {
    domain.data <- rbind(object$data[ , object$coord.names], attr(object$data, "PA")[ , object$coord.names])
  } else {
    if (!all(object$coord.names %in% colnames(domain.data))) {
      stop(paste0("Model coord.names, ", object$coord.names, ", not found in 'domain.data' frame provided"))
    }
  }
  # use three quarters the number of observations as the maximum k to try if a maximum is not provided
  if (missing(max.basis.functions)) {
    # set the max number of basis functions to half the number of presences
    max.basis.functions <- 0.75 * max(c(sum(attr(object$data, "PA")[ , all.vars(object$formula[[2]])] > 0), sum(object$data[, all.vars(object$formula[[2]])])))
  }

  # create the base Binomial Model without SRE
  base.m <- do.call("update", list(object, include.sre = F, model.type = "PO", sre.approx = which.approx))

  ## Set up initial model comparison to initialise while loop ##

  # initialise storage objects
  tmp.nodes <- 0
  tmp.k <- 0
  tmp.radius <- NA
  tmp.ll <- logLik(base.m)
  tmp.bic <- BIC(base.m)
  tmp.time <- base.m$cpu[1]
  tmp.conv <- base.m$convergence
  basis.functions.list <- list()
  # first iteration
  m1 <- base.m
  tmp.nodes <- c(tmp.nodes, start.nodes)
  # simple basis function configuration
  tmp.bfs <- simple_basis(tmp.nodes[2], data = domain.data, radius.type = radius.type)
  basis.functions.list[[1]] <- tmp.bfs
  m2 <- do.call("update", list(base.m, include.sre = T, basis.functions = tmp.bfs, starting.pars = base.m))

  # store info
  tmp.k <- c(tmp.k, sum(m2$basis.per.res))
  tmp.radius <- c(tmp.radius, tmp.bfs$scale[1])
  tmp.ll <- c(tmp.ll, logLik(m2))
  tmp.bic <- c(tmp.bic, BIC(m2))
  tmp.time <- c(tmp.time, m2$cpu[1])
  tmp.conv <- c(tmp.conv, m2$convergence)
  counter <- 2
  while (BIC(m2) < BIC(m1) & nrow(tmp.bfs) < max.basis.functions) {
    # set the current model as the base
    m1 <- m2
    rm(m2)
    # create new basis config.
    tmp.nodes <- c(tmp.nodes, tmp.nodes[[counter]] + search.rate)
    tmp.bfs <- simple_basis(tmp.nodes[counter + 1], data = domain.data, radius.type = radius.type)
    # update the model
    m2 <- do.call("update", list(base.m, include.sre = T, basis.functions = tmp.bfs, starting.pars = base.m))
    # store in basis function list
    basis.functions.list[[counter]] <- tmp.bfs
    # record the number of basis functions, log-likelihood and BIC, etc.
    tmp.k <- c(tmp.k, sum(m2$basis.per.res))
    tmp.radius <- c(tmp.radius, tmp.bfs$scale[1])
    tmp.ll <- c(tmp.ll, logLik(m2))
    tmp.bic <- c(tmp.bic, BIC(m2))
    tmp.time <- c(tmp.time, m2$cpu[1])
    tmp.conv <- c(tmp.conv, m2$convergence)
    # increase the counter (loop control variable)
    counter <- counter + 1
  }

  res <- data.frame(nodes = tmp.nodes,
                    k = tmp.k,
                    radius = tmp.radius,
                    ll = tmp.ll,
                    BIC = tmp.bic,
                    cpu = tmp.time,
                    convergence = tmp.conv
  )

  # determine which object is to be returned
  if (return.model) {
    # add in the basis search computation time
    m1$cpu["basis.search"] <- sum(res$cpu)
    ret.obj <- m1
    attr(ret.obj, "search.res") <- res
  } else {
    ret.obj <- res
  }
  attr(ret.obj, "bfs") <- basis.functions.list
  return(ret.obj)
}
