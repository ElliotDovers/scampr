#' Search algorithm for simple 2D basis functions configurations on scampr IDMs
#'
#' @description This function takes in a scampr model and calculates likelihoods and BIC for the list of basis functions supplied. If none are supplied then the algorithm fits increasingly dense regular grids of basis functions (of the type created by scampr::simple_basis). The algorithm starts with an IPP (i.e. zero basis functions) and increases to <= 'max.basis.functions'. If 'po.fold.id' and/or 'pa.fold.id' are supplied then the function will perform a k-fold hold-one-out cross validation to calculate out-of-sample likelihoods (conditional on the latent field).
#'
#' @param object a scampr model: object of class 'scampr' that provides the framework for the search algorithm. Recommended that an IPP model of the appropriate type is used.
#' @param search.rate an integer determining the rate of increasingly dense basis function configurations trialled. Default is \code{search.rate = 1}, however this can be increased to reduce computation time (at the expense of how fine-scale the search will be).
#' @param return.model a logical indicating whether to return the model with the lowest BIC found through the search. Default is \code{FALSE}, meaning the full search results are returned.
#' @param max.basis.functions Optional. An integer describing a rough upper limit to the number of basis functions to search. Defaults to half the number of presences in the data sets.
#' @param radius.type a character string describing the type of radius length to use. One of 'diag' = diagonal dist. between nodes or 'limiting' = sqrt(Domain Area)/log(k).
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to use sparse or dense matrix computations for the basis functions created.
#' @param domain.data Optional. A data frame of columns 'coord.names' that contains at least the extremities of the domain of interest. Useful to ensure the same basis function configurations are created by 'simple_basis' if comparing to various searches.
#' @param start.nodes an integer determining the effective number of basis functions to start the search from (\code{k = start.nodes^2} on a square domain). Default is \code{start.nodes = 4}, however this can be increased so that the search is started from a denser basis function configuration (and will likely increase computation time).
#'
#' @return a data.frame with columns including- 'nodes': number used in scampr::simple_basis to create basis configuration. 'k': the number of basis functions. 'radius': the radius of the basis function configuration. 'k_bias': the number of basis functions used for the PO biasing latent field. 'radius_bias': the radius of the PO biasing basis function configuration.'ll': the fitting marginal log-likelihood. 'BIC': the corresponding Bayesian Info. Crit. 'cpu': the computation time for the model fits. 'convergence': indicator for whether the model converged properly (0 = convergence).#' @export
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
basis.search.idm <- function(object, search.rate = 1, return.model = FALSE, max.basis.functions, radius.type = c("diag", "limiting"), bf.matrix.type = c("sparse", "dense"), domain.data, start.nodes = 4) {

  if (object$model.type != "IDM") {
    stop("Model provided must be of type 'IDM'")
  }

  # use three quarters the number of observations as the maximum k to try if a maximum is not provided
  if (missing(max.basis.functions)) {
    # set the max number of basis functions to half the number of presences
    max.basis.functions <- 0.75 * max(c(sum(attr(object$data, "PA")[ , all.vars(object$formula[[2]])] > 0), sum(object$data[, all.vars(object$formula[[2]])])))
  }

  # get the call
  mc <- match.call() # gets the arguments (must be updated to the alterations above)
  call.list <- as.list(mc)
  # remove the call function
  call.list[[1]] <- NULL
  # ensure the PA search returns the res table
  call.list$return.model = F

  # perform the presence/absence data search
  res_pa <- do.call("basis.search.pa", call.list)

  # extract the last two basis configurations
  bf.list <- attr(res_pa, "bfs")
  if (length(bf.list) - 1 > 1) {
    bf1 <- bf.list[[length(bf.list) - 1]]
    bf2 <- bf.list[[length(bf.list)]]
  } else { # this covers the case where the PA model does not require ANY SRE - could change this once I update the IDM to handle no SRE of the joint model but allow for biasing SRE
    bf1 <- bf.list[[length(bf.list)]]
    bf2 <- bf.list[[length(bf.list)]]
  }

  # update the supplied model
  mA <- do.call("update", list(object, include.sre = T, model.type = "IDM", basis.functions = bf1, latent.po.biasing = T, po.biasing.basis.functions = bf2))
  mB <- do.call("update", list(object, include.sre = T, model.type = "IDM", basis.functions = bf1, latent.po.biasing = T, po.biasing.basis.functions = bf1))

  # compare the models
  direction <- c("denser", "coarser")[which.max(logLik(mA, mB)$logLik)]

  # start a new search based on the direction
  if (direction == "denser") {
    ## Set up initial model comparison to initialise while loop ##
    m1 <- mA
    # initialise storage objects
    tmp.nodes <- res_pa[c(nrow(res_pa) - 1, nrow(res_pa)), "nodes"]
    tmp.k <- res_pa[c(nrow(res_pa) - 1, nrow(res_pa)), "k"]
    tmp.radius <- res_pa[c(nrow(res_pa) - 1, nrow(res_pa)), "radius"]
    tmp.ll <- c(logLik(mB), logLik(m1))
    tmp.bic <- BIC(mB, m1)$BIC
    tmp.time <- c(mB$cpu[1], m1$cpu[1])
    tmp.conv <- c(mB$convergence, m1$convergence)
    basis.functions.list <- list(bf1, bf2)
    # first iteration
    counter <- 2
    tmp.nodes <- c(tmp.nodes, tmp.nodes[[counter]] + search.rate)
    # simple basis function configuration
    tmp.bfs <- simple_basis(tmp.nodes[counter + 1], data = domain.data, radius.type = radius.type)
    basis.functions.list[[counter + 1]] <- tmp.bfs
    m2 <- do.call("update", list(m1, po.biasing.basis.functions = tmp.bfs))

    # store info
    tmp.k <- c(tmp.k, nrow(m2$random.bias.effects))
    tmp.radius <- c(tmp.radius, tmp.bfs$scale[1])
    tmp.ll <- c(tmp.ll, logLik(m2))
    tmp.bic <- c(tmp.bic, BIC(m2))
    tmp.time <- c(tmp.time, m2$cpu[1])
    tmp.conv <- c(tmp.conv, m2$convergence)
    counter <- 3
    while (BIC(m2) < BIC(m1) & nrow(tmp.bfs) < max.basis.functions) {
      # set the current model as the base
      m1 <- m2
      rm(m2)
      # create new basis config.
      tmp.nodes <- c(tmp.nodes, tmp.nodes[[counter]] + search.rate)
      tmp.bfs <- simple_basis(tmp.nodes[counter + 1], data = domain.data, radius.type = radius.type)
      # update the model
      m2 <- do.call("update", list(m1, po.biasing.basis.functions = tmp.bfs))
      # store in basis function list
      basis.functions.list[[counter + 1]] <- tmp.bfs
      # record the number of basis functions, log-likelihood and BIC, etc.
      tmp.k <- c(tmp.k, nrow(m2$random.bias.effects))
      tmp.radius <- c(tmp.radius, tmp.bfs$scale[1])
      tmp.ll <- c(tmp.ll, logLik(m2))
      tmp.bic <- c(tmp.bic, BIC(m2))
      tmp.time <- c(tmp.time, m2$cpu[1])
      tmp.conv <- c(tmp.conv, m2$convergence)
      # increase the counter (loop control variable)
      counter <- counter + 1
    }
  } else {
    ## Set up initial model comparison to initialise while loop ##
    m1 <- mB
    # initialise storage objects
    tmp.nodes <- res_pa[c(nrow(res_pa), nrow(res_pa) - 1), "nodes"]
    tmp.k <- res_pa[c(nrow(res_pa), nrow(res_pa) - 1), "k"]
    tmp.radius <- res_pa[c(nrow(res_pa), nrow(res_pa) - 1), "radius"]
    tmp.ll <- c(logLik(mA), logLik(m1))
    tmp.bic <- BIC(mA, m1)$BIC
    tmp.time <- c(mA$cpu[1], m1$cpu[1])
    tmp.conv <- c(mA$convergence, m1$convergence)
    basis.functions.list <- list(bf2, bf1)
    # first iteration
    counter <- 2
    if (tmp.nodes[[counter]] - search.rate > 1) {
      tmp.nodes <- c(tmp.nodes, tmp.nodes[[counter]] - search.rate)
      # simple basis function configuration
      tmp.bfs <- simple_basis(tmp.nodes[counter + 1], data = domain.data, radius.type = radius.type)
      basis.functions.list[[counter + 1]] <- tmp.bfs
      m2 <- do.call("update", list(m1, po.biasing.basis.functions = tmp.bfs))

      # store info
      tmp.k <- c(tmp.k, nrow(m2$random.bias.effects))
      tmp.radius <- c(tmp.radius, tmp.bfs$scale[1])
      tmp.ll <- c(tmp.ll, logLik(m2))
      tmp.bic <- c(tmp.bic, BIC(m2))
      tmp.time <- c(tmp.time, m2$cpu[1])
      tmp.conv <- c(tmp.conv, m2$convergence)
      counter <- 3
      while (BIC(m2) < BIC(m1) & tmp.nodes[[counter]] - search.rate > 1) {
        # set the current model as the base
        m1 <- m2
        rm(m2)
        # create new basis config.
        tmp.nodes <- c(tmp.nodes, tmp.nodes[[counter]] - search.rate)
        tmp.bfs <- simple_basis(tmp.nodes[counter + 1], data = domain.data, radius.type = radius.type)
        # update the model
        m2 <- do.call("update", list(m1, po.biasing.basis.functions = tmp.bfs))
        # store in basis function list
        basis.functions.list[[counter + 1]] <- tmp.bfs
        # record the number of basis functions, log-likelihood and BIC, etc.
        tmp.k <- c(tmp.k, nrow(m2$random.bias.effects))
        tmp.radius <- c(tmp.radius, tmp.bfs$scale[1])
        tmp.ll <- c(tmp.ll, logLik(m2))
        tmp.bic <- c(tmp.bic, BIC(m2))
        tmp.time <- c(tmp.time, m2$cpu[1])
        tmp.conv <- c(tmp.conv, m2$convergence)
        # increase the counter (loop control variable)
        counter <- counter + 1
      }
    }
  }
  # collate the results
  res <- data.frame(nodes = tmp.nodes,
                    k = res_pa[nrow(res_pa) - 1, "k"],
                    radius = res_pa[nrow(res_pa) - 1, "radius"],
                    k_bias = tmp.k,
                    radius_bias = tmp.radius,
                    ll = tmp.ll,
                    BIC = tmp.bic,
                    cpu = tmp.time,
                    convergence = tmp.conv
  )
  # alter the PA search results to combine
  res_pa$k_bias <- NA
  res_pa$radius_bias <- NA
  res <- rbind(res_pa[ , c("nodes", "k", "radius", "k_bias", "radius_bias", "ll", "BIC", "cpu", "convergence")], res)
  attr(res, "bfs") <- basis.functions.list

  # determine which object is to be returned
  if (return.model) {
    # add in the basis search computation time
    m1$cpu["basis.search"] <- sum(res$cpu) + sum(res_pa$cpu)
    ret.obj <- m1
    attr(ret.obj, "search.res") <- res
  } else {
    ret.obj <- res
  }
  # add in the PA data search
  attr(ret.obj, "pa.search") <- res_pa
  attr(ret.obj, "bfs") <- basis.functions.list
  return(ret.obj)
}
