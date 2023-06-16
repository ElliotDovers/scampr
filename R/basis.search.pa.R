#' Search algorithm for simple 2D basis functions configurations on scampr PA models
#'
#' @description This function takes in a scampr model and fits increasingly dense regular grids of basis functions (of the type created by \code{scampr::simple_basis}) to find an optimal basis function configuration according to either log-likelihood, AIC or BIC.
#'
#' @param object a scampr model: object of class 'scampr' that provides the framework for the search algorithm. It is recommended that the model does not include spatial random effects to save the computational burden of fitting such a model first.
#' @param metric a character string describing the metric upon which to choose the optimal basis function configuration. One of 'll' (log-Likelihood), 'aic' (Akaike Information Criterion), 'bic' (Bayesian Information Criterion).
#' @param return.model a logical indicating whether to return the model with optimal basis function configuration according to \code{metric}. Default is \code{TRUE}, meaning search results are returned are attached to the model as an attribute, i.e. accessible via \code{attr(., "search.res")}.
#' @param max.basis.functions Optional. An integer describing a rough upper limit to the number of basis functions to search. Defaults to half the number of presences in the data sets.
#' @param radius.type a character string describing the type of radius length to use. One of 'diag' = diagonal dist. between nodes or 'limiting' = sqrt(Domain Area)/log(k).
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to use sparse or dense matrix computations for the basis functions created.
#' @param domain.data Optional. A data frame of columns 'coord.names' that contains at least the extremities of the domain of interest. Useful to ensure the same basis function configurations are created by 'simple_basis' if comparing to various searches.
#' @param start.nodes an integer determining the effective number of basis functions to start the search from (\code{k = start.nodes^2} on a square domain). Default is \code{start.nodes = 4}, however this can be increased so that the search is started from a denser basis function configuration (this will likely increase computation time).
#' @param search.rate an integer determining the rate of increasingly dense basis function configurations trialled. Default is \code{search.rate = 1}, however this can be increased to reduce computation time (at the expense of how fine-scale the search will be).
#' @param metric.tol a numeric describing the tolerance level for the search stopping rule. Specifically, the proportion of the metric (calculated from \code{object}).
#' @param lag an integer determining the lag/window length for the moving average of the selection metric. Default is 3.
#'
#' @return Depends on \code{return.model}: If \code{TRUE} then the \code{scampr} model fitted with the optimised basis function configuration. If \code{FALSE} then a data.frame with columns including- 'nodes': number used in scampr::simple_basis to create basis configuration. 'k': the number of basis functions. 'radius': the radius of the basis function configuration. 'll': the fitting marginal log-likelihood. 'BIC': the corresponding Bayesian Info. Crit. 'cpu': the computation time for the model fits. 'convergence': indicator for whether the model converged properly (0 = convergence).
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
#' # perform search
#' res <- basis.search(m.ipp)
#' }
basis.search.pa <- function(object, metric = c("ll", "aic", "bic"), return.model = TRUE, max.basis.functions, radius.type = c("diag", "limiting"), bf.matrix.type = c("sparse", "dense"), domain.data, start.nodes = 4, search.rate = 1, metric.tol = 0, lag = 3) {

  if (object$model.type == "PO") {
    stop("Model provided must be of type 'PA' or 'IDM'")
  }

  # checks not covered by model fitting
  metric <- match.arg(metric)
  radius.type <- match.arg(radius.type)
  bf.matrix.type <- match.arg(bf.matrix.type)

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
  if (object$model.type == "IDM") {
    base.m <- do.call("update", list(object, formula = object$formula, include.sre = F, model.type = "PA", data = attr(object$data, "PA")))
  } else {
    base.m <- do.call("update", list(object, include.sre = F, model.type = "PA"))
  }

  ## Set up initial model comparison to initialise while loop ##

  # initialise storage objects
  mod.list <- list()
  basis.functions.list <- list()
  tmp.nodes <- 0
  tmp.k <- 0
  tmp.radius <- NA
  tmp.ll <- logLik(base.m)
  tmp.aic <- AIC(base.m)
  tmp.bic <- BIC(base.m)
  tmp.time <- base.m$cpu[1]
  tmp.conv <- base.m$convergence
  basis.functions.list[[1]] <- NA

  # first iteration
  m1 <- base.m
  # store model in list
  mod.list[[1]] <- m1
  tmp.nodes <- c(tmp.nodes, start.nodes)
  # simple basis function configuration
  tmp.bfs <- simple_basis(tmp.nodes[2], data = domain.data, radius.type = radius.type)
  basis.functions.list[[2]] <- tmp.bfs
  m2 <- do.call("update", list(base.m, include.sre = T, basis.functions = tmp.bfs, starting.pars = base.m))
  # store model in list
  mod.list[[2]] <- m2
  # store info
  tmp.k <- c(tmp.k, sum(m2$basis.per.res))
  tmp.radius <- c(tmp.radius, tmp.bfs$scale[1])
  tmp.ll <- c(tmp.ll, logLik(m2))
  tmp.aic <- c(tmp.aic, AIC(m2))
  tmp.bic <- c(tmp.bic, BIC(m2))
  tmp.time <- c(tmp.time, m2$cpu[1])
  tmp.conv <- c(tmp.conv, m2$convergence)
  counter <- 2
  # determine the metrics for algorithm stopping
  met1 <- switch(metric,
                 ll = -tmp.ll[1:(length(tmp.ll) - 1)],
                 aic = min(tmp.aic[1:(length(tmp.aic) - 1)]),
                 bic = min(tmp.bic[1:(length(tmp.bic) - 1)])
  )
  met2 <- switch(metric,
                 ll = -tmp.ll[length(tmp.ll)],
                 aic = tmp.aic[length(tmp.aic)],
                 bic = tmp.bic[length(tmp.bic)]
  )
  # calculate the tolerance level from the base model
  base.tol <- switch(metric,
                 ll = abs(-logLik(base.m)) * metric.tol,
                 aic = abs(AIC(base.m)) * metric.tol,
                 bic = abs(BIC(base.m)) * metric.tol
  )
  # create the stopping condition
  continue.if <- if (length(met1) < lag) {
    met2 < met1[length(met1)] + base.tol | met2 < mean(met1) + base.tol
  } else {
    met2 < met1[length(met1)] + base.tol | met2 < mean(met1[length(met1):(length(met1) - (lag - 1))]) + base.tol
  }
  while (continue.if & nrow(tmp.bfs) < max.basis.functions) {
    # set the current model as the base
    m1 <- m2
    rm(m2)
    # create new basis config.
    tmp.nodes <- c(tmp.nodes, tmp.nodes[[counter]] + search.rate)
    tmp.bfs <- simple_basis(tmp.nodes[counter + 1], data = domain.data, radius.type = radius.type)
    # update the model
    m2 <- do.call("update", list(base.m, include.sre = T, basis.functions = tmp.bfs, starting.pars = base.m))
    # store model in list
    mod.list[[counter + 1]] <- m2
    # store in basis function list
    basis.functions.list[[counter + 1]] <- tmp.bfs
    # record the number of basis functions, log-likelihood and BIC, etc.
    tmp.k <- c(tmp.k, sum(m2$basis.per.res))
    tmp.radius <- c(tmp.radius, tmp.bfs$scale[1])
    tmp.ll <- c(tmp.ll, logLik(m2))
    tmp.aic <- c(tmp.aic, AIC(m2))
    tmp.bic <- c(tmp.bic, BIC(m2))
    tmp.time <- c(tmp.time, m2$cpu[1])
    tmp.conv <- c(tmp.conv, m2$convergence)
    # determine the metrics for algorithm stopping
    met1 <- switch(metric,
                   ll = -tmp.ll[1:(length(tmp.ll) - 1)],
                   aic = min(tmp.aic[1:(length(tmp.aic) - 1)]),
                   bic = min(tmp.bic[1:(length(tmp.bic) - 1)])
    )
    met2 <- switch(metric,
                   ll = -tmp.ll[length(tmp.ll)],
                   aic = tmp.aic[length(tmp.aic)],
                   bic = tmp.bic[length(tmp.bic)]
    )
    # re-calculate the stopping condition
    continue.if <- if (length(met1) < lag) {
      met2 < met1[length(met1)] + base.tol | met2 < mean(met1) + base.tol
    } else {
      met2 < met1[length(met1)] + base.tol | met2 < mean(met1[length(met1):(length(met1) - (lag - 1))]) + base.tol
    }
    # increase the counter (loop control variable)
    counter <- counter + 1
  }

  # get the optimised model
  opt.mod.id <- length(tmp.ll) - 1

  res <- data.frame(nodes = tmp.nodes,
                    k = tmp.k,
                    radius = tmp.radius,
                    ll = tmp.ll,
                    AIC = tmp.aic,
                    BIC = tmp.bic,
                    cpu = tmp.time,
                    convergence = tmp.conv
  )

  # add in the model selected
  res$opt <- (1:nrow(res)) == opt.mod.id

  # determine which object is to be returned
  if (return.model) {
    ret.obj <- mod.list[[opt.mod.id]]
    # add in the basis search computation time
    ret.obj$cpu["basis.search"] <- sum(res$cpu)
    attr(ret.obj, "search.res") <- res
  } else {
    ret.obj <- res
  }
  attr(ret.obj, "bfs") <- basis.functions.list
  return(ret.obj)
}
