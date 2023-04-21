#' Search algorithm for simple 2D basis functions configurations on scampr PA models
#'
#' @description This function takes in a scampr model and calculates likelihoods and AIC for the list of basis functions supplied. If none are supplied then the algorithm fits increasingly dense regular grids of basis functions (of the type created by scampr::simple_basis). The algorithm starts with an IPP (i.e. zero basis functions) and increases to <= 'max.basis.functions'. If 'po.fold.id' and/or 'pa.fold.id' are supplied then the function will perform a k-fold hold-one-out cross validation to calculate out-of-sample likelihoods (conditional on the latent field).
#'
#' @param pres a scampr model: object of class 'scampr' that provides the framework for the search algorithm. Recommended that an IPP model of the appropriate type is used.
#' @param quad an integer determining the rate of increasingly dense basis function configurations trialled. Default is \code{search.rate = 1}, however this can be increased to reduce computation time (at the expense of how fine-scale the search will be).
#' @param formula a logical indicating whether to return the model with the lowest BIC found through the search. Default is \code{FALSE}, meaning the full search results are returned.
#' @param coord.names a vector of character strings describing the column names of the coordinates in both data frames.
#' @param quad.weights.name a character string of the column name of quadrature weights in the data.
#' @param quad.size a numeric indicating the common quadrature size of the full quadrature points provided in \code{quad}.
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to use sparse or dense matrix computations for the basis functions created.
#' @param domain.data Optional. A data frame of columns 'coord.names' that contains at least the extremities of the domain of interest. Useful to ensure the same basis function configurations are created by 'simple_basis' if comparing to various searches.
#' @param start.nodes an integer determining the effective number of basis functions to start the search from (\code{k = start.nodes^2} on a square domain). Default is \code{start.nodes = 4}, however this can be increased so that the search is started from a denser basis function configuration (and will likely increase computation time).
#'
#' @return a data.frame with columns including- 'nodes': number used in scampr::simple_basis to create basis configuration. 'k': the number of basis functions. 'radius': the radius of the basis function configuration. 'll': the fitting marginal log-likelihood. 'BIC': the corresponding Bayesian Info. Crit. 'cpu': the computation time for the model fits. 'convergence': indicator for whether the model converged properly (0 = convergence).
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
make.scampr.data <- function(pres, quad, formula = ~1, coord.names = c("x", "y"), quad.weights.name = "quad.size", smallest.proportion = 0.01, return.model = FALSE, max.basis.functions, radius.type = c("diag", "limiting"), bf.matrix.type = c("sparse", "dense"), domain.data, start.nodes = 4) {

  if ((!all(coord.names %in% colnames(pres)))) {
    stop(paste0("'coord.names', ", coord.names, ", not found in 'pres'"))
  }
  if ((!all(coord.names %in% colnames(quad)))) {
    stop(paste0("'coord.names', ", coord.names, ", not found in 'quad'"))
  }
  if ((!quad.weights.name %in% colnames(quad))) {
    warning(paste0("'quad.weights.name', ", quad.weights.name, ", not found in 'quad'. Defaulting to 'quad.size' assumed to be 1 distance unit^2."))
    eval(parse(text = paste0("quad$", quad.weights.name, " <- 1")))
  }
  if ((!quad.weights.name %in% colnames(pres))) {
    eval(parse(text = paste0("pres$", quad.weights.name, " <- 0")))
  }
  # get the response provided
  resp <- all.vars(formula[[2]])
  if (length(resp) == 0) {
    resp <- "presence"
    formula <- eval(parse(text = paste0("update(formula, ", resp, " ~ .)")))
    update(formula)
  }
  if ((!resp %in% colnames(quad))) {
    eval(parse(text = paste0("quad$", resp, " <- 0")))
  }
  if ((!resp %in% colnames(pres))) {
    eval(parse(text = paste0("pres$", resp, " <- 1")))
  }

  lls <- NULL
  for (i in 1:10) {
    # take an initial sample
    q.id <- sample(1:nrow(quad), nrow(quad) * smallest.proportion)
    subquad <- quad[q.id, ]
    # merge the presence records and quadrature points on common names
    dat <- merge(pres, subquad, by = intersect(names(pres), names(subquad)), all.x = T)
    # fit the model
    m <- scampr(formula, data = dat, model.type = "PO", include.sre = F)
    lls[i] <- logLik(m)
  }
  # take an initial sample
  q.id <- sample(1:nrow(quad), nrow(quad) * smallest.proportion)
  subquad <- quad[q.id, ]
  # merge the presence records and quadrature points on common names
  dat <- merge(pres, subquad, by = intersect(names(pres), names(subquad)), all.x = T)
  # fit the model
  m <- scampr(formula, data = dat, model.type = "PO", include.sre = F)



  if (object$model.type == "PO") {
    stop("Model provided must be of type 'PA' or 'IDM'")
  }
  get.design.matrix(formula, data)
  # checks not covered by model fitting
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
