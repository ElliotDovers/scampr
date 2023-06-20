#' Search algorithm for simple 2D basis functions configurations on scampr models
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
#' @param which.approx a character string indicating the type of approximation to be used to marginalise over the spatial random effects. May be one of 'laplace' or 'variational' - the latter should result in faster search times.
#' @param start.nodes an integer determining the effective number of basis functions to start the search from (\code{k = start.nodes^2} on a square domain). Default is \code{start.nodes = 4}, however this can be increased so that the search is started from a denser basis function configuration (and will likely increase computation time).
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
#' # Search through an increasingly dense regular grid of basis functions
#' res <- basis.search(m.ipp)
#' }
basis.search <- function(object, metric = c("ll", "aic", "bic"), return.model = TRUE, max.basis.functions, radius.type = c("diag", "limiting"), bf.matrix.type = c("sparse", "dense"), domain.data, which.approx = c("variational", "laplace"), start.nodes = 4, search.rate = 1, metric.tol = 0, lag = 3) {

  call.list <- as.list(match.call())
  call.list[[1]] <- NULL

  if (object$model.type == "IDM") {
    # search basis functions on the PA data to inform the shared field approximation
    res_pa <- do.call("basis.search.pa", call.list)
    # search basis functions on the PO data to inform the PO biasing field approximation
    res_po <- do.call("basis.search.po", call.list)
    if (return.model) {
      # fit the optimised IDM
      res <- do.call("update", list(object,
                             include.sre = if (res_pa$approx.type == "not_sre") {FALSE} else {TRUE},
                             basis.functions = res_pa$basis.functions,
                             latent.po.biasing = if (res_po$approx.type == "not_sre") {FALSE} else {TRUE},
                             po.biasing.basis.functions = res_po$basis.functions,
                             starting.pars = object))
      # collate the different searches to return
      attr(res_pa, "search.res")$data = "PA"
      attr(res_po, "search.res")$data = "PO"
      attr(res, "search.res") <- rbind(attr(res_pa, "search.res"), attr(res_po, "search.res"))
    } else {
      # collate the different searches to return
      res_pa$data = "PA"
      res_po$data = "PO"
      res <- rbind(res_pa, res_po)
    }

  } else {
    res <- switch(object$model.type,
           PA = do.call("basis.search.pa", call.list),
           PO = do.call("basis.search.po", call.list)
    )
  }
  return(res)
}
