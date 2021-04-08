#' Search algorithim for simple 2D basis functions configurations on scampr models
#'
#' @description This function takes in a scampr model (ipp, lgcp, po, pa, popa) and calculates likelihoods and AIC for increasingly dense regular grids of basis functions (of the type created by scampr::simple_basis). The algorithm starts with an IPP (i.e. zero basis functions) and increases to <= 'max.basis.functions'. If 'po.fold.id' and/or 'pa.fold.id' are supplied then the function will perform a k-fold hold-one-out cross validation to calculate out-of-sample likelihoods (conditional on the latent field).
#'
#' @param scampr.model a scampr model: object of class 'scampr' that provides the framework for the search alogrithm. Recommended that an IPP model of the appropriate type is used.
#' @param po.fold.id Optional for cross-validation. An integer or factor vector the same length as the po.data that describes the CV fold that each location falls into.
#' @param pa.fold.id Optional for cross-validation. An integer or factor vector the same length as the pa.data that describes the CV fold that each location falls into.
#' @param max.basis.functions Optional. An integer describing a rough upper limit to the number of basis functions to search. Defaults to half the number of presences in the data sets.
#' @param radius.type a character string describing the type of radius length to use. One of 'diag' = diagonal dist. between nodes or 'limiting' = sqrt(Domain Area)/log(k).
#' @param domain.data Optional. A data frame of columns 'coord.names' that contains at least the extremities of the domain of interest. Useful to ensure the same basis function configurations are created by 'simple_basis' if comparing to various searches.
#' @param loglik.type a charater string indicating the type of approximation to use for the intractible marginalisation. One of 'variational', 'laplace' or 'both'. Currently variational (and hence both) is only available to PO models.
#' @param trunc.pa.prob Optional. A small positve number by which the predicted probability of presence is truncated. This can be used to ensure infinite values are avoiding within the cross-validation.
#'
#' @return a data.frame with columns including- 'nodes.on.long.edge': number used in scampr::simple_basis to create basis configuration. 'bf': the number of basis functions. 'loglik': the fitting marginal log-likelihood. 'aic': the corresponding AIC. Optionally, 'predicted_cll_po': the conditional (on the latent field) Presence-only likelihood. 'predicted_cll_pa': the conditional (on the latent field) Presence/Absence likelihood. 'roc_auc': Area under the ROC curve on the Presence/Absence data. Optional columns are the results from a cross-validation described by 'po.fold.id' and/or 'pa.fold.id'. (_va or _lp subscript for approx. type if both are calculated).
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
simple_basis_search <- function(scampr.model, po.fold.id, pa.fold.id, max.basis.functions, radius.type = c("diag", "limiting"), domain.data, loglik.type = c("laplace", "variational", "both"), trunc.pa.prob = 1e-7) {

  # Determine model and search scenario
  mod.type <- c("po", "pa", "popa")[which(scampr.model$data.model.type == c("po", "pa", "popa"))]
  ll.type <- match.arg(loglik.type)
  switch.id <- paste(mod.type, ll.type, sep = ".")
  result <- switch(switch.id,
                   po.laplace = simple_basis_search_po(po.formula = scampr.model$formula, po.data = scampr.model$data, po.fold.id, max.basis.functions = max.basis.functions, approx.with = "laplace", coord.names = scampr.model$coord.names, quad.weights.name = scampr.model$quad.weights.name, radius.type = radius.type, bf.matrix.type = "sparse", domain.data = domain.data),
                   pa.laplace = simple_basis_search_pa(pa.formula = scampr.model$formula, pa.data = scampr.model$data, pa.fold.id, max.basis.functions = max.basis.functions, coord.names = scampr.model$coord.names, radius.type = radius.type, bf.matrix.type = "sparse", domain.data, approx.with = "laplace", trunc.pa.prob = trunc.pa.prob),
                   popa.laplace = simple_basis_search_popa(po.formula = scampr.model$formula, pa.formula = attr(scampr.model$formula, "pa"), po.data = scampr.model$data, pa.data = attr(scampr.model$data, "pa"), po.fold.id, pa.fold.id, max.basis.functions = max.basis.functions, coord.names = scampr.model$coord.names, quad.weights.name = scampr.model$quad.weights.name, radius.type = radius.type, bf.matrix.type = "sparse", domain.data = domain.data, approx.with = "laplace", trunc.pa.prob = trunc.pa.prob),
                   po.variational = simple_basis_search_po(po.formula = scampr.model$formula, po.data = scampr.model$data, po.fold.id, max.basis.functions = max.basis.functions, approx.with = "variational", coord.names = scampr.model$coord.names, quad.weights.name = scampr.model$quad.weights.name, radius.type = radius.type, bf.matrix.type = "sparse", domain.data = domain.data),
                   pa.variational = stop("No compatibility for 'variational' approx. with PA models"),
                   popa.variational =  stop("No compatibility for 'variational' approx. with POPA models"),
                   po.both = simple_basis_search_po_both_approx(po.formula = scampr.model$formula, po.data = scampr.model$data, po.fold.id, max.basis.functions = max.basis.functions, coord.names = scampr.model$coord.names, quad.weights.name = scampr.model$quad.weights.name, radius.type = radius.type, bf.matrix.type = "sparse", domain.data = domain.data),
                   pa.both = stop("No compatibility for 'variational' approx. with PA models"),
                   popa.both = stop("No compatibility for 'variational' approx. with POPA models")
  )
  return(result)
}
