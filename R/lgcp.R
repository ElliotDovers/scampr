#' Approximate log-Gaussian Cox Process Models for Point Patterns
#'
#' @description Fit a log-Gaussian Cox process model to a point pattern using numerical quadrature (provided with the data, see e.g. scampr::gorillas) and one of either a Laplace or variational approximation to marginalise over the latent field.
#'
#' @param formula formula describing fixed effects of the linear predictor. Response must be the presence/quadrature binary (1/0) column name.
#' @param data data frame containing predictors at both presence-records and quadrature
#' @param coord.names vector of character strings describing the column names of the coordinates in data
#' @param quad.weights.name  charater string of the column name of quadrature weights in data
#' @param FRK.basis.functions object of class 'Basis' from FRK package. If neither 'FRK.basis.functions' nor 'simple.basis' is specified will use default FRK::auto_basis with 2 spatial resolutions
#' @param simple.basis Alternative to 'FRK.basis.functions', a data.frame of basis functions information created by 'simple_basis()'.
#' @param approx.with charater indicating the type of approximation to use for the intractible marginalisation: variational or Laplace
#' @param se logical indicating whether standard errors should be calculated
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to use sparse or dense matrix computations for the basis functions created.
#' @param starting.pars optional named list or scampr model that gives warm starts for the parameters of the model
#' @param subset an optional vector describing a subset of the data to be used.
#'
#' @return a scampr model object.
#' @export
#'
#' @examples
#' #' # Get the gorilla nesting data
#' dat <- gorillas
#'
#' # Standardise the elevation covariate
#' dat$elev.std <- scale(dat$elevation)
#'
#' # Set up a simple 2D grid of basis functions to fit a LGCP model to the data
#' bfs <- simple_basis(nodes.on.long.edge = 9, data = dat)
#'
#' # Fit a LGCP model using variational approximation
#' m.lgcp_va <- lgcp(pres ~ elev.std, data = dat, approx.with = "variational", simple.basis = bfs)
#'
#' \dontrun{
#' # Fit a LGCP model using Laplace approximation
#' m.lgcp_lp <- lgcp(pres ~ elev.std, data = dat, approx.with = "laplace", simple.basis = bfs)
#' }
lgcp <- function(formula, data, coord.names = c("x", "y"), quad.weights.name = "quad.size", FRK.basis.functions, simple.basis, approx.with = c("laplace", "variational"), se = TRUE, bf.matrix.type = c("sparse", "dense"), starting.pars, subset) {

  approx.with <- match.arg(approx.with)
  bf.matrix.type <- match.arg(bf.matrix.type)
  # Use the model for presence-only data with particular parameters hard-coded
  mod <- po(formula, data, coord.names = coord.names, quad.weights.name = quad.weights.name, FRK.basis.functions = FRK.basis.functions, simple.basis = simple.basis, model.type = approx.with, bf.matrix.type = bf.matrix.type, se = se, starting.pars = starting.pars, subset = subset)

  return(mod)
}
