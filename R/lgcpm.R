#' Approximate log-Gaussian Cox Process Model for a Point Pattern
#'
#' @description Blah blah blah
#'
#' @param formula formula describing fixed effects of the linear predictor. Response must be the presence/quadrature binary (1/0) column name.
#' @param data data frame containing predictors at both presence-records and quadrature
#' @param coord.names vector of character strings describing the column names of the coordinates in data
#' @param quad.weights.name  charater string of the column name of quadrature weights in data
#' @param FRK.basis.functions object of class 'Basis' from FRK package. If neither 'FRK.basis.functions' nor 'simple.basis' is specified will use default FRK::auto_basis with 2 spatial resolutions
#' @param simple.basis Alternative to 'FRK.basis.functions', a data.frame of basis functions information created by 'simple_basis()'.
#' @param approx.with charater indicating the type of approximation to use for the intractible marginalisation: variational or Laplace
#' @param se logical indicating whether standard errors should be calculated
#' @param starting.pars optional named list or scampr model that gives warm starts for the parameters of the model
#' @param subset optional subset of the data
#' @param na.action optional way of handling NA's in the data, default is omit
#'
#' @return scampr model object
#' @export
#'
#' @examples
lgcpm <- function(formula, data, coord.names = c("x", "y"), quad.weights.name = "quad.size", FRK.basis.functions, simple.basis, approx.with = c("laplace", "variational"), se = TRUE, starting.pars, subset, na.action) {

  approx.with <- match.arg(approx.with)
  # Use the model for presence-only data with particular parameters hard-coded
  mod <- po(formula, data, coord.names = coord.names, quad.weights.name = quad.weights.name, FRK.basis.functions = FRK.basis.functions, simple.basis = simple.basis, model.type = approx.with, starting.pars = starting.pars, se = se, subset = subset, na.action = na.action)

  return(mod)
}
