#' Inhomogeneous Poisson Process Model for Point Patterns
#'
#' @description Fit an inhomogeneous Poisson process model to a point pattern using numerical quadrature (provided with the data, see e.g. scampr::gorillas)
#'
#' @param formula formula describing fixed effects of the linear predictor. Response must be the presence/quadrature identifier.
#' @param data data frame containing predictors at both presence-records and quadrature
#' @param starting.pars optional named list or scampr model that gives warm starts for the parameters of the model
#' @param se logical indicating whether standard errors should be calculated
#' @param coord.names vector of character strings describing the column names of the coordinates in data
#' @param quad.weights.name charater string of the column name of quadrature weights in data
#' @param subset an optional vector describing a subset of the data to be used.
#'
#' @return scampr model object
#' @export
#'
#' @examples
#' # Get the gorilla nesting data
#' dat <- gorillas
#'
#' # Standardise the elevation covariate
#' dat$elev.std <- scale(dat$elevation)
#'
#' # Fit an IPP model to the point pattern
#' m.ipp <- ipp(pres ~ elev.std, data = dat)
ipp <- function(formula, data, coord.names = c("x", "y"), quad.weights.name = "quad.size", starting.pars, se = TRUE, subset) {

  # Use the model for presence-only data with particular parameters hard-coded
  mod <- po(formula, data, coord.names = coord.names, quad.weights.name = quad.weights.name, model.type = "ipp", starting.pars = starting.pars, se = se, subset = subset)

  return(mod)
}
