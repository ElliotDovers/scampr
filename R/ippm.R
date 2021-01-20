#' Inhomogeneous Poisson Process Model for Point Patterns
#'
#' @description Fit an inhomogeneous Poisson process model to presence records using numerical quadrature
#'
#' @param formula formula describing fixed effects of the linear predictor. Response must be the presence/quadrature identifier.
#' @param data data frame containing predictors at both presence-records and quadrature
#' @param starting.pars optional named list or scampr model that gives warm starts for the parameters of the model
#' @param se logical indicating whether standard errors should be calculated
#' @param coord.names vector of character strings describing the column names of the coordinates in data
#' @param quad.weights.name charater string of the column name of quadrature weights in data
#'
#' @return scampr model object
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
#' m.ipp <- ippm(pres ~ elev.std, data = dat)
ippm <- function(formula, data, coord.names = c("x", "y"), quad.weights.name = "quad.size", starting.pars, se = TRUE) {

  # Use the model for presence-only data with particular parameters hard-coded
  mod <- po(formula, data, coord.names = coord.names, quad.weights.name = quad.weights.name, model.type = "ipp", starting.pars = starting.pars, se = se)

  return(mod)
}
