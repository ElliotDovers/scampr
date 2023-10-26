#' Internal scampr function that gets warm starting parameters to pass to TMB
#'
#' @param new.start.pars either a fitted model or named list of starting parameters
#' @param old.start.pars a named list of starting parameters
#' @param target.approx.type a character string, one of 'laplace', 'variational' or 'not_sre'.
#'
#' @return a data.frame (sparse or dense depending on parameter bf.matrix.type)
#' @export
#'
#' @examples
#' # Get the gorilla nesting data
#' data(gorillas, package = "scampr")
#' dat <- gorillas
#'
#' # Standardise the elevation covariate
#' dat$elev.std <- scale(dat$elevation)
#'
#' # Fit a scampr model to the point pattern
#' m <- scampr(pres ~ elev.std, data = dat, include.sre = F)
#'
#' # Get new start parameters
#' update.startingp.parameters(m, list(fixed = rep(0, ncol(dat.list$X_PO_pres)),
#'   random = rep(0, ncol(dat.list$Z_PO_pres)),
#'   log_variance_component = var.starts))
update.starting.parameters <- function(new.start.pars, old.start.pars, target.approx.type = c("laplace", "variational", "not_sre")) {
  target.approx.type <- match.arg(target.approx.type)
  # check if the new starting parameters are a fitted scampr model
  if (is(new.start.pars, "scampr")) {
    # collect the model object
    tmp.m <- new.start.pars
    new.start.pars <- lapply(split(tmp.m$par, names(tmp.m$par)), unname)
    # check the model isn't an IPP
    if (tmp.m$approx.type != "not_sre") {
      if (target.approx.type == "laplace" & tmp.m$approx.type == "variational") {
        # make appropriate change to the variance parameter if going from VA to Laplace
        new.start.pars$log_variance_component <- unname(log(sqrt(tmp.m$prior.variances)))
      } else if (target.approx.type == "variational" & tmp.m$approx.type == "laplace") {
        # make appropriate change to the variance parameter going from Laplace to VA
        new.start.pars$log_variance_component <- NULL
      }
      # need to add the random parameters if the existing model is laplace
      if (tmp.m$approx.type == "laplace") {
        new.start.pars$random <- unname(tmp.m$random.effects[, 1L])
      }
    }
    rm(tmp.m)
  }

  # loop through the parameter names to match and replace
  for (n in names(new.start.pars)) {
    # check the name is in the new starting parameter list
    if (n %in% names(old.start.pars)) {
      # check that the lengths are correct
      if (length(new.start.pars[[n]]) != length(old.start.pars[[n]])) {
        warning(paste0("'", n, "' starting parameters not used - the number provided does not match the proposed model."))
      } else {
        old.start.pars[[n]] <- new.start.pars[[n]]
      }
    }
  }

  return(old.start.pars)
}
