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
#' dat <- gorillas
update.starting.parameters <- function(new.start.pars, old.start.pars, target.approx.type = c("laplace", "variational", "not_sre")) {
  target.approx.type <- match.arg(target.approx.type)
  # check if the new starting parameters are a fitted scampr model
  if (is(new.start.pars, "scampr")) {
    # collect the model object
    tmp.m <- new.start.pars
    new.start.pars <- lapply(split(tmp.m$par, names(tmp.m$par)), unname)
    # check the model isn't an IPP
    if (!is.na(tmp.m$approx.type)) {
      if (target.approx.type == "laplace" & tmp.m$approx.type == "variational") {
        # make appropriate change to the variance parameter if going from VA to Laplace
        new.start.pars$log_variance_component <- unname(log(sqrt(tmp.m$random.effects[grepl("Prior Var ", rownames(tmp.m$random.effects), fixed = T), 1L])))
      } else if (target.approx.type == "variational" & tmp.m$approx.type == "laplace") {
        # make appropriate change to the variance parameter going from Laplace to VA
        new.start.pars$log_variance_component <- NULL
      }
      # need to add the random parameters if the existing model is laplace
      if (tmp.m$approx.type == "laplace") {
        new.start.pars$random <- unname(tmp.m$random.effects[grepl("LP Posterior Mean", rownames(tmp.m$random.effects), fixed = T), 1L])
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
        stop(paste0("The number of '", n, "' starting parameters provided does not match the proposed model"))
      }
      old.start.pars[[n]] <- new.start.pars[[n]]
    }
  }

  return(old.start.pars)
}
