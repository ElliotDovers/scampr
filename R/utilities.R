# Utilities

scampr2startpars <- function(scampr.model, for.type = c("laplace", "variational", "ipp")) {
  starting.pars <- lapply(split(scampr.model$par, names(scampr.model$par)), unname)
  # check the model isn't an IPP
  if (!is.na(scampr.model$approx.type)) {
    if (for.type == "laplace" & scampr.model$approx.type == "variational") {
      # make appropriate change to the variance parameter if going from VA to Laplace
      starting.pars$log_variance_component <- unname(log(sqrt(scampr.model$random.effects[grepl("Prior Var ", rownames(scampr.model$random.effects), fixed = T), 1])))
    } else if (for.type == "variational" & scampr.model$approx.type == "laplace") {
      # make appropriate change to the variance parameter going from Laplace to VA
      starting.pars$log_variance_component <- NULL
    }
    # need to add the random parameters if the existing model is laplace
    if (scampr.model$approx.type == "laplace") {
      starting.pars$random <- unname(scampr.model$random.effects[grepl("LP Posterior Mean", rownames(scampr.model$random.effects), fixed = T), 1L])
    }
  }
  return(starting.pars)
}

get.single.model.aic <- function(object, k = 2) {
  if (class(object) != "scampr") {
    stop(paste0(deparse(substitute(object)), " is not a scampr model"))
  }
  # get an identifier for the model type
  mod.id <- NULL
  if (is.na(object$approx.type)) {
    mod.id <- "ipp"
  } else {
    mod.id <- object$approx.type
  }
  # Need to get the random effect coefficients from Laplace models
  add.coef <- switch(mod.id,
                     ipp = 0,
                     variational = length(object$basis.per.res) -2 * sum(object$basis.per.res), # adjusts for posterior pars included in $coefficients
                     laplace = length(object$basis.per.res)
  )
  aic <- -2*logLik.scampr(object) + k*length(object$coefficients) + k*add.coef
  return(aic)
}
