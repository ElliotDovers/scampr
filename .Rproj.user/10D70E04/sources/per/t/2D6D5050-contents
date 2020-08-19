#' Inhomogeneous Poisson Process
#'
#' @description Fit an inhomogeneous Poisson process model to presence records using numerical quadrature
#'
#' @param formula formula describing fixed effects of the linear predictor. Response must be the presence/quadrature identifier.
#' @param data data frame containing predictors at both presence-records and quadrature
#' @param starting.pars optional named list or scampr model that gives warm starts for the parameters of the model
#' @param se logical indicating whether standard errors should be calculated
#'
#' @return scampr model object
#' @export
#'
#' @examples
ipp <- function(formula, data, starting.pars = NULL, se = TRUE) {

  ## checks ##
  # starting parameters are of the correct type
  if (!is.null(starting.pars)) {
    if (typeof(starting.pars) == "character") {
      if (starting.pars != "ipp") {
        stop("starting.pars must be one of 'ipp', a scampr model object, a named list or NULL")
      }
    } else if (typeof(starting.pars) != "list") {
      stop("starting.pars must be one of 'ipp', a scampr model object, a named list or NULL")
    }
  }
  # se is of the correct type
  if (!is.logical(se)) {
    stop(paste0("'se' must be a logcial indicating whether or not to calculate standard erros"))
  }
  ############################################################

  # Need to adjust for different input (presence pts. and quadrature)

  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(x = c("formula", "data", "subset", "weights", "na.action", "offset"),
             table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(expr = mf, envir = parent.frame())

  mt <- attr(x = mf, which = "terms")
  des.mat <- model.matrix(object = mt, data = mf, contrasts.arg = contrasts)
  pt.quad.id <- mf[ , 1L]
  fixed.names <- colnames(des.mat)

  # TMB Setup
  dat.list <- list(
    X_pres = des.mat[pt.quad.id == 1, ],
    X_quad = des.mat[pt.quad.id == 0, ],
    quad_size = dat$quad.size[pt.quad.id == 0] # Need to work this out
  )

  # Link C++ file
  # Hopefully I can store the pre-compiled .dll/.so files (so OS dictates use) in PACKAGE/inst
  # Otherwise compile on package install? This would require TMB pre installed
  path2file <- paste0(path.package("scampr"), "/", TMB::dynlib("ipp"))
  dyn.load(path2file)

  start.pars <- list(fixed = rep(0, ncol(des.mat)))
  if (!is.null(starting.pars)) {
    for (n in names(starting.pars)) {
      if (n %in% names(start.pars)) {
        if (length(starting.pars[[n]]) != length(start.pars[[n]])) {
          stop(paste0("The number of '", n, "' starting parameters provided does not match the proposed model"))
        }
        start.pars[[n]] <- starting.pars[[n]]
      }
    }
  }
  obj <- TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "ipp", silent = T)
  res <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS", control = list(maxit = 1000))
  if (se) {
    tmp.estimates <- summary(TMB::sdreport(obj))
  } else {
    tmp.estimates <- cbind(Estimate = res$pars, `Std. Error` = rep(NA, length(res$pars)))
  }
  res$fixed.effects <- tmp.estimates[1:length(fixed.names), ]
  res$random.effects <- NA
  row.names(res$fixed.effects) <- fixed.names
  res$starting.pars <- start.pars
  res$data <- data
  res$formula <- formula
  res$approx_type <- NA
  res$basis_per_res <- NA
  class(res) <- "scampr"
  return(res)
}
