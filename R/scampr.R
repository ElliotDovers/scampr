#' Spatially Correlated, Aprroximate Modelling of Presence Records
#'
#' @description Blah blah blah
#'
#' @param formula formula describing fixed effects of the linear predictor. Response must be the presence/quadrature identifier.
#' @param data data frame containing predictors at both presence-records and quadrature
#' @param coord.names vector of character strings describing the column names of the coordinates in data
#' @param quad.weights.name  charater string of the column name of quadrature weights in data
#' @param bf.matrix matrix of basis functions evaulated at the locations in data
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
scampr <- function(formula, data, coord.names = c("x", "y"), quad.weights.name = "quad.size", bf.matrix = NULL, approx.with = c("laplace", "variational"), se = TRUE, starting.pars = NULL, subset, na.action) {

  ## checks ##
  # coordinate names exist in the data
  if (!all(coord.names %in% colnames(data))) {
    stop(paste0("coord.names, ", coord.names, " not found in the data provided"))
  }
  # quadrature weights column exists in the data
  if (!quad.weights.name %in% colnames(data)) {
    stop(paste0("quad.weights.name, ", quad.weights.nam, " not found in the data provided"))
  }
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
  # approx.with is one of those available
  match.arg(approx.with)
  ############################################################

  # Need to adjust for different input (presence pts. and quadrature)?

  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(x = c("formula", "data", "subset", "na.action"),
             table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(expr = mf, envir = parent.frame())

  mt <- attr(x = mf, which = "terms")
  des.mat <- model.matrix(object = mt, data = mf, contrasts.arg = contrasts)
  pt.quad.id <- mf[ , 1L]
  fixed.names <- colnames(des.mat)

  # When a basis function matrix is not provided use defaults
  if (is.null(bf.matrix)) {
    bf.matrix <- make_FRK_bf_matrix(data = data, coord.names = coord.names)
  }

  # based on the type of starting.pars adjust to req. list structure
  if (!is.null(starting.pars)) {
    if (typeof(starting.pars) == "character") {
      tmp.m <- ipp(formula, data)
      starting.pars <- lapply(split(tmp.m$par, names(tmp.m$par)), unname)
      rm(tmp.m)
    } else if (class(starting.pars) == "scampr") {
      starting.pars <- lapply(split(starting.pars$par, names(starting.pars$par)), unname)
    }
  }

  # TMB required data setup
  dat.list <- list(
    X_pres = des.mat[pt.quad.id == 1, ],
    Z_pres = bf.matrix[pt.quad.id == 1, ],
    X_quad = des.mat[pt.quad.id == 0, ],
    Z_quad = bf.matrix[pt.quad.id == 0, ],
    quad_size = data[ , quad.weights.name][pt.quad.id == 0],
    bf_per_res = as.numeric(table(attr(bf.matrix, "resolution_id")))
  )

  # Link C++ file
  # Hopefully I can store the pre-compiled .dll/.so files (so OS dictates use) in PACKAGE/inst
  # Otherwise compile on package install? This would require TMB pre installed
  path2file <- paste0(path.package("scampr"), "/", TMB::dynlib(paste0("scampr_", approx.with[1], "_loglik")))
  dyn.load(path2file)

  if (approx.with[1] == "variational") {
    start.pars <- list(fixed = rep(0, ncol(dat.list$X_pres)), random = rep(0, ncol(bf.matrix)), logPosteriorVar = rep(0, ncol(bf.matrix)))
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
    obj <- TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = paste0("scampr_", approx.with[1], "_loglik"), silent = T)
  } else {
    start.pars <- list(fixed = rep(0, ncol(dat.list$X_pres)), random = rep(0, ncol(bf.matrix)), logPriorSD = rep(0, length(dat.list$bf_per_res)))
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
    obj <- TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = paste0("scampr_", approx.with[1], "_loglik"), silent = T)
  }
  # optimise the parameters
  res <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS", control = list(maxit = 1000))
  # get standard errors if required
  if (se) {
    tmp.estimates <- summary(TMB::sdreport(obj))
  } else {
    tmp.estimates <- cbind(Estimate = res$par, `Std. Error` = rep(NA, length(res$par)))
  }

  # unload the function lib
  # dyn.unload(TMB::dynlib(paste0("scampr_", approx.with[1], "_loglik")))

  # add required information to the results list
  res$fixed.effects <- tmp.estimates[1:length(fixed.names), ]
  res$random.effects <- tmp.estimates[(length(fixed.names) + 1):nrow(tmp.estimates), ]
  row.names(res$fixed.effects) <- fixed.names
  res$starting.pars <- start.pars
  # a.list <- lapply(split(res$par, names(res$par)), unname)
  res$data <- data
  res$formula <- formula
  res$approx_type <- approx.with
  res$basis_per_res <- dat.list$bf_per_res
  res$coordinates <- coord.names
  class(res) <- "scampr"
  return(res)
}
