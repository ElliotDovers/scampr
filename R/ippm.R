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
ippm <- function(formula, data, coord.names = c("x", "y"), quad.weights.name = "quad.size", starting.pars, se = TRUE) {

  ## checks ##
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
  des.mat <- model.matrix(object = mt, data = mf)
  pt.quad.id <- mf[ , 1L]
  fixed.names <- colnames(des.mat)

  # if starting.pars provided is a scampr model adjust to req. list structure
  if (!missing(starting.pars)) {
    if (class(starting.pars) == "scampr") {
      starting.pars <- lapply(split(starting.pars$par, names(starting.pars$par)), unname)
    }
  }
  # TMB required data setup
  dat.list <- list(
    X_pres = des.mat[pt.quad.id == 1, ],
    Z_pres = as(matrix(rep(0, 4), nrow = 2), "sparseMatrix"),
    X_quad = des.mat[pt.quad.id == 0, ],
    Z_quad = as(matrix(rep(0, 4), nrow = 2), "sparseMatrix"),
    quad_size = data[ , quad.weights.name][pt.quad.id == 0],
    bf_per_res = 2,
    mod_type = 0
  )
  # obtain warm starts for parameters if provided
  start.pars <- list(fixed = rep(0, ncol(dat.list$X_pres)), random = rep(0, ncol(dat.list$Z_pres)), log_variance_component = rep(0, ncol(dat.list$Z_pres)))
  if (!missing(starting.pars)) {
    for (n in names(starting.pars)) {
      if (n %in% names(start.pars)) {
        if (length(starting.pars[[n]]) != length(start.pars[[n]])) {
          stop(paste0("The number of '", n, "' starting parameters provided does not match the proposed model"))
        }
        start.pars[[n]] <- starting.pars[[n]]
      }
    }
  }
  # set up the objective function w.r.t. approx.with
  obj <- TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", map = list(random = factor(rep(NA, ncol(dat.list$Z_pres))), log_variance_component = factor(rep(NA, ncol(dat.list$Z_pres)))), silent = T)
  # optimise the parameters
  res <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS", control = list(maxit = 1000))
  # get standard errors if required
  if (se) {
    tmp.estimates <- summary(TMB::sdreport(obj))
  } else {
    tmp.estimates <- cbind(Estimate = res$pars, `Std. Error` = rep(NA, length(res$pars)))
  }
  # add required information to the results list
  res$coefficients <- res$par
  names(res$coefficients) <- fixed.names
  res$fixed.effects <- tmp.estimates[1:length(fixed.names), ]
  res$random.effects <- NA
  row.names(res$fixed.effects) <- fixed.names
  res$starting.pars <- start.pars
  res$data <- data
  res$fitted.values <- as.vector(des.mat %*% res$fixed.effects[ , 1])
  res$formula <- formula
  res$coord.names <- coord.names
  res$quad.weights.name <- quad.weights.name
  res$pt.quad.id <- pt.quad.id
  res$approx.type <- NA
  res$basis.per.res <- NA
  res$FRK.basis.functions <- NULL
  res$basis.fn.info <- NULL
  class(res) <- "scampr"
  return(res)
}
