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

  ## checks ##
  # coordinate names exist in the data
  if (!all(coord.names %in% colnames(data))) {
    stop(paste0("coord.names, ", coord.names, " not found in the data provided"))
  }
  # quadrature weights column exists in the data
  if (!quad.weights.name %in% colnames(data)) {
    stop(paste0("quad.weights.name, ", quad.weights.nam, " not found in the data provided"))
  }
  # se is of the correct type
  if (!is.logical(se)) {
    stop(paste0("'se' must be a logcial indicating whether or not to calculate standard erros"))
  }
  # approx.with is one of those available
  approx.with <- match.arg(approx.with)
  ############################################################

  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(x = c("formula", "data", "subset", "na.action"),
             table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(expr = mf, envir = parent.frame())

  mt <- attr(x = mf, which = "terms")
  des.mat <- model.matrix(object = mt, data = mf)
  pt.quad.id <- mf[ , 1L]
  fixed.names <- colnames(des.mat)

  # When a basis function matrix is not provided use defaults (with 2 spatial resolutions)
  if (missing(simple.basis)) {
    if (missing(FRK.basis.functions)) {
      # create a spatial pixels data frame as required by FRK::auto_basis
      sp.data <- data
      coordinates(sp.data) <- coord.names
      FRK.basis.functions <- FRK::auto_basis(data = sp.data, nres = 2)

    }
    bf.matrix <- FRK::eval_basis(basis = FRK.basis.functions, as.matrix(data[ , coord.names]))
    bf.info <- FRK.basis.functions@df
    colnames(bf.info)[grepl("loc", colnames(bf.info), fixed = T)] <- coord.names
    class(bf.info) <- c(class(bf.info), "bf.df")
  } else {
    bf.matrix <- scampr:::get.bf.matrix(simple.basis, data[ , coord.names])
    bf.info <- simple.basis
    FRK.basis.functions <- NULL
  }

  # if starting.pars provided is a scampr model adjust to req. list structure
  if (!missing(starting.pars)) {
    if (class(starting.pars) == "scampr") {
      starting.pars <- lapply(split(starting.pars$par, names(starting.pars$par)), unname)
    }
  }
  # TMB required data setup
  dat.list <- list(
    X_pres = des.mat[pt.quad.id == 1, ],
    Z_pres = as(bf.matrix[pt.quad.id == 1, ], "sparseMatrix"),
    X_quad = des.mat[pt.quad.id == 0, ],
    Z_quad = as(bf.matrix[pt.quad.id == 0, ], "sparseMatrix"),
    quad_size = data[ , quad.weights.name][pt.quad.id == 0],
    bf_per_res = as.numeric(table(bf.info$res)),
    mod_type = as.integer(which(approx.with[1] == c("ipp", "variational", "laplace")) - 1)
  )
  # create the appropraite start parameters for the variance component w.r.t. approx. type
  var.starts <- switch(approx.with,
         variational = rep(0, ncol(dat.list$Z_pres)),
         laplace = rep(0, length(dat.list$bf_per_res))
         )

  # initilise starting parameters at 0
  start.pars <- list(fixed = rep(0, ncol(dat.list$X_pres)), random = rep(0, ncol(dat.list$Z_pres)), log_variance_component = var.starts)
  # obtain warm starts for parameters if provided
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
  obj <- switch(approx.with,
                variational = TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", silent = T),
                laplace = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", silent = T)
  )
  # optimise the parameters
  res <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS", control = list(maxit = 1000))
  # get standard errors if required
  if (se) {
    tmp.estimates <- summary(TMB::sdreport(obj))
  } else {
    tmp.estimates <- cbind(Estimate = res$par, `Std. Error` = rep(NA, length(res$par)))
  }
  # get the random component names
  random.nos <- NULL
  if (length(dat.list$bf_per_res) == 1L) {
    random.nos <- 1L:dat.list$bf_per_res
  } else {
    for (lvl in 1L:length(dat.list$bf_per_res)) {
      random.nos <- c(random.nos, paste(lvl, 1L:dat.list$bf_per_res[lvl], sep = "."))
    }
  }

  # add required information to the results list
  res$coefficients <- res$par
  coef.names <- switch(approx.with,
         variational = c(fixed.names, paste0("VA Posterior Mean (bf ", random.nos, ")"), paste0("VA Posterior log sd (bf ", random.nos, ")")),
         laplace = c(fixed.names, paste0("Prior log sd (res. ", 1:length(dat.list$bf_per_res), ")"))
  )
  names(res$coefficients) <- coef.names
  res$fixed.effects <- tmp.estimates[1:length(fixed.names), ]
  rownames(res$fixed.effects) <- fixed.names
  res$random.effects <- tmp.estimates[(length(fixed.names) + 1):nrow(tmp.estimates), ]
  res$random.effects <- res$random.effects[!grepl("log_", rownames(res$random.effects), fixed = T), ]
  rand.names <- switch(approx.with,
                       variational = c(paste0("VA Posterior Mean (bf ", random.nos, ")"), paste0("VA Posterior Var (bf ", random.nos, ")"), paste0("Prior Var (res. ", 1:length(dat.list$bf_per_res), ")")),
                       laplace = c(paste0("LP Posterior Mean (bf ", random.nos, ")"), paste0("Prior Var (res. ", 1:length(dat.list$bf_per_res), ")"))
  )
  rownames(res$random.effects) <- rand.names
  res$starting.pars <- start.pars
  res$data <- data
  res$fitted.values <- as.vector(des.mat %*% res$fixed.effects[ , 1] + bf.matrix %*% res$random.effects[grepl(" Mean ", rownames(res$random.effects), fixed = T), 1])
  res$formula <- formula
  res$coord.names <- coord.names
  res$quad.weights.name <- quad.weights.name
  res$pt.quad.id <- pt.quad.id
  res$approx.type <- approx.with
  res$basis.per.res <- dat.list$bf_per_res
  res$FRK.basis.functions <- FRK.basis.functions
  res$basis.fn.info <- bf.info
  class(res) <- "scampr"
  return(res)
}
