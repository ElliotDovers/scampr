#' Point Process Model for Presence-only Data
#'
#' @description Fit either an inhomogeneous Poisson, or log-Gaussian Cox, process model to presence records. This uses numerical quadrature  (provided with the data, see e.g. scampr::gorillas) to approx. the spatial integral. If fitting a LGCP, uses one of either a Laplace or variational approximation to marginalise over the latent field.
#'
#' @param po.formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the data model to be fitted. The 'response' must be a binary that indicates whether a datum is a presence or quadrature point. See GLM function for further formula details.
#' @param po.data a data frame containing predictors at both presence-records and quadrature as well as the po.formula 'response'.
#' @param coord.names a vector of character strings describing the column names of the coordinates in both data frames
#' @param quad.weights.name a charater string of the column name of quadrature weights in the po.data
#' @param FRK.basis.functions an optional object of class 'Basis' from FRK package. If neither 'FRK.basis.functions' nor 'simple.basis' is specified will use default FRK::auto_basis with 2 spatial resolutions.
#' @param simple.basis an alternative to 'FRK.basis.functions': a data.frame of basis functions information created by 'simple_basis()'.
#' @param model.type a character string indicating the type of model to be used. May be one of 'laplace' or 'variational' for Cox Processes involving spatially correlated errors or 'ipp' for a model that follows an inhomgeneous Poisson process.
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to use sparse or dense matrix computations for the basis functions created.
#' @param se a logical indicating whether standard errors should be calculated.
#' @param starting.pars an optional named list or scampr model object that gives warm starting values for the parameters of the model.
#' @param subset an optional vector describing a subset of the data to be used.
#'
#' @return a scampr model object
#' @noRd
#'
#' @importFrom methods as
#' @importFrom stats optim
#' @importFrom sp coordinates
#' @importFrom FRK auto_basis eval_basis
#' @importFrom TMB MakeADFun sdreport
#'
#' @examples
#' # Get the gorilla nesting data
#' dat <- gorillas
#'
#' # Standardise the elevation covariate
#' dat$elev.std <- scale(dat$elevation)
#'
#' # Fit an IPP model to the point pattern
#' m.ipp <- po(pres ~ elev.std, po.data = dat, model.type = "ipp")
#'
#' # Set up a simple 2D grid of basis functions to fit a LGCP model to the data
#' bfs <- simple_basis(nodes.on.long.edge = 9, data = dat)
#'
#' # Fit a LGCP model using variational approximation
#' m.lgcp_va <- po(pres ~ elev.std, po.data = dat, model.type = "variational", simple.basis = bfs)
#'
#' \dontrun{
#' # Fit a LGCP model using Laplace approximation
#' m.lgcp_lp <- po(pres ~ elev.std, po.data = dat, model.type = "laplace", simple.basis = bfs)
#' }
po <- function(po.formula, po.data, coord.names = c("x", "y"), quad.weights.name = "quad.size", FRK.basis.functions, simple.basis, model.type = c("variational", "laplace", "ipp"), bf.matrix.type = c("sparse", "dense"), se = TRUE, starting.pars, subset) {

  # CAN'T JUST GIVE A BASIS FUNCTION MATRIX BECAUSE THEN YOU CAN'T PREDICT ETC. AS WE DON'T KNOW ENOUGH ABOUT THE FUNCTIONS

  # Get the response and predictor names
  po.resp <- all.vars(po.formula[[2]])
  po.pred <- all.vars(po.formula[[3]])

  ## checks ##

  if (length(po.resp) != 1) {
    stop("Formula can only take a single response")
  }
  if (!all(c(po.resp, po.pred) %in% colnames(po.data))) {
    stop("data does not contain the formula terms")
  }
  if (!all(coord.names %in% colnames(po.data))) {
    stop(paste0("coord.names, ", coord.names, ", not found in data set provided"))
  }
  if (!quad.weights.name %in% colnames(po.data)) {
    stop(paste0("quad.weights.name, ", quad.weights.name, ", not found in the data provided"))
  }
  if (!is.logical(se)) {
    stop(paste0("'se' must be a logcial indicating whether or not to calculate standard errors"))
  }
  # parameters of restricted strings
  model.type <- match.arg(model.type)
  bf.matrix.type <- match.arg(bf.matrix.type)

  ############################################################

  # Apply subsetting to data if supplied (with checks)
  if (!missing(subset)) {
    if (!is.vector(subset)) {
      stop("subset must be a vector")
    } else {
      if (class(subset) == "logical") {
        if (length(subset) != nrow(po.data)) {
          stop("Logical subset must be of same dimension as data provided")
        } else {
          po.data <- po.data[subset, ]
        }
      } else if (class(subset) == "integer" | class(subset) == "numeric") {
        if (!all(subset %in% 1:nrow(po.data))) {
          stop("numeric subset must include row numbers of the data provided")
        }
        po.data <- po.data[subset, ]
      } else {
        stop("subset has an incorrect format")
      }
    }
  }

  # default na.action is to remove any data rows with na (for terms involved in the model)
  rm.rows <- attr(na.omit(po.data[ , c(coord.names, quad.weights.name, po.resp, po.pred)]), "na.action")
  if (!is.null(rm.rows)) {
    po.data <- po.data[-rm.rows, ]
  }

  # Get the design matrix
  po.des.mat <- get.design.matrix(po.formula, po.data)
  fixed.names <- colnames(po.des.mat)

  # Get the presence point/ quadrature point identifier (use get.design.matrix to replicate subsets, etc.)
  pt.quad.id <- po.data[ , po.resp]

  # Set the fixed effect names
  fixed.names <- colnames(po.des.mat)

  # Determine the basis functions to be used
  if (model.type != "ipp") {
    if (missing(simple.basis)) {
      if (missing(FRK.basis.functions)) { # When none is provided use FRK defaults (with 2 spatial resolutions)
        # create a spatial pixels data frame as required by FRK::auto_basis
        sp.data <- po.data[ ,coord.names]
        sp::coordinates(sp.data) <- coord.names
        FRK.basis.functions <- FRK::auto_basis(data = sp.data, max_basis = sum(pt.quad.id)*0.25)
      }
      po.bf.matrix <- FRK::eval_basis(basis = FRK.basis.functions, as.matrix(po.data[ , coord.names]))
      bf.info <- FRK.basis.functions@df
      colnames(bf.info)[grepl("loc", colnames(bf.info), fixed = T)] <- coord.names
      class(bf.info) <- c(class(bf.info), "bf.df")
    } else { # Otherwise use the provided simple basis
      po.bf.matrix <- get.bf.matrix(simple.basis, po.data[ , coord.names])
      bf.info <- simple.basis
      FRK.basis.functions <- NULL
    }
  } else {
    po.bf.matrix <- matrix(rep(0, nrow(po.des.mat)), ncol = 1)
    bf.info <- cbind.data.frame(x = NA, y = NA, scale = NA, res = 1)
  }

  # if starting.pars provided is a scampr model adjust to req. list structure
  if (!missing(starting.pars)) {
    if (class(starting.pars) == "scampr") {
      tmp.m <- starting.pars
      starting.pars <- lapply(split(starting.pars$par, names(starting.pars$par)), unname)
      # check the model isn't an IPP
      if (!is.na(tmp.m$approx.type)) {
        # make appropriate change to the variance parameter if going from VA to Laplace
        if (model.type == "laplace" & tmp.m$approx.type == "variational") {
          starting.pars$log_variance_component <- unname(log(sqrt(tmp.m$random.effects[grepl("Prior Var ", rownames(tmp.m$random.effects), fixed = T), 1])))
        }
        # make appropriate change to the variance parameter going from Laplace to VA
        if (model.type == "variational" & tmp.m$approx.type == "laplace") {
          starting.pars$log_variance_component <- NULL
        }
        # need to add the random parameters if the existing model is laplace
        if (tmp.m$approx.type == "laplace") {
          starting.pars$random <- unname(tmp.m$random.effects[grepl("LP Posterior Mean", rownames(tmp.m$random.effects), fixed = T), 1L])
        }
      }
      rm(tmp.m)
    }
  }
  # TMB required data setup
  dat.list <- list(
    X_PO_pres = as.matrix(po.des.mat[pt.quad.id == 1, ]),
    B_PO_pres = matrix(rep(0, sum(pt.quad.id == 1)), ncol = 1),
    X_PO_quad = as.matrix(po.des.mat[pt.quad.id == 0, ]),
    B_PO_quad = matrix(rep(0, sum(pt.quad.id == 0)), ncol = 1),
    X_PA = matrix(0, ncol = 1),
    Z_PO_pres = if(bf.matrix.type == "sparse") {
      methods::as(po.bf.matrix[pt.quad.id == 1, ], "sparseMatrix")
    } else {
      as.matrix(po.bf.matrix[pt.quad.id == 1, ])
    },
    Z_PO_quad = if(bf.matrix.type == "sparse") {
      methods::as(po.bf.matrix[pt.quad.id == 0, ], "sparseMatrix")
    } else {
      as.matrix(po.bf.matrix[pt.quad.id == 0, ])
    },
    Z_PA = if(bf.matrix.type == "sparse") {
      methods::as(matrix(0, ncol = 1), "sparseMatrix")
    } else {
      matrix(0, ncol = 1)
    },
    quad_size = po.data[pt.quad.id == 0, quad.weights.name],
    Y = 0,
    bf_per_res = as.numeric(table(bf.info$res)),
    mod_type = as.integer(which(model.type == c("ipp", "variational", "laplace")) - 1),
    data_type = 0,
    bf_matrix_type = bf.matrix.type
  )

  # create the appropraite start parameters for the variance component w.r.t. approx. type
  var.starts <- switch(model.type,
                       ipp = rep(0, ncol(dat.list$Z_PO_pres)),
                       variational = rep(0, ncol(dat.list$Z_PO_pres)),
                       laplace = rep(0, length(dat.list$bf_per_res))
  )

  # initilise starting parameters at 0
  start.pars <- list(fixed = rep(0, ncol(dat.list$X_PO_pres)),
                     bias = rep(0, ncol(dat.list$B_PO_pres)),
                     random = rep(0, ncol(dat.list$Z_PO_pres)),
                     log_variance_component = var.starts
  )

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

  # set up the objective function w.r.t. model.type
  obj <- switch(model.type,
                ipp = TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", map = list(bias = as.factor(rep(NA, ncol(dat.list$B_PO_pres))), random = factor(rep(NA, ncol(dat.list$Z_PA))), log_variance_component = factor(rep(NA, length(var.starts)))), silent = T),
                variational = TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", map = list(bias = as.factor(rep(NA, ncol(dat.list$B_PO_pres)))), silent = T),
                laplace = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", map = list(bias = as.factor(rep(NA, ncol(dat.list$B_PO_pres)))), silent = T)
  )
  # optimise the parameters
  res <- stats::optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS", control = list(maxit = 1000))
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
  coef.names <- switch(model.type,
                       ipp = fixed.names,
                       variational = c(fixed.names, paste0("VA Posterior Mean (bf ", random.nos, ")"), paste0("VA Posterior log Var (bf ", random.nos, ")")),
                       laplace = c(fixed.names, paste0("Prior log sd (res. ", 1:length(dat.list$bf_per_res), ")"))
  )
  names(res$coefficients) <- coef.names
  res$fixed.effects <- tmp.estimates[1:length(fixed.names), ]
  rownames(res$fixed.effects) <- fixed.names
  if (model.type != "ipp") {
    res$random.effects <- tmp.estimates[(length(fixed.names) + 1):nrow(tmp.estimates), ]
    res$random.effects <- res$random.effects[!grepl("log_", rownames(res$random.effects), fixed = T), ]
    rand.names <- switch(model.type,
                         variational = c(paste0("VA Posterior Mean (bf ", random.nos, ")"), paste0("VA Posterior Var (bf ", random.nos, ")"), paste0("Prior Var (res. ", 1:length(dat.list$bf_per_res), ")")),
                         laplace = c(paste0("LP Posterior Mean (bf ", random.nos, ")"), paste0("Prior Var (res. ", 1:length(dat.list$bf_per_res), ")"))
    )
    rownames(res$random.effects) <- rand.names
    res$basis.per.res <- dat.list$bf_per_res
    res$FRK.basis.functions <- FRK.basis.functions
    res$basis.fn.info <- bf.info
    res$approx.type <- model.type
    res$fitted.values <- as.vector(po.des.mat %*% res$fixed.effects[ , 1] + po.bf.matrix %*% res$random.effects[grepl(" Mean ", rownames(res$random.effects), fixed = T), 1])
  } else {
    res$random.effects <- NA
    res$basis.per.res <- NA
    res$FRK.basis.functions <- NULL
    res$basis.fn.info <- NULL
    res$approx.type <- NA
    res$fitted.values <- as.vector(po.des.mat %*% res$fixed.effects[ , 1])
  }
  res$starting.pars <- start.pars
  res$data <- po.data
  res$formula <- po.formula
  res$coord.names <- coord.names
  res$quad.weights.name <- quad.weights.name
  res$pt.quad.id <- pt.quad.id
  res$data.model.type <- "po"
  res$bf.matrix.type <- bf.matrix.type
  class(res) <- "scampr"
  return(res)
}
