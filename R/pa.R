#' Model for Presence/Absence Data
#'
#' @description Fits a binary regression model to presence/absence data using a complimentary log-log link function. Can accomodate an approx. latent field as spatial random effects.
#'
#' @param pa.formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the Presence/Absence data model to be fitted. The 'response' must be a must be either the site abundance or a binary indicating whether there is a presence or not. See GLM function for further formula details.
#' @param pa.data a data frame containing predictors and response for the pa.formula.
#' @param coord.names a vector of character strings describing the column names of the coordinates in both data frames
#' @param FRK.basis.functions an optional object of class 'Basis' from FRK package. If neither 'FRK.basis.functions' nor 'simple.basis' is specified will use default FRK::auto_basis with 2 spatial resolutions.
#' @param simple.basis an alternative to 'FRK.basis.functions': a data.frame of basis functions information created by 'simple_basis()'.
#' @param model.type a character string indicating the type of model to be used. May be one of 'laplace' or 'variational' for Cox Processes involving spatially correlated errors or 'ipp' for a model that follows an inhomgeneous Poisson process.
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to use sparse or dense matrix computations for the basis functions created.
#' @param se a logical indicating whether standard errors should be calculated.
#' @param starting.pars an optional named list or scampr model object that gives warm starting values for the parameters of the model.
#' @param subset an optional subset of the data to be used.
#' @param na.action an optional way of handling NA's in the data, default is omit.
#'
#' @return a scampr model object
#' @export
#'
#' @importFrom methods as
#' @importFrom stats optim
#' @importFrom sp coordinates
#' @importFrom FRK auto_basis eval_basis
#' @importFrom TMB MakeADFun sdreport
#'
#' @examples
#' # Fit without a shared latent field
#' m1 <- pa(Y ~ TMP_MIN, pa.data = eucalypt[["pa"]], model.type = "ipp")
#'
#' # Fit with a shared latent field using FRK package default basis functions
#' m2 <- pa(Y ~ TMP_MIN, eucalypt[["pa"]], model.type = "laplace")
pa <- function(pa.formula, pa.data, coord.names = c("x", "y"), FRK.basis.functions, simple.basis, model.type = c("laplace", "variational", "ipp"), bf.matrix.type = c("sparse", "dense"), se = TRUE, starting.pars, subset, na.action) {

  # CAN'T JUST GIVE A BASIS FUNCTION MATRIX BECAUSE THEN YOU CAN'T PREDICT ETC. AS WE DON'T KNOW ENOUGH ABOUT THE FUNCTIONS

  # Get the response and predictor names
  pa.resp <- all.vars(pa.formula[[2]])
  pa.pred <- all.vars(pa.formula[[3]])

  ## checks ##

  if (length(pa.resp) != 1) {
    stop("Formula can only take a single response")
  }
  if (!all(c(pa.resp, pa.pred) %in% colnames(pa.data))) {
    stop("PA data does not contain the formula terms")
  }
  if (!all(coord.names %in% colnames(pa.data))) {
    stop(paste0("coord.names, ", coord.names, ", not found in data set provided"))
  }
  if (!is.logical(se)) {
    stop(paste0("'se' must be a logcial indicating whether or not to calculate standard errors"))
  }
  # parameters of restricted strings
  model.type <- match.arg(model.type)
  bf.matrix.type <- match.arg(bf.matrix.type)

  ############################################################

  # Get the PA design matrix
  pa.des.mat <- get.desgin.matrix(pa.formula, pa.data)
  fixed.names <- colnames(pa.des.mat)

  # Determine the basis functions to be used
  if (model.type != "ipp") {
    if (missing(simple.basis)) {
      if (missing(FRK.basis.functions)) { # When none is provided use FRK defaults (with 2 spatial resolutions)
        # create a spatial pixels data frame as required by FRK::auto_basis
        sp.data <- pa.data[ , coord.names]
        sp::coordinates(sp.data) <- coord.names
        FRK.basis.functions <- FRK::auto_basis(data = sp.data, nres = 2)
      }
      pa.bf.matrix <- FRK::eval_basis(basis = FRK.basis.functions, as.matrix(pa.data[ , coord.names]))
      bf.info <- FRK.basis.functions@df
      colnames(bf.info)[grepl("loc", colnames(bf.info), fixed = T)] <- coord.names
      class(bf.info) <- c(class(bf.info), "bf.df")
    } else { # Otherwise use the provided simple basis
      pa.bf.matrix <- get.bf.matrix(simple.basis, pa.data[ , coord.names])
      bf.info <- simple.basis
      FRK.basis.functions <- NULL
    }
  } else {
    pa.bf.matrix <- matrix(0, nrow = 1)
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
    X_PO_pres = matrix(0, ncol = 1),
    B_PO_pres = matrix(0, ncol = 1),
    X_PO_quad = matrix(0, ncol = 1),
    B_PO_quad = matrix(0, ncol = 1),
    X_PA = as.matrix(pa.des.mat),
    Z_PO_pres = if(bf.matrix.type == "sparse") {
      methods::as(matrix(0, ncol = 1), "sparseMatrix")
    } else {
      matrix(0, ncol = 1)
    },
    Z_PO_quad = if(bf.matrix.type == "sparse") {
      methods::as(matrix(0, ncol = 1), "sparseMatrix")
    } else {
      matrix(0, ncol = 1)
    },
    Z_PA = if(bf.matrix.type == "sparse") {
      methods::as(pa.bf.matrix, "sparseMatrix")
    } else {
      as.matrix(pa.bf.matrix)
    },
    quad_size = 0,
    Y = pa.data[ , pa.resp] > 0, # corrects in the case of abundance
    bf_per_res = as.numeric(table(bf.info$res)),
    mod_type = as.integer(which(model.type == c("ipp", "variational", "laplace")) - 1),
    data_type = 1,
    bf_matrix_type = bf.matrix.type
  )

  # # AT THIS STAGE CAN ALONE PERFORM LAPLAC APPROX.
  # # create the appropraite start parameters for the variance component
  # #   w.r.t. approx. type and data type
  # var.starts <- switch(model.type,
  #                      variational = rep(0, ncol(dat.list$Z_PO_pres)),
  #                      laplace = rep(0, length(dat.list$bf_per_res))
  # )
  var.starts <- rep(0, length(dat.list$bf_per_res))
  start.pars <- list(fixed = rep(0, ncol(dat.list$X_PA)),
                     bias = 0,
                     random = rep(0, ncol(dat.list$Z_PA)),
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

  # # AT THIS STAGE CAN ONLY PERFORM LAPLAC APPROX.
  # set up the objective function w.r.t. model.type
  obj <- switch(model.type,
                ipp = TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", map = list(bias = as.factor(rep(NA, ncol(dat.list$B_PO_pres))), random = factor(rep(NA, ncol(dat.list$Z_PA))), log_variance_component = factor(rep(NA, length(var.starts)))), silent = T),
                variational = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", map = list(bias = as.factor(rep(NA, ncol(dat.list$B_PO_pres)))), silent = T),
                laplace = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", map = list(bias = as.factor(rep(NA, ncol(dat.list$B_PO_pres)))), silent = T)
  )
  # obj <- TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", silent = T)
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
                       variational = c(fixed.names, paste0("VA Posterior Mean (bf ", random.nos, ")"), paste0("VA Posterior log sd (bf ", random.nos, ")")),
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
    res$approx.type <- "laplace" #model.type until PA models can handle VA
    res$fitted.values <- as.vector(pa.des.mat %*% res$fixed.effects[ , 1] + pa.bf.matrix %*% res$random.effects[grepl(" Mean ", rownames(res$random.effects), fixed = T), 1])
  } else {
    res$random.effects <- NA
    res$basis.per.res <- NA
    res$FRK.basis.functions <- NULL
    res$basis.fn.info <- NULL
    res$approx.type <- NA
    res$fitted.values <- as.vector(pa.des.mat %*% res$fixed.effects[ , 1])
  }
  res$starting.pars <- start.pars
  res$data <- pa.data
  res$formula <- pa.formula
  res$coord.names <- coord.names
  res$quad.weights.name <- NA
  res$pt.quad.id <- NA
  res$data.model.type <- "pa"
  res$bf.matrix.type <- bf.matrix.type
  class(res) <- "scampr"
  return(res)
}
