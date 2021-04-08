#' Internal scampr function that creates a list of data and starting parameters for scampr models
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the Presence-only data model to be fitted. The 'response' must be a binary that indicates whether a datum is a presence or quadrature point. See GLM function for further formula details.
#' @param pa.formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the Presence/Absence data model to be fitted. The 'response' must be a must be either the site abundance or a binary indicating whether there is a presence or not. See GLM function for further formula details.
#' @param data a data frame containing predictors at both presence-records and quadrature as well as the formula 'response'.
#' @param pa.data a data frame containing predictors and response for the pa.formula.
#' @param coord.names a vector of character strings describing the column names of the coordinates in both data frames.
#' @param quad.weights.name a charater string of the column name of quadrature weights in the data.
#' @param FRK.basis.functions an optional object of class 'Basis' from FRK package. If neither 'FRK.basis.functions' nor 'simple.basis' is specified will use default FRK::auto_basis with 2 spatial resolutions.
#' @param simple.basis an alternative to 'FRK.basis.functions': a data.frame of basis functions information created by 'simple_basis()'.
#' @param model.type a character string indicating the type of model to be used. May be one of 'laplace' or 'variational' for Cox Processes involving spatially correlated errors or 'ipp' for a model that follows an inhomgeneous Poisson process.
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to use sparse or dense matrix computations for the basis functions created.
#' @param data.type a logical indicating whether standard errors should be calculated.
#' @param starting.pars an optional named list or scampr model object that gives warm starting values for the parameters of the model.
#'
#' @return list of elements required for TMB::MakeADFun
#' @noRd
#'
#' @importFrom methods as
#' @importFrom stats as.formula
#' @importFrom sp coordinates
#' @importFrom FRK auto_basis eval_basis
#'
#' @examples
#' # Get the flora data for one of the species
#' dat_po <- flora$po$sp1
#' dat_pa <- flora$pa
#'
#' # Attach the quadrature to the PO data
#' dat_po <- rbind.data.frame(dat_po, flora$quad)
#'
#' # Get the TMB data lists for a combined data model without latent field
#' tmb.input <- scampr:::get.TMB.data.input(pres ~ MNT + D.Main, sp1 ~ MNT, po.data = dat_po, pa.data = dat_pa, model.type = "ipp")
#' str(tmp.input)
get.TMB.data.input <- function(formula, pa.formula, data, pa.data, coord.names = c("x", "y"), quad.weights.name = "quad.size", FRK.basis.functions, simple.basis, model.type = c("laplace", "variational", "ipp"), bf.matrix.type = c("sparse", "dense"), data.type = c("po", "pa", "popa"), starting.pars) {

  # parameters of restricted strings
  data.type <- match.arg(data.type)
  model.type <- match.arg(model.type)
  bf.matrix.type <- match.arg(bf.matrix.type)

  if (data.type == "popa") {

    pa.resp <- all.vars(pa.formula[[2]])
    pa.pred <- all.vars(pa.formula[[3]])
    po.resp <- all.vars(formula[[2]])
    po.pred <- all.vars(formula[[3]])

    ## checks ##

    if (length(pa.resp) != 1 & length(po.resp) != 1) {
      stop("Both formulae can only take a single response")
    }
    if (!all(c(pa.resp, pa.pred) %in% colnames(pa.data))) {
      stop("PA data does not contain the formula terms")
    }
    if (!all(c(po.resp, po.pred) %in% colnames(data))) {
      stop("PO data does not contain the formula terms")
    }
    if (!all(coord.names %in% colnames(pa.data)) & !all(coord.names %in% colnames(data))) {
      stop(paste0("coord.names, ", coord.names, ", not found in both data sets provided"))
    }
    if (!quad.weights.name %in% colnames(data)) {
      stop(paste0("quad.weights.name, ", quad.weights.name, ", not found in the PO data provided"))
    }
    if (!all(pa.pred %in% po.pred)) {
      stop("Not all predictors in the PA formula are found in the PO formula")
    }
    if (all(po.pred %in% pa.pred)) {
      warning("Predictors are the same for PO and PA.\nUnless the presence records are the result of a perfectly observed domain, this is misspecified")
    }
    if (model.type == "variational") {
      stop("Combined Data Model cannot handle variational approx. at this stage. Try model.type = 'laplace'")
    }
    ############################################################

    # Create a vector of row names to track changes data manipulations (for po data)
    row.id <- 1:nrow(data)

    # Get the PA design matrix
    pa.des.mat <- get.design.matrix(pa.formula, pa.data)
    # Get the full PO design matrix
    full.po.des.mat <- get.design.matrix(formula, data)
    # Get the seperate predictors (and remove intercepts if present - to be dealt with later)
    po.preds <- colnames(full.po.des.mat)[colnames(full.po.des.mat) %in% colnames(pa.des.mat)]
    po.preds <- po.preds[!grepl("Intercept", po.preds, fixed = T)]
    bias.preds <- colnames(full.po.des.mat)[!colnames(full.po.des.mat) %in% colnames(pa.des.mat)]
    bias.preds <- bias.preds[!grepl("Intercept", bias.preds, fixed = T)]
    # Add/remove intercept term as required
    if (any(grepl("Intercept", colnames(full.po.des.mat), fixed = T))) {
      if (!any(grepl("Intercept", colnames(pa.des.mat), fixed = T))) {
        stop("Intercept term found in formula and not in pa.formula")
      }
      po.pred.formula <- stats::as.formula(paste0(po.resp, " ~ ", paste(po.preds, collapse = " + ")))
      po.bias.formula <- stats::as.formula(paste0(po.resp, " ~ ", paste(bias.preds, collapse = " + ")))
    } else {
      if (any(grepl("Intercept", colnames(pa.des.mat), fixed = T))) {
        stop("Intercept term found in pa.formula and not in formula")
      }
      po.pred.formula <- stats::as.formula(paste0(po.resp, " ~ - 1 + ", paste(po.preds, collapse = " + ")))
      po.bias.formula <- stats::as.formula(paste0(po.resp, " ~ - 1 + ", paste(bias.preds, collapse = " + ")))
    }
    # Get the PO design matrices
    po.des.mat <- get.design.matrix(po.pred.formula, data)
    po.bias.des.mat <- get.design.matrix(po.bias.formula, data)
    # Adjust the Intercept name if required
    if (any(grepl("Intercept", colnames(po.bias.des.mat), fixed = T))) {
      colnames(po.bias.des.mat)[grepl("Intercept", colnames(po.bias.des.mat), fixed = T)] <- "(Bias Intercept)"
    }
    # Get the presence point/ quadrature point identifier
    pt.quad.id <- data[ , po.resp]
    # Set the fixed effect names
    fixed.names <- c(colnames(po.des.mat), colnames(po.bias.des.mat))
    fixed.names.bias.id <- c(rep(F, ncol(po.des.mat)), rep(T, ncol(po.bias.des.mat)))

    # Determine the basis functions to be used
    if (model.type != "ipp") {
      if (missing(simple.basis)) {
        if (missing(FRK.basis.functions)) { # When none is provided use FRK defaults (with 2 spatial resolutions)
          # create a spatial pixels data frame as required by FRK::auto_basis
          sp.data <- rbind(data[ ,coord.names], pa.data[ , coord.names])
          sp::coordinates(sp.data) <- coord.names
          FRK.basis.functions <- FRK::auto_basis(data = sp.data, max_basis = sum(pt.quad.id)*0.25)
        }
        po.bf.matrix <- FRK::eval_basis(basis = FRK.basis.functions, as.matrix(data[ , coord.names]))
        pa.bf.matrix <- FRK::eval_basis(basis = FRK.basis.functions, as.matrix(pa.data[ , coord.names]))
        bf.info <- FRK.basis.functions@df
        colnames(bf.info)[grepl("loc", colnames(bf.info), fixed = T)] <- coord.names
        class(bf.info) <- c(class(bf.info), "bf.df")
      } else { # Otherwise use the provided simple basis
        po.bf.matrix <- get.bf.matrix(simple.basis, data[ , coord.names], bf.matrix.type = bf.matrix.type)
        pa.bf.matrix <- get.bf.matrix(simple.basis, pa.data[ , coord.names], bf.matrix.type = bf.matrix.type)
        bf.info <- simple.basis
        FRK.basis.functions <- NULL
      }
    } else {
      po.bf.matrix <- matrix(rep(0, nrow(po.des.mat)), ncol = 1)
      pa.bf.matrix <- matrix(0, nrow = 1)
      bf.info <- cbind.data.frame(x = NA, y = NA, scale = NA, res = 1)
      FRK.basis.functions <- NULL
    }

    # TMB required data setup
    dat.list <- list(
      X_PO_pres = as.matrix(po.des.mat[pt.quad.id == 1, ]),
      B_PO_pres = as.matrix(po.bias.des.mat[pt.quad.id == 1, ]),
      X_PO_quad = as.matrix(po.des.mat[pt.quad.id == 0, ]),
      B_PO_quad = as.matrix(po.bias.des.mat[pt.quad.id == 0, ]),
      X_PA = as.matrix(pa.des.mat),
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
        methods::as(pa.bf.matrix, "sparseMatrix")
      } else {
        as.matrix(pa.bf.matrix)
      },
      quad_size = data[pt.quad.id == 0, quad.weights.name],
      Y = pa.data[ , pa.resp] > 0, # corrects in the case of abundance
      bf_per_res = as.numeric(table(bf.info$res)),
      mod_type = as.integer(which(model.type == c("ipp", "variational", "laplace")) - 1),
      data_type = 2,
      bf_matrix_type = bf.matrix.type
    )
    # track row number movements
    pres.rows <- row.id[pt.quad.id == 1]
    quad.rows <- row.id[pt.quad.id == 0]

    # Starting Parameters #

    # if starting.pars provided is a scampr model adjust to req. list structure
    if (!missing(starting.pars)) {
      if (class(starting.pars) == "scampr") {
        tmp.m <- starting.pars
        starting.pars <- lapply(split(starting.pars$par, names(starting.pars$par)), unname)
        # check the model isn't an IPP
        if (!is.na(tmp.m$approx.type)) {
          if (model.type == "laplace" & tmp.m$approx.type == "variational") {
            # make appropriate change to the variance parameter if going from VA to Laplace
            starting.pars$log_variance_component <- unname(log(sqrt(tmp.m$random.effects[grepl("Prior Var ", rownames(tmp.m$random.effects), fixed = T), 1])))
          } else if (model.type == "variational" & tmp.m$approx.type == "laplace") {
            # make appropriate change to the variance parameter going from Laplace to VA
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
    # create the appropraite start parameters for the variance component w.r.t. approx. type
    var.starts <- rep(0, length(dat.list$bf_per_res))
    # initilise starting parameters at 0
    start.pars <- list(fixed = rep(0, ncol(dat.list$X_PA)),
                       bias = rep(0, ncol(dat.list$B_PO_pres)),
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

    # Collate info to be returned
    return.info <- list(tmb.data = dat.list, tmb.pars = start.pars, po.info = cbind(pt.quad.id, row.id = order(c(pres.rows, quad.rows))), fixed.names = fixed.names, fixed.names.bias.id = fixed.names.bias.id, bf.info = bf.info, FRK.basis.functions = FRK.basis.functions)

    ###################################################################
  } else if (data.type == "po") {

    po.resp <- all.vars(formula[[2]])
    po.pred <- all.vars(formula[[3]])

    ## checks ##

    if (length(po.resp) != 1) {
      stop("Formula can only take a single response")
    }
    if (!all(c(po.resp, po.pred) %in% colnames(data))) {
      stop("PO data does not contain the formula terms")
    }
    if (!all(coord.names %in% colnames(data))) {
      stop(paste0("coord.names, ", coord.names, ", not found in data provided"))
    }
    if (!quad.weights.name %in% colnames(data)) {
      stop(paste0("quad.weights.name, ", quad.weights.name, ", not found in data provided"))
    }
    ############################################################

    # Create a vector of row names to track changes data manipulations (for po data)
    row.id <- 1:nrow(data)

    # Get the PO design matrix
    po.des.mat <- get.design.matrix(formula, data)
    # Get the presence point/ quadrature point identifier
    pt.quad.id <- data[ , po.resp]
    # Set the fixed effect names
    fixed.names <- colnames(po.des.mat)

    # Determine the basis functions to be used
    if (model.type != "ipp") {
      if (missing(simple.basis)) {
        if (missing(FRK.basis.functions)) { # When none is provided use FRK defaults (with 2 spatial resolutions)
          # create a spatial pixels data frame as required by FRK::auto_basis
          sp.data <- data[ ,coord.names]
          sp::coordinates(sp.data) <- coord.names
          FRK.basis.functions <- FRK::auto_basis(data = sp.data, max_basis = sum(pt.quad.id)*0.25)
        }
        po.bf.matrix <- FRK::eval_basis(basis = FRK.basis.functions, as.matrix(data[ , coord.names]))
        bf.info <- FRK.basis.functions@df
        colnames(bf.info)[grepl("loc", colnames(bf.info), fixed = T)] <- coord.names
        class(bf.info) <- c(class(bf.info), "bf.df")
      } else { # Otherwise use the provided simple basis
        po.bf.matrix <- get.bf.matrix(simple.basis, data[ , coord.names], bf.matrix.type = bf.matrix.type)
        bf.info <- simple.basis
        FRK.basis.functions <- NULL
      }
    } else {
      po.bf.matrix <- matrix(rep(0, nrow(po.des.mat)), ncol = 1)
      bf.info <- cbind.data.frame(x = NA, y = NA, scale = NA, res = 1)
      FRK.basis.functions <- NULL
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
      quad_size = data[pt.quad.id == 0, quad.weights.name],
      Y = 0,
      bf_per_res = as.numeric(table(bf.info$res)),
      mod_type = as.integer(which(model.type == c("ipp", "variational", "laplace")) - 1),
      data_type = 0,
      bf_matrix_type = bf.matrix.type
    )
    # track row number movements
    pres.rows <- row.id[pt.quad.id == 1]
    quad.rows <- row.id[pt.quad.id == 0]

    # Starting Parameters #
    # if starting.pars provided is a scampr model adjust to req. list structure
    if (!missing(starting.pars)) {
      if (class(starting.pars) == "scampr") {
        tmp.m <- starting.pars
        starting.pars <- lapply(split(starting.pars$par, names(starting.pars$par)), unname)
        # check the model isn't an IPP
        if (!is.na(tmp.m$approx.type)) {
          if (model.type == "laplace" & tmp.m$approx.type == "variational") {
            # make appropriate change to the variance parameter if going from VA to Laplace
            starting.pars$log_variance_component <- unname(log(sqrt(tmp.m$random.effects[grepl("Prior Var ", rownames(tmp.m$random.effects), fixed = T), 1])))
          } else if (model.type == "variational" & tmp.m$approx.type == "laplace") {
            # make appropriate change to the variance parameter going from Laplace to VA
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

    # Collate info to be returned
    return.info <- list(tmb.data = dat.list, tmb.pars = start.pars, po.info = cbind(pt.quad.id, row.id = order(c(pres.rows, quad.rows))), fixed.names = fixed.names, fixed.names.bias.id = NULL, bf.info = bf.info, FRK.basis.functions = FRK.basis.functions)


    ###################################################################
  } else if (data.type == "pa") {

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
      stop(paste0("coord.names, ", coord.names, ", not found in data provided"))
    }
    if (model.type == "variational") {
      stop("PA Data Model cannot handle variational approx. at this stage. Try model.type = 'laplace'")
    }
    ############################################################

    # Get the PA design matrix
    pa.des.mat <- get.design.matrix(pa.formula, pa.data)
    fixed.names <- colnames(pa.des.mat)

    # Determine the basis functions to be used
    if (model.type != "ipp") {
      if (missing(simple.basis)) {
        if (missing(FRK.basis.functions)) { # When none is provided use FRK defaults (with 2 spatial resolutions)
          # create a spatial pixels data frame as required by FRK::auto_basis
          sp.data <- pa.data[ , coord.names]
          sp::coordinates(sp.data) <- coord.names
          FRK.basis.functions <- FRK::auto_basis(data = sp.data, max_basis = sum(pt.quad.id)*0.25)
        }
        pa.bf.matrix <- FRK::eval_basis(basis = FRK.basis.functions, as.matrix(pa.data[ , coord.names]))
        bf.info <- FRK.basis.functions@df
        colnames(bf.info)[grepl("loc", colnames(bf.info), fixed = T)] <- coord.names
        class(bf.info) <- c(class(bf.info), "bf.df")
      } else { # Otherwise use the provided simple basis
        pa.bf.matrix <- get.bf.matrix(simple.basis, pa.data[ , coord.names], bf.matrix.type = bf.matrix.type)
        bf.info <- simple.basis
        FRK.basis.functions <- NULL
      }
    } else {
      pa.bf.matrix <- matrix(0, nrow = 1)
      bf.info <- cbind.data.frame(x = NA, y = NA, scale = NA, res = 1)
      FRK.basis.functions <- NULL
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

    # Starting Parameters #

    # if starting.pars provided is a scampr model adjust to req. list structure
    if (!missing(starting.pars)) {
      if (class(starting.pars) == "scampr") {
        tmp.m <- starting.pars
        starting.pars <- lapply(split(starting.pars$par, names(starting.pars$par)), unname)
        # check the model isn't an IPP
        if (!is.na(tmp.m$approx.type)) {
          if (model.type == "laplace" & tmp.m$approx.type == "variational") {
            # make appropriate change to the variance parameter if going from VA to Laplace
            starting.pars$log_variance_component <- unname(log(sqrt(tmp.m$random.effects[grepl("Prior Var ", rownames(tmp.m$random.effects), fixed = T), 1])))
          } else if (model.type == "variational" & tmp.m$approx.type == "laplace") {
            # make appropriate change to the variance parameter going from Laplace to VA
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
    # create the appropraite start parameters for the variance component w.r.t. approx. type
    var.starts <- rep(0, length(dat.list$bf_per_res))
    # initilise starting parameters at 0
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

    # Collate info to be returned
    return.info <- list(tmb.data = dat.list, tmb.pars = start.pars, po.info = NULL, fixed.names = fixed.names, fixed.names.bias.id = NULL, bf.info = bf.info, FRK.basis.functions = FRK.basis.functions)

  } else {
    stop("Incorrect 'data.type' provided. Must be one of 'po', 'pa' or 'popa'")
  }
  return(return.info)
}
