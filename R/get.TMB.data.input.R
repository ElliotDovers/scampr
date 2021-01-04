#' Internal scampr function that creates a list of data and starting parameters for scampr models
#'
#' @param po.formula
#' @param pa.formula
#' @param po.data
#' @param pa.data
#' @param coord.names
#' @param quad.weights.name
#' @param FRK.basis.functions
#' @param simple.basis
#' @param model.type
#' @param bf.matrix.type
#' @param data.type
#' @param starting.pars
#'
#' @return
#'
#' @examples
get.TMB.data.input <- function(po.formula, pa.formula, po.data, pa.data, coord.names = c("x", "y"), quad.weights.name = "quad.size", FRK.basis.functions, simple.basis, model.type = c("laplace", "variational", "ipp"), bf.matrix.type = c("sparse", "dense"), data.type = c("po", "pa", "popa"), starting.pars) {

  # parameters of restricted strings
  data.type <- match.arg(data.type)
  model.type <- match.arg(model.type)
  bf.matrix.type <- match.arg(bf.matrix.type)

  if (data.type == "popa") {

    pa.resp <- all.vars(pa.formula[[2]])
    pa.pred <- all.vars(pa.formula[[3]])
    po.resp <- all.vars(po.formula[[2]])
    po.pred <- all.vars(po.formula[[3]])

    ## checks ##

    if (length(pa.resp) != 1 & length(po.resp) != 1) {
      stop("Both formulae can only take a single response")
    }
    if (!all(c(pa.resp, pa.pred) %in% colnames(pa.data))) {
      stop("PA data does not contain the formula terms")
    }
    if (!all(c(po.resp, po.pred) %in% colnames(po.data))) {
      stop("PO data does not contain the formula terms")
    }
    if (!all(coord.names %in% colnames(pa.data)) & !all(coord.names %in% colnames(po.data))) {
      stop(paste0("coord.names, ", coord.names, ", not found in both data sets provided"))
    }
    if (!quad.weights.name %in% colnames(po.data)) {
      stop(paste0("quad.weights.name, ", quad.weights.name, ", not found in the PO data provided"))
    }
    if (!all(pa.pred %in% po.pred)) {
      stop("Not all predictors in the PA formula are found in the PO formula")
    }
    if (all(po.pred %in% pa.pred)) {
      warning("Predictors are the same for PO and PA.\nUnless the presence records are the result of a perfectly observed domain, this is misspecified")
    }
    ############################################################

    # Get the PA design matrix
    pa.des.mat <- scampr:::get.desgin.matrix(pa.formula, pa.data)
    # Determine the bias predictors as those in PO formula and not in PA formula
    bias.preds <- po.pred[!po.pred %in% pa.pred]
    # Separate the PO formulae
    po.pred.formula <- as.formula(paste0(po.resp, " ~ ", paste(po.pred[!po.pred %in% bias.preds], collapse = " + ")))
    po.bias.formula <- as.formula(paste0(po.resp, " ~ ", paste(bias.preds, collapse = " + ")))
    # Get the PO design matrices
    po.des.mat <- scampr:::get.desgin.matrix(po.pred.formula, po.data)
    po.bias.des.mat <- scampr:::get.desgin.matrix(po.bias.formula, po.data)
    # Get the presence point/ quadrature point identifier
    pt.quad.id <- po.data[ , po.resp]
    # Re-adjust the bias intercept name
    colnames(po.bias.des.mat)[grepl("(Intercept)", colnames(po.bias.des.mat), fixed = T)] <- "(Bias Intercept)"
    # Set the fixed effect names
    fixed.names <- c(colnames(po.des.mat), colnames(po.bias.des.mat))
    fixed.names.bias.id <- c(rep(F, ncol(po.des.mat)), rep(T, ncol(po.bias.des.mat)))

    # Determine the basis functions to be used
    if (model.type != "ipp") {
      if (missing(simple.basis)) {
        if (missing(FRK.basis.functions)) { # When none is provided use FRK defaults (with 2 spatial resolutions)
          # create a spatial pixels data frame as required by FRK::auto_basis
          sp.data <- rbind(po.data[ ,coord.names], pa.data[ , coord.names])
          sp::coordinates(sp.data) <- coord.names
          FRK.basis.functions <- FRK::auto_basis(data = sp.data, nres = 2)
        }
        po.bf.matrix <- FRK::eval_basis(basis = FRK.basis.functions, as.matrix(po.data[ , coord.names]))
        pa.bf.matrix <- FRK::eval_basis(basis = FRK.basis.functions, as.matrix(pa.data[ , coord.names]))
        bf.info <- FRK.basis.functions@df
        colnames(bf.info)[grepl("loc", colnames(bf.info), fixed = T)] <- coord.names
        class(bf.info) <- c(class(bf.info), "bf.df")
      } else { # Otherwise use the provided simple basis
        po.bf.matrix <- scampr:::get.bf.matrix(simple.basis, po.data[ , coord.names])
        pa.bf.matrix <- scampr:::get.bf.matrix(simple.basis, pa.data[ , coord.names])
        bf.info <- simple.basis
        FRK.basis.functions <- NULL
      }
    } else {
      po.bf.matrix <- matrix(rep(0, nrow(po.des.mat)), ncol = 1)
      pa.bf.matrix <- matrix(0, nrow = 1)
      bf.info <- cbind.data.frame(x = NA, y = NA, scale = NA, res = 1)
    }

    # TMB required data setup
    dat.list <- list(
      X_PO_pres = as.matrix(po.des.mat[pt.quad.id == 1, ]),
      B_PO_pres = as.matrix(po.bias.des.mat[pt.quad.id == 1, ]),
      X_PO_quad = as.matrix(po.des.mat[pt.quad.id == 0, ]),
      B_PO_quad = as.matrix(po.bias.des.mat[pt.quad.id == 0, ]),
      X_PA = as.matrix(pa.des.mat),
      Z_PO_pres = if(bf.matrix.type == "sparse") {
        as(po.bf.matrix[pt.quad.id == 1, ], "sparseMatrix")
      } else {
        as.matrix(po.bf.matrix[pt.quad.id == 1, ])
      },
      Z_PO_quad = if(bf.matrix.type == "sparse") {
        as(po.bf.matrix[pt.quad.id == 0, ], "sparseMatrix")
      } else {
        as.matrix(po.bf.matrix[pt.quad.id == 0, ])
      },
      Z_PA = if(bf.matrix.type == "sparse") {
        as(pa.bf.matrix, "sparseMatrix")
      } else {
        as.matrix(pa.bf.matrix)
      },
      quad_size = po.data[pt.quad.id == 0, quad.weights.name],
      Y = pa.data[ , pa.resp] > 0, # corrects in the case of abundance
      bf_per_res = as.numeric(table(bf.info$res)),
      mod_type = as.integer(which(model.type == c("ipp", "variational", "laplace")) - 1),
      data_type = 2,
      bf_matrix_type = bf.matrix.type
    )

    # Starting Parameters #

    # if warm starts are provided
    if (!missing(starting.pars)) {
      if (class(starting.pars) == "scampr") {
        mod.obj <- starting.pars
        starting.pars <- lapply(split(mod.obj$par, names(mod.obj$par)), unname)
        if (model.type == "laplace") {
          starting.pars$random <- unname(mod.obj$random.effects[grepl("LP Posterior Mean", rownames(mod.obj$random.effects), fixed = T), 1L])
        }
      }
    }
    var.starts <- rep(0, length(dat.list$bf_per_res))
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
    ###################################################################
  } else if (data.type == "po") {

    po.resp <- all.vars(po.formula[[2]])
    po.pred <- all.vars(po.formula[[3]])

    ## checks ##

    if (length(po.resp) != 1) {
      stop("Formula can only take a single response")
    }
    if (!all(c(po.resp, po.pred) %in% colnames(po.data))) {
      stop("PO data does not contain the formula terms")
    }
    if (!all(coord.names %in% colnames(po.data))) {
      stop(paste0("coord.names, ", coord.names, ", not found in both data sets provided"))
    }
    if (!quad.weights.name %in% colnames(po.data)) {
      stop(paste0("quad.weights.name, ", quad.weights.name, ", not found in the PO data provided"))
    }
    ############################################################

    # Get the PO design matrix
    po.des.mat <- scampr:::get.desgin.matrix(po.formula, po.data)
    # Get the presence point/ quadrature point identifier
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
          FRK.basis.functions <- FRK::auto_basis(data = sp.data, nres = 2)
        }
        po.bf.matrix <- FRK::eval_basis(basis = FRK.basis.functions, as.matrix(po.data[ , coord.names]))
        bf.info <- FRK.basis.functions@df
        colnames(bf.info)[grepl("loc", colnames(bf.info), fixed = T)] <- coord.names
        class(bf.info) <- c(class(bf.info), "bf.df")
      } else { # Otherwise use the provided simple basis
        po.bf.matrix <- scampr:::get.bf.matrix(simple.basis, po.data[ , coord.names])
        bf.info <- simple.basis
        FRK.basis.functions <- NULL
      }
    } else {
      po.bf.matrix <- matrix(rep(0, nrow(po.des.mat)), ncol = 1)
      bf.info <- cbind.data.frame(x = NA, y = NA, scale = NA, res = 1)
    }
    # TMB required data setup
    dat.list <- list(
      X_PO_pres = as.matrix(po.des.mat[pt.quad.id == 1, ]),
      B_PO_pres = matrix(rep(0, sum(pt.quad.id == 1)), ncol = 1),
      X_PO_quad = as.matrix(po.des.mat[pt.quad.id == 0, ]),
      B_PO_quad = matrix(rep(0, sum(pt.quad.id == 0)), ncol = 1),
      X_PA = matrix(0, ncol = 1),
      Z_PO_pres = if(bf.matrix.type == "sparse") {
        as(po.bf.matrix[pt.quad.id == 1, ], "sparseMatrix")
      } else {
        as.matrix(po.bf.matrix[pt.quad.id == 1, ])
      },
      Z_PO_quad = if(bf.matrix.type == "sparse") {
        as(po.bf.matrix[pt.quad.id == 0, ], "sparseMatrix")
      } else {
        as.matrix(po.bf.matrix[pt.quad.id == 0, ])
      },
      Z_PA = if(bf.matrix.type == "sparse") {
        as(matrix(0, ncol = 1), "sparseMatrix")
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

    # Starting Parameters #

    # if warm starts are provided
    if (!missing(starting.pars)) {
      if (class(starting.pars) == "scampr") {
        mod.obj <- starting.pars
        starting.pars <- lapply(split(mod.obj$par, names(mod.obj$par)), unname)
        if (model.type == "laplace") {
          starting.pars$random <- unname(mod.obj$random.effects[grepl("LP Posterior Mean", rownames(mod.obj$random.effects), fixed = T), 1L])
        }
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
      stop(paste0("coord.names, ", coord.names, ", not found in data set provided"))
    }
    ############################################################

    # Get the PA design matrix
    pa.des.mat <- scampr:::get.desgin.matrix(pa.formula, pa.data)
    fixed.names <- colnames(pa.des.mat)

    # Determine the basis functions to be used
    if (model.type != "ipp") {
      if (missing(simple.basis)) {
        if (missing(FRK.basis.functions)) { # When none is provided use FRK defaults (with 2 spatial resolutions)
          # create a spatial pixels data frame as required by FRK::auto_basis
          sp.data <- pa.data[ , coord.names]
          coordinates(sp.data) <- coord.names
          FRK.basis.functions <- FRK::auto_basis(data = sp.data, nres = 2)
        }
        pa.bf.matrix <- FRK::eval_basis(basis = FRK.basis.functions, as.matrix(pa.data[ , coord.names]))
        bf.info <- FRK.basis.functions@df
        colnames(bf.info)[grepl("loc", colnames(bf.info), fixed = T)] <- coord.names
        class(bf.info) <- c(class(bf.info), "bf.df")
      } else { # Otherwise use the provided simple basis
        pa.bf.matrix <- scampr:::get.bf.matrix(simple.basis, pa.data[ , coord.names])
        bf.info <- simple.basis
        FRK.basis.functions <- NULL
      }
    } else {
      pa.bf.matrix <- matrix(0, nrow = 1)
      bf.info <- cbind.data.frame(x = NA, y = NA, scale = NA, res = 1)
    }
    # TMB required data setup
    dat.list <- list(
      X_PO_pres = matrix(0, ncol = 1),
      B_PO_pres = matrix(0, ncol = 1),
      X_PO_quad = matrix(0, ncol = 1),
      B_PO_quad = matrix(0, ncol = 1),
      X_PA = as.matrix(pa.des.mat),
      Z_PO_pres = if(bf.matrix.type == "sparse") {
        as(matrix(0, ncol = 1), "sparseMatrix")
      } else {
        matrix(0, ncol = 1)
      },
      Z_PO_quad = if(bf.matrix.type == "sparse") {
        as(matrix(0, ncol = 1), "sparseMatrix")
      } else {
        matrix(0, ncol = 1)
      },
      Z_PA = if(bf.matrix.type == "sparse") {
        as(pa.bf.matrix, "sparseMatrix")
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

    # if warm starts are provided
    if (!missing(starting.pars)) {
      if (class(starting.pars) == "scampr") {
        mod.obj <- starting.pars
        starting.pars <- lapply(split(mod.obj$par, names(mod.obj$par)), unname)
        if (model.type == "laplace") {
          starting.pars$random <- unname(mod.obj$random.effects[grepl("LP Posterior Mean", rownames(mod.obj$random.effects), fixed = T), 1L])
        }
      }
    }
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
    ###################################################################
  } else {
    stop("Incorrect 'data.type' provided. Must be one of 'po', 'pa' or 'popa'")
  }
  return(list(tmb.data = dat.list, tmb.pars = start.pars))
}
