#' Combined Model for Presence-only and Presence-Absence Data
#'
#' @param po.formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the Presence-only data model to be fitted. The 'response' must be a binary that indicates whether a datum is a presence or quadrature point. See GLM function for further formula details.
#' @param pa.formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the Presence/Absence data model to be fitted. The 'response' must be a must be either the site abundance or a binary indicating whether there is a presence or not. See GLM function for further formula details.
#' @param po.data a data frame containing predictors at both presence-records and quadrature as well as the po.formula 'response'.
#' @param pa.data a data frame containing predictors and response for the pa.formula.
#' @param coord.names a vector of character strings describing the column names of the coordinates in both data frames
#' @param quad.weights.name a charater string of the column name of quadrature weights in the po.data
#' @param FRK.basis.functions an optional object of class 'Basis' from FRK package. If neither 'FRK.basis.functions' nor 'simple.basis' is specified will use default FRK::auto_basis with 2 spatial resolutions
#' @param simple.basis an alternative to 'FRK.basis.functions': a data.frame of basis functions information created by 'simple_basis()'.
#' @param model.type a character string indicating the type of model to be used. May be one of 'laplace' or 'variational' for Cox Processes involving spatially correlated errors or 'ipp' for a model that follows an inhomgeneous Poisson process.
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to use sparse or dense matrix computations for the basis functions created.
#' @param se a logical indicating whether standard errors should be calculated
#' @param starting.pars an optional named list or scampr model that gives warm starting values for the parameters of the model.
#' @param subset an optional subset of the data to be used.
#' @param na.action an optional way of handling NA's in the data, default is omit.
#'
#' @return
#' @export
#'
#' @examples
popa <- function(po.formula, pa.formula, po.data, pa.data, coord.names = c("x", "y"), quad.weights.name = "quad.size", FRK.basis.functions, simple.basis, model.type = c("laplace", "variational", "ipp"), bf.matrix.type = c("sparse", "dense"), se = TRUE, starting.pars, subset, na.action) {

  # CAN'T JUST GIVE A BASIS FUNCTION MATRIX BECAUSE THEN YOU CAN'T PREDICT ETC. AS WE DON'T KNOW ENOUGH ABOUT THE FUNCTIONS

  # Get the response and predictor names
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
  if (!is.logical(se)) {
    stop(paste0("'se' must be a logcial indicating whether or not to calculate standard erros"))
  }
  # parameters of restricted strings
  model.type <- match.arg(model.type)
  bf.matrix.type <- match.arg(bf.matrix.type)

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

  # if starting.pars provided is a scampr model adjust to req. list structure
  if (!missing(starting.pars)) {
    if (class(starting.pars) == "scampr") {
      starting.pars <- lapply(split(starting.pars$par, names(starting.pars$par)), unname)
    }
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

  # # AT THIS STAGE CAN ALONE PERFORM LAPLAC APPROX.
  # # create the appropraite start parameters for the variance component
  # #   w.r.t. approx. type and data type
  # var.starts <- switch(model.type,
  #                      variational = rep(0, ncol(dat.list$Z_PO_pres)),
  #                      laplace = rep(0, length(dat.list$bf_per_res))
  # )
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

  # # AT THIS STAGE CAN ONLY PERFORM LAPLAC APPROX.
  # set up the objective function w.r.t. model.type
  obj <- switch(model.type,
                ipp = TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", map = list(random = factor(rep(NA, ncol(dat.list$Z_PA))), log_variance_component = factor(rep(NA, length(var.starts)))), silent = T),
                variational = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", silent = T),
                laplace = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", silent = T)
  )
  # obj <- TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", silent = T)
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
    res$approx.type <- model.type
    res$fitted.values <- as.vector(po.des.mat %*% res$fixed.effects[!fixed.names.bias.id, 1] + po.bias.des.mat %*% res$fixed.effects[fixed.names.bias.id, 1] + po.bf.matrix %*% res$random.effects[grepl(" Mean ", rownames(res$random.effects), fixed = T), 1])
    attr(res$fitted.values, "abundance") <- as.vector(pa.des.mat %*% res$fixed.effects[!fixed.names.bias.id, 1] + pa.bf.matrix %*% res$random.effects[grepl(" Mean ", rownames(res$random.effects), fixed = T), 1])
  } else {
    res$random.effects <- NA
    res$basis.per.res <- NA
    res$FRK.basis.functions <- NULL
    res$basis.fn.info <- NULL
    res$approx.type <- NA
    res$fitted.values <- as.vector(po.des.mat %*% res$fixed.effects[!fixed.names.bias.id, 1] + po.bias.des.mat %*% res$fixed.effects[fixed.names.bias.id, 1])
    attr(res$fitted.values, "abundance") <- as.vector(pa.des.mat %*% res$fixed.effects[!fixed.names.bias.id, 1])
  }
  res$starting.pars <- start.pars
  res$data <- po.data
  attr(res$data, "pa") <- pa.data
  res$formula <- paste0(po.resp, " ~ ", paste(po.pred, collapse = " + "), " |&| ", pa.resp, " ~ ", paste(pa.pred, collapse = " + "))
  res$coord.names <- coord.names
  res$quad.weights.name <- quad.weights.name
  res$pt.quad.id <- pt.quad.id
  res$data.model.type <- "popa"
  res$bf.matrix.type <- bf.matrix.type
  class(res) <- "scampr"
  return(res)
}
