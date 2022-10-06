#' Combined Model for Presence-only and Presence/Absence Data
#'
#' @description Jointly fits a model to presence-only and presence/absence data as linked by response to environmental predictors provided in each formula. The presence-only formula must also contain biasing predictors to account for opportunistic collection.
#'
#' @param po.formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the Presence-only data model to be fitted. The 'response' must be a binary that indicates whether a datum is a presence or quadrature point. See GLM function for further formula details.
#' @param pa.formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the Presence/Absence data model to be fitted. All predictor terms must also be included in po.formula. The 'response' must be a must be either the site abundance or a binary indicating whether there is a presence or not. See GLM function for further formula details.
#' @param po.data a data frame containing predictors at both presence-records and quadrature as well as the po.formula 'response'.
#' @param pa.data a data frame containing predictors and response for the pa.formula.
#' @param coord.names a vector of character strings describing the column names of the coordinates in both data frames.
#' @param quad.weights.name a character string of the column name of quadrature weights in the po.data.
#' @param FRK.basis.functions an optional object of class 'Basis' from FRK package. If neither 'FRK.basis.functions' nor 'simple.basis' is specified will use default FRK::auto_basis with 2 spatial resolutions.
#' @param simple.basis an alternative to 'FRK.basis.functions': a data.frame of basis functions information created by 'simple_basis()'.
#' @param model.type a character string indicating the type of model to be used. May be one of 'laplace' or 'variational' for Cox Processes involving spatially correlated errors or 'ipp' for a model that follows an inhomgeneous Poisson process.
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to use sparse or dense matrix computations for the basis functions created.
#' @param se a logical indicating whether standard errors should be calculated.
#' @param starting.pars an optional named list or scampr model object that gives warm starting values for the parameters of the model.
#' @param additional.latent.field a logical indicating whether an additional latent field should be included in the Presence-only data model to capture additional detection bias. Only used in the integrated data model.
#'
#' @return a scampr model object
#' @noRd
#'
#' @importFrom methods as
#' @importFrom stats as.formula optim
#' @importFrom sp coordinates
#' @importFrom FRK auto_basis eval_basis
#' @importFrom TMB MakeADFun sdreport
#'
#' @examples
#' # Get the flora data for one of the species
#' dat_po <- flora$po$sp1
#' dat_pa <- flora$pa
#'
#' # Attach the quadrature to the PO data
#' dat_po <- rbind.data.frame(dat_po, flora$quad)
#'
#' # Fit without a shared latent field
#' m1 <- popa(pres ~ MNT + D.Main, sp1 ~ MNT,
#' po.data = dat_po, pa.data = dat_pa, model.type = "ipp")
#'
#' # Set up a simple 2D grid of basis functions to fit a LGCP model to the data
#' bfs <- simple_basis(nodes.on.long.edge = 9, data = dat_po)
#'
#' \dontrun{
#' # Fit with a shared latent field
#' m2 <- popa(pres ~ MNT + D.Main, Y ~ MNT,
#' po.data = dat_po, pa.data = dat_pa, simple.basis = bfs)
#' }
popa <- function(po.formula, pa.formula, po.data, pa.data, coord.names = c("x", "y"), quad.weights.name = "quad.size", FRK.basis.functions, simple.basis, model.type = c("laplace", "variational", "ipp"), bf.matrix.type = c("sparse", "dense"), se = TRUE, starting.pars, additional.latent.field = FALSE) {

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
  # if (!all(pa.pred %in% po.pred)) {
  #   stop("Not all predictors in the PA formula are found in the PO formula")
  # }
  if (!is.logical(se)) {
    stop(paste0("'se' must be a logcial indicating whether or not to calculate standard errors"))
  }
  # parameters of restricted strings
  model.type <- match.arg(model.type)
  if (model.type == "variational") {
    stop("Combined Data Model cannot handle variational approx. at this stage. Try model.type = 'laplace'")
  }
  bf.matrix.type <- match.arg(bf.matrix.type)

  ############################################################

  # default na.action is to remove any data rows with na (for terms involved in the model)
  rm.rows <- attr(na.omit(po.data[ , c(coord.names, quad.weights.name, po.resp, po.pred)]), "na.action")
  if (!is.null(rm.rows)) {
    po.data <- po.data[-rm.rows, ]
  }

  # default na.action is to remove any data rows with na (for terms involved in the model)
  rm.rows <- attr(na.omit(pa.data[ , c(coord.names, pa.resp, pa.pred)]), "na.action")
  if (!is.null(rm.rows)) {
    pa.data <- pa.data[-rm.rows, ]
  }

  # Get the PA design matrix
  pa.des.mat <- get.design.matrix(pa.formula, pa.data)

  # Need to adjust the breakdown of formulas depending on whether bias correcting predictors have been included for the PO data
  # create a switch object for this:
  if (all(po.pred %in% pa.pred)) {
    biasing.predictors <- "absent"
  } else {
    biasing.predictors <- "present"
  }

  # Get the full PO design matrix
  full.po.des.mat <- get.design.matrix(po.formula, po.data)
  # Get the separate predictors (and remove intercepts if present - to be dealt with later)
  po.preds <- colnames(full.po.des.mat)[colnames(full.po.des.mat) %in% colnames(pa.des.mat)]
  po.preds <- po.preds[!grepl("Intercept", po.preds, fixed = T)]
  bias.preds <- colnames(full.po.des.mat)[!colnames(full.po.des.mat) %in% colnames(pa.des.mat)]
  bias.preds <- bias.preds[!grepl("Intercept", bias.preds, fixed = T)]

  # set a switch for when the model uses biasing terms including an intercept (needed for TMB objective function)
  biasing.terms <- "present" # will assume present until switched off in the event that there are no intercepts, nor bias predictors

  # Add/remove intercept term as required - unless po.formula and pa.formula DO NOT contain an intercept, the PO data will be modelled with a separate intercept
  if (any(grepl("Intercept", colnames(full.po.des.mat), fixed = T))) {
    if (!any(grepl("Intercept", colnames(pa.des.mat), fixed = T))) {
      po.pred.formula <- stats::as.formula(paste0(po.resp, " ~ - 1 + ", paste(po.preds, collapse = " + ")))
      po.bias.formula <- switch(biasing.predictors,
                                present = stats::as.formula(paste0(po.resp, " ~ ", paste(bias.preds, collapse = " + "))),
                                absent = stats::as.formula(paste0(po.resp, " ~ 1"))
      )
    } else {
      po.pred.formula <- stats::as.formula(paste0(po.resp, " ~ ", paste(po.preds, collapse = " + ")))
      po.bias.formula <- switch(biasing.predictors,
                                present = stats::as.formula(paste0(po.resp, " ~ ", paste(bias.preds, collapse = " + "))),
                                absent = stats::as.formula(paste0(po.resp, " ~ 1"))
      )
    }
    # Get the PO design matrices in this scenario
    po.des.mat <- get.design.matrix(po.pred.formula, po.data)
    po.bias.des.mat <- get.design.matrix(po.bias.formula, po.data)
    # Set an identifier for the biasing predictors in this scenario
    # fixed.names.bias.id <- switch(biasing.predictors,
    #                              present = c(rep(F, ncol(po.des.mat)), rep(T, ncol(po.bias.des.mat))),
    #                              absent = c(rep(F, ncol(po.des.mat)), T)
    # )
  } else {
    # if (any(grepl("Intercept", colnames(pa.des.mat), fixed = T))) {
    #   stop("Intercept term found in pa.formula and not in po.formula")
    # }
    po.pred.formula <- stats::as.formula(paste0(po.resp, " ~ - 1 + ", paste(po.preds, collapse = " + ")))
    po.bias.formula <- switch(biasing.predictors,
                              present = stats::as.formula(paste0(po.resp, " ~ - 1 + ", paste(bias.preds, collapse = " + "))),
                              absent = NA
    )
    # Get the PO design matrices in this scenario
    po.des.mat <- get.design.matrix(po.pred.formula, po.data)
    po.bias.des.mat <- switch(biasing.predictors,
                              present = get.design.matrix(po.bias.formula, po.data),
                              absent = matrix(rep(0, nrow(po.des.mat)), ncol = 1)
    )
    # Set an identifier for the biasing predictors in this scenario
    # fixed.names.bias.id <- switch(biasing.predictors,
    #                               present = c(rep(F, ncol(po.des.mat)), rep(T, ncol(po.bias.des.mat))),
    #                               absent = rep(F, ncol(po.des.mat))
    # )
    # this scenarios flicks the switch for mapping biasing terms (within this scope there are none if there are no biasing predictors)
    biasing.terms <- biasing.predictors

  }

  # Adjust the Intercept name if required
  if (any(grepl("Intercept", colnames(po.bias.des.mat), fixed = T))) {
    colnames(po.bias.des.mat)[grepl("Intercept", colnames(po.bias.des.mat), fixed = T)] <- "(Bias Intercept)"
  }
  # Get the presence point/ quadrature point identifier
  pt.quad.id <- po.data[ , po.resp]
  # Set the fixed effect names
  fixed.names <- c(colnames(po.des.mat), colnames(pa.des.mat), colnames(po.bias.des.mat))
  fixed.names.bias.id <- c(rep(F, ncol(po.des.mat)), rep(F, ncol(pa.des.mat)), rep(T, ncol(po.bias.des.mat)))
  remove.terms.id <- !duplicated(fixed.names)
  fixed.names <- fixed.names[remove.terms.id]
  fixed.names.bias.id <- fixed.names.bias.id[remove.terms.id]

  # Determine the basis functions to be used
  if (model.type != "ipp") {
    if (missing(simple.basis)) {
      if (missing(FRK.basis.functions)) { # When none is provided use FRK defaults (with 2 spatial resolutions)
        # create a spatial pixels data frame as required by FRK::auto_basis
        sp.data <- rbind(po.data[ ,coord.names], pa.data[ , coord.names])
        sp::coordinates(sp.data) <- coord.names
        FRK.basis.functions <- FRK::auto_basis(data = sp.data, max_basis = sum(pt.quad.id)*0.25)
      }
      po.bf.matrix <- FRK::eval_basis(basis = FRK.basis.functions, as.matrix(po.data[ , coord.names]))
      pa.bf.matrix <- FRK::eval_basis(basis = FRK.basis.functions, as.matrix(pa.data[ , coord.names]))
      bf.info <- FRK.basis.functions@df
      colnames(bf.info)[grepl("loc", colnames(bf.info), fixed = T)] <- coord.names
      class(bf.info) <- c(class(bf.info), "bf.df")
    } else { # Otherwise use the provided simple basis
      po.bf.matrix <- get.bf.matrix(simple.basis, po.data[ , coord.names])
      pa.bf.matrix <- get.bf.matrix(simple.basis, pa.data[ , coord.names])
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
    quad_size = po.data[pt.quad.id == 0, quad.weights.name],
    Y = pa.data[ , pa.resp] > 0, # corrects in the case of abundance
    bf_per_res = as.numeric(table(bf.info$res)),
    mod_type = as.integer(which(model.type == c("ipp", "variational", "laplace")) - 1),
    data_type = 2,
    bf_matrix_type = bf.matrix.type
  )

  # # AT THIS STAGE CAN ONLY PERFORM LAPLAC APPROX.
  # # create the appropriate start parameters for the variance component
  # #   w.r.t. approx. type and data type
  # var.starts <- switch(model.type,
  #                      variational = rep(0, ncol(dat.list$Z_PO_pres)),
  #                      laplace = rep(0, length(dat.list$bf_per_res))
  # )
  var.starts <- rep(0, length(dat.list$bf_per_res))
  start.pars <- list(fixed = rep(0, ncol(dat.list$X_PA)),
                     bias = rep(0, ncol(dat.list$B_PO_pres)),
                     random = rep(0, ncol(dat.list$Z_PA)),
                     random_bias = rep(0, ncol(dat.list$Z_PO_pres)),
                     log_variance_component = var.starts,
                     log_variance_component_bias = var.starts
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

  # AT THIS STAGE CAN ONLY PERFORM LAPLACE APPROX.
  # set up the objective function w.r.t. model.type
  if (additional.latent.field) {
    if (biasing.terms == "present") {
      obj <- switch(model.type,
                    ipp = TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", map = list(random = factor(rep(NA, length(start.pars$random))), random_bias = factor(rep(NA, length(start.pars$random_bias))), log_variance_component = factor(rep(NA, length(start.pars$log_variance_component))), log_variance_component_bias = factor(rep(NA, length(start.pars$log_variance_component_bias)))), silent = T),
                    variational = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = c("random", "random_bias"), DLL = "scampr", silent = T),
                    laplace = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = c("random", "random_bias"), DLL = "scampr", silent = T)
      )
    } else {
      # ensure bias parameters are at zero for map
      start.pars$bias <- rep(0, ncol(dat.list$B_PO_pres))
      obj <- switch(model.type,
                    ipp = TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", map = list(bias = factor(rep(NA, length(start.pars$bias))), random = factor(rep(NA, length(start.pars$random))), log_variance_component = factor(rep(NA, length(start.pars$log_variance_component)))), silent = T),
                    variational = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = c("random", "random_bias"), DLL = "scampr", map = list(bias = factor(rep(NA, length(start.pars$bias)))), silent = T),
                    laplace = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = c("random", "random_bias"), DLL = "scampr", map = list(bias = factor(rep(NA, length(start.pars$bias)))), silent = T)
      )
    }
  } else {
    # ensure random_bias parameters are at zero for map
    start.pars$random_bias <- rep(0, ncol(dat.list$Z_PO_pres))
    # start.pars$log_variance_component_bias <- -1e6 # set near enough to zero on the exponential-scale
    if (biasing.terms == "present") {
      obj <- switch(model.type,
                    ipp = TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", map = list(random = factor(rep(NA, length(start.pars$random))), random_bias = factor(rep(NA, length(start.pars$random_bias))), log_variance_component = factor(rep(NA, length(start.pars$log_variance_component))), log_variance_component_bias = factor(rep(NA, length(start.pars$log_variance_component_bias)))), silent = T),
                    variational = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", map = list(random_bias = factor(rep(NA, length(start.pars$random_bias))), log_variance_component_bias = factor(rep(NA, length(start.pars$log_variance_component_bias)))), silent = T),
                    laplace = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", map = list(random_bias = factor(rep(NA, length(start.pars$random_bias))), log_variance_component_bias = factor(rep(NA, length(start.pars$log_variance_component_bias)))), silent = T)
      )
    } else {
      # ensure bias parameters are at zero for map
      start.pars$bias <- rep(0, ncol(dat.list$B_PO_pres))
      obj <- switch(model.type,
                    ipp = TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", map = list(bias = factor(rep(NA, length(start.pars$bias))), random = factor(rep(NA, length(start.pars$random))), random_bias = factor(rep(NA, length(start.pars$random_bias))), log_variance_component = factor(rep(NA, length(start.pars$log_variance_component_bias))), log_variance_component_bias = factor(rep(NA, length(start.pars$log_variance_component_bias)))), silent = T),
                    variational = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", map = list(bias = factor(rep(NA, length(start.pars$bias))), random_bias = factor(rep(NA, length(start.pars$random_bias))), log_variance_component_bias = factor(rep(NA, length(start.pars$log_variance_component_bias)))), silent = T),
                    laplace = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", map = list(bias = factor(rep(NA, length(start.pars$bias))), random_bias = factor(rep(NA, length(start.pars$random_bias))), log_variance_component_bias = factor(rep(NA, length(start.pars$log_variance_component_bias)))), silent = T)
      )
    }
  }

  # # AT THIS STAGE CAN ONLY PERFORM LAPLACE APPROX.
  # # set up the objective function w.r.t. model.type
  # if (additional.latent.field) {
  #   if (biasing.predictors == "present") {
  #     obj <- switch(model.type,
  #                   ipp = TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", map = list(random = factor(rep(NA, length(start.pars$random))), random_bias = factor(rep(NA, length(start.pars$random_bias))), log_variance_component = factor(rep(NA, length(start.pars$log_variance_component))), log_variance_component_bias = factor(rep(NA, length(start.pars$log_variance_component_bias)))), silent = T),
  #                   variational = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = c("random", "random_bias"), DLL = "scampr", silent = T),
  #                   laplace = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = c("random", "random_bias"), DLL = "scampr", silent = T)
  #     )
  #   } else {
  #     # ensure bias parameters are at zero for map
  #     start.pars$bias <- rep(0, ncol(dat.list$B_PO_pres))
  #     obj <- switch(model.type,
  #                   ipp = TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", map = list(bias = factor(rep(NA, length(start.pars$bias))), random = factor(rep(NA, length(start.pars$random))), log_variance_component = factor(rep(NA, length(start.pars$log_variance_component)))), silent = T),
  #                   variational = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = c("random", "random_bias"), DLL = "scampr", map = list(bias = factor(rep(NA, length(start.pars$bias)))), silent = T),
  #                   laplace = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = c("random", "random_bias"), DLL = "scampr", map = list(bias = factor(rep(NA, length(start.pars$bias)))), silent = T)
  #     )
  #   }
  # } else {
  #   # ensure random_bias parameters are at zero for map
  #   start.pars$random_bias <- rep(0, ncol(dat.list$Z_PO_pres))
  #   if (biasing.predictors == "present") {
  #     obj <- switch(model.type,
  #                   ipp = TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", map = list(random = factor(rep(NA, length(start.pars$random))), random_bias = factor(rep(NA, length(start.pars$random_bias))), log_variance_component = factor(rep(NA, length(start.pars$log_variance_component))), log_variance_component_bias = factor(rep(NA, length(start.pars$log_variance_component_bias)))), silent = T),
  #                   variational = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", map = list(random_bias = factor(rep(NA, length(start.pars$random_bias))), log_variance_component_bias = factor(rep(NA, length(start.pars$log_variance_component_bias)))), silent = T),
  #                   laplace = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", map = list(random_bias = factor(rep(NA, length(start.pars$random_bias))), log_variance_component_bias = factor(rep(NA, length(start.pars$log_variance_component_bias)))), silent = T)
  #     )
  #   } else {
  #     # ensure bias parameters are at zero for map
  #     start.pars$bias <- rep(0, ncol(dat.list$B_PO_pres))
  #     obj <- switch(model.type,
  #                   ipp = TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", map = list(bias = factor(rep(NA, length(start.pars$bias))), random = factor(rep(NA, length(start.pars$random))), random_bias = factor(rep(NA, length(start.pars$random_bias))), log_variance_component = factor(rep(NA, length(start.pars$log_variance_component_bias))), log_variance_component_bias = factor(rep(NA, length(start.pars$log_variance_component_bias)))), silent = T),
  #                   variational = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", map = list(bias = factor(rep(NA, length(start.pars$bias))), random_bias = factor(rep(NA, length(start.pars$random_bias))), log_variance_component_bias = factor(rep(NA, length(start.pars$log_variance_component_bias)))), silent = T),
  #                   laplace = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", map = list(bias = factor(rep(NA, length(start.pars$bias))), random_bias = factor(rep(NA, length(start.pars$random_bias))), log_variance_component_bias = factor(rep(NA, length(start.pars$log_variance_component_bias)))), silent = T)
  #     )
  #   }
  # }

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
  if (additional.latent.field) {
    coef.names <- switch(model.type,
                         ipp = fixed.names,
                         variational = c(fixed.names, paste0("VA Posterior Mean (bf ", random.nos, ")"), paste0("VA Posterior log sd (bf ", random.nos, ")")),
                         laplace = c(fixed.names, paste0("Prior log sd (res. ", 1:length(dat.list$bf_per_res), ")"), paste0("Bias Fld. Prior log sd (res. ", 1:length(dat.list$bf_per_res), ")"))
    )
  } else {
    coef.names <- switch(model.type,
                         ipp = fixed.names,
                         variational = c(fixed.names, paste0("VA Posterior Mean (bf ", random.nos, ")"), paste0("VA Posterior log sd (bf ", random.nos, ")")),
                         laplace = c(fixed.names, paste0("Prior log sd (res. ", 1:length(dat.list$bf_per_res), ")"))
    )
  }

  names(res$coefficients) <- coef.names
  # check for a single fixed effect to adjust the resulting data frame
  if (length(fixed.names) == 1) {
    res$fixed.effects <- data.frame(t(tmp.estimates[1:length(fixed.names), ]))
    colnames(res$fixed.effects)[2] <- "Std. Error" # need to correct the white space in name
  } else {
    res$fixed.effects <- tmp.estimates[1:length(fixed.names), ]
  }
  rownames(res$fixed.effects) <- fixed.names
  if (model.type != "ipp") {
    res$random.effects <- tmp.estimates[(length(fixed.names) + 1):nrow(tmp.estimates), ]
    res$random.effects <- res$random.effects[!grepl("log_", rownames(res$random.effects), fixed = T), ]
    res$random.effects <- res$random.effects[!grepl("_bias", rownames(res$random.effects), fixed = T), ] # additionally remove the bias field coefficients and variance
    rand.names <- switch(model.type,
                         variational = c(paste0("VA Posterior Mean (bf ", random.nos, ")"), paste0("VA Posterior Var (bf ", random.nos, ")"), paste0("Prior Var (res. ", 1:length(dat.list$bf_per_res), ")")),
                         laplace = c(paste0("LP Posterior Mean (bf ", random.nos, ")"), paste0("Prior Var (res. ", 1:length(dat.list$bf_per_res), ")"))
    )
    if (additional.latent.field) {
      res$bias.field <- tmp.estimates[grepl("_bias", rownames(tmp.estimates), fixed = T), ]
      res$bias.field <- res$bias.field[!grepl("log_", rownames(res$bias.field), fixed = T), ]
      bias.rand.names <- switch(model.type,
                           variational = c(paste0("VA Posterior Mean (bf ", random.nos, ")"), paste0("VA Posterior Var (bf ", random.nos, ")"), paste0("Prior Var (res. ", 1:length(dat.list$bf_per_res), ")")),
                           laplace = c(paste0("Bias Fld. LP Posterior Mean (bf ", random.nos, ")"), paste0("Bias Fld. Prior Var (res. ", 1:length(dat.list$bf_per_res), ")"))
      )
      rownames(res$bias.field) <- bias.rand.names
    }
    rownames(res$random.effects) <- rand.names
    res$basis.per.res <- dat.list$bf_per_res
    res$FRK.basis.functions <- FRK.basis.functions
    res$basis.fn.info <- bf.info
    res$approx.type <- "laplace" #model.type until POPA models can handle VA
    if (additional.latent.field) {
      res$fitted.values <- switch(biasing.terms,
                                  present = as.vector(po.des.mat %*% res$fixed.effects[na.omit(match(fixed.names, colnames(po.des.mat))), 1] + po.bias.des.mat %*% res$fixed.effects[na.omit(match(fixed.names, colnames(po.bias.des.mat))), 1] + po.bf.matrix %*% res$random.effects[grepl(" Mean ", rownames(res$random.effects), fixed = T), 1] + po.bf.matrix %*% res$bias.field[grepl(" Mean ", rownames(res$bias.field), fixed = T), 1]),
                                  absent = as.vector(po.des.mat %*% res$fixed.effects[na.omit(match(fixed.names, colnames(po.des.mat))), 1] + po.bf.matrix %*% res$random.effects[grepl(" Mean ", rownames(res$random.effects), fixed = T), 1] + po.bf.matrix %*% res$bias.field[grepl(" Mean ", rownames(res$bias.field), fixed = T), 1])
      )
      attr(res$fitted.values, "abundance") <- as.vector(pa.des.mat %*% res$fixed.effects[na.omit(match(fixed.names, colnames(pa.des.mat))), 1] + pa.bf.matrix %*% res$random.effects[grepl(" Mean ", rownames(res$random.effects), fixed = T), 1])
    } else {
      res$fitted.values <- switch(biasing.terms,
                                  present = as.vector(po.des.mat %*% res$fixed.effects[na.omit(match(fixed.names, colnames(po.des.mat))), 1] + po.bias.des.mat %*% res$fixed.effects[na.omit(match(fixed.names, colnames(po.bias.des.mat))), 1] + po.bf.matrix %*% res$random.effects[grepl(" Mean ", rownames(res$random.effects), fixed = T), 1]),
                                  absent = as.vector(po.des.mat %*% res$fixed.effects[na.omit(match(fixed.names, colnames(po.des.mat))), 1] + po.bf.matrix %*% res$random.effects[grepl(" Mean ", rownames(res$random.effects), fixed = T), 1])
      )
      attr(res$fitted.values, "abundance") <- as.vector(pa.des.mat %*% res$fixed.effects[na.omit(match(fixed.names, colnames(pa.des.mat))), 1] + pa.bf.matrix %*% res$random.effects[grepl(" Mean ", rownames(res$random.effects), fixed = T), 1])
    }
  } else {
    res$random.effects <- NA
    res$basis.per.res <- NA
    res$FRK.basis.functions <- NULL
    res$basis.fn.info <- NULL
    res$approx.type <- NA
    res$fitted.values <- switch(biasing.terms,
                                present = as.vector(po.des.mat %*% res$fixed.effects[na.omit(match(fixed.names, colnames(po.des.mat))), 1] + po.bias.des.mat %*% res$fixed.effects[na.omit(match(fixed.names, colnames(po.bias.des.mat))), 1]),
                                absent = as.vector(po.des.mat %*% res$fixed.effects[na.omit(match(fixed.names, colnames(po.des.mat))), 1])
    )
    attr(res$fitted.values, "abundance") <- as.vector(pa.des.mat %*% res$fixed.effects[na.omit(match(fixed.names, colnames(pa.des.mat))), 1])
  }
  res$starting.pars <- start.pars
  res$data <- po.data
  attr(res$data, "pa") <- pa.data
  res$formula <- po.formula
  attr(res$formula, "pa") <- pa.formula
  res$coord.names <- coord.names
  res$quad.weights.name <- quad.weights.name
  res$pt.quad.id <- pt.quad.id
  res$data.model.type <- "popa"
  res$bf.matrix.type <- bf.matrix.type
  class(res) <- "scampr"
  return(res)
}
