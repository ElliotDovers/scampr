#' Internal scampr function that creates a list of data and starting parameters for scampr models
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the fixed effects of the model to be fitted. The 'response' must be a binary that indicates whether a datum is a presence or: quadrature point (for point process models)/ absence (for binary models). See GLM function for further formula details.
#' @param data a data frame containing response and predictors within \code{formula}.
#' @param bias.formula an object of class "formula" (or one that can be coerced to that class) OR the character string "latent". In the formula case, this is a symbolic description of the predictors included to account for bias in the presence-only data (no response term is needed). In the case of fitting an integrated data model, \code{bias.formula = "latent"} will fit an approximate latent Gaussian field to account for the bias.
#' @param IDM.presence.absence.df an optional data frame. When fitting an integrated data model use this to pass in the presence/absence data.
#' @param coord.names a vector of character strings describing the column names of the coordinates in both data frames.
#' @param quad.weights.name a charater string of the column name of quadrature weights in the data.
#' @param approx.type a character string indicating the type of model to be used. May be one of 'laplace' or 'variational' for models involving spatial random effects, OR 'not_sre' for a fixed effect model.
#' @param model.type a character string indicating the type of data to be used. May be one of 'PO' (for a presence-only PPM) or 'PA' (for a presence/absence Binary GLM) or 'IDM' (for an integrated data model).
#' @param basis.functions an optional object of class 'Basis' created by \code{FRK::auto_basis()} or 'bf.df' created by \code{scampr::simple_basis()}. Either object describes a set of basis functions for approximating the latent Gaussian field. If NULL the model will use default \code{FRK::auto_basis()} with \code{max_basis = 0.25 * # of points}.
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to use sparse or dense matrix computations for the basis functions created.
#' @param starting.pars an optional named list or scampr model object that gives warm starting values for the parameters of the model.
#' @param latent.po.biasing a logical indicating whether biasing in the presence-only data should be accounted for via an additional latent field. Applies to IDM only.
#' @param po.biasing.basis.functions an optional extra set of basis functions that can be used when \code{latent.po.biasing = TRUE}, otherwise \code{basis.functions} are used.
#' @param prune.bfs an integer indicating the number of presence-only records required within a basis function's radius for it NOT to be pruned. Applies to the PO and IDM model (additionally, within the presence-only biasing basis functions in the IDM case) to assist with stability in model convergence. Default is zero, i.e. no pruning.
#'
#' @return list of elements required for TMB::MakeADFun
#' @export
#'
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
#' tmb.input <- scampr:::get.TMB.data.input(pres ~ MNT + D.Main, sp1 ~ MNT, po.data = dat_po, pa.data = dat_pa, approx.type = "not_sre")
#' str(tmb.input)
get.TMB.data.input <- function(formula, data, bias.formula, IDM.presence.absence.df, coord.names = c("x", "y"), quad.weights.name = "quad.size", approx.type = c("variational", "laplace", "not_sre"), model.type = c("PO", "PA", "IDM"), basis.functions, bf.matrix.type = c("sparse", "dense"), starting.pars, latent.po.biasing = TRUE, po.biasing.basis.functions, prune.bfs = 4) {

  # checks for parameters of restricted strings
  model.type <- match.arg(model.type)
  approx.type <- match.arg(approx.type)
  bf.matrix.type <- match.arg(bf.matrix.type)

  # store the arguments
  # arg.info <- as.list(environment())
  arg.info <- as.list(match.call())

  # approach based on the data type #

  if (model.type == "PO") { # Presence-only data

    # get the presence point/ quadrature point identifier
    pt.quad.id <- as.numeric(data[ , all.vars(formula[[2]])])

    # track row number splits on presence and quadrature points
    row.id <- 1:nrow(data)
    pres.rows <- row.id[pt.quad.id == 1]
    quad.rows <- row.id[pt.quad.id == 0]

    # split the data into presence points and quadrature
    data_pres <- data[pt.quad.id == 1, ]
    data_quad <- data[pt.quad.id == 0, ]

    ## Fixed Effects ###########################################################

    # get the fixed effect design matrix
    des.mat <- get.design.matrix(formula, data)
    # need to adjust for a single column design matrix
    if (ncol(des.mat) == 1) {
      des.mat_pres <- as.matrix(data.frame(des.mat[pt.quad.id == 1, ]))
      des.mat_quad <- as.matrix(data.frame(des.mat[pt.quad.id == 0, ]))
      colnames(des.mat_pres) <- colnames(des.mat)
      colnames(des.mat_quad) <- colnames(des.mat)
    } else {
      # need to adjust for a single row design matrix
      if (sum(pt.quad.id == 1) == 1) {
        des.mat_pres <- as.matrix(as.data.frame(t(des.mat[pt.quad.id == 1, ])))
      } else {
        des.mat_pres <- des.mat[pt.quad.id == 1, ]
      }
      if (sum(pt.quad.id == 0) == 1) {
        des.mat_quad <- as.matrix(as.data.frame(t(des.mat[pt.quad.id == 0, ])))
      } else {
        des.mat_quad <- des.mat[pt.quad.id == 0, ]
      }
    }

    # get the bias predictor design matrix
    if (missing(bias.formula)) {
      fixed.bias.type <- "missing"
    } else if (is(bias.formula, "formula")) {
      bias.des.mat <- get.design.matrix(bias.formula, data)
      # need to adjust for a single column design matrix
      if (ncol(bias.des.mat) == 1) {
        bias.des.mat_pres <- as.matrix(data.frame(bias.des.mat[pt.quad.id == 1, ]))
        bias.des.mat_quad <- as.matrix(data.frame(bias.des.mat[pt.quad.id == 0, ]))
        colnames(bias.des.mat_pres) <- colnames(bias.des.mat)
        colnames(bias.des.mat_quad) <- colnames(bias.des.mat)
      } else {
        # need to adjust for a single row design matrix
        if (sum(pt.quad.id == 1) == 1) {
          bias.des.mat_pres <- as.matrix(as.data.frame(t(bias.des.mat[pt.quad.id == 1, ])))
        } else {
          bias.des.mat_pres <- bias.des.mat[pt.quad.id == 1, ]
        }
        if (sum(pt.quad.id == 0) == 1) {
          bias.des.mat_quad <- as.matrix(as.data.frame(t(bias.des.mat[pt.quad.id == 0, ])))
        } else {
          bias.des.mat_quad <- bias.des.mat[pt.quad.id == 0, ]
        }
      }
      # Adjust the Intercept names if required
      if (any(grepl("(Intercept)", colnames(bias.des.mat_pres), fixed = T))) {
        colnames(bias.des.mat_pres)[grepl("(Intercept)", colnames(bias.des.mat_pres), fixed = T)] <- "(Bias Intercept)"
      }
      if (any(grepl("(Intercept)", colnames(bias.des.mat_quad), fixed = T))) {
        colnames(bias.des.mat_quad)[grepl("(Intercept)", colnames(bias.des.mat_quad), fixed = T)] <- "(Bias Intercept)"
      }
      fixed.bias.type <- "covariates"
    } else {
      stop(paste0("'bias.formula' provided is not of class formula"))
    }
    ############################################################################

    ## Random Effects ##########################################################

    # Determine the basis functions to be used
    if (approx.type != "not_sre") {
      # when no basis functions are provided use a simple basis default
      if (missing(basis.functions)) {
        # get a rough guide for the number of basis functions (to be 50% of the presence points)
        sqrt_number_bfs <- sqrt(sum(pt.quad.id)*0.5)
        # set the basis function
        basis.functions <- simple_basis(sqrt_number_bfs, data, coord.names = coord.names)
      }
      # calculate the basis function matrices
      po.bf.matrix_pres <- get.bf.matrix(basis.functions, point.locations = data_pres[ , coord.names], bf.matrix.type = bf.matrix.type)
      po.bf.matrix_quad <- get.bf.matrix(basis.functions, point.locations = data_quad[ , coord.names], bf.matrix.type = bf.matrix.type)
      # store the basis function information
      bf.info <- attr(po.bf.matrix_pres, "bf.df")

      # prune the basis functions if required
      if (prune.bfs != 0) {
        # determine basis functions that do not intersect any presence points
        prune.id <- do.call("apply", list(po.bf.matrix_pres, MARGIN = 2, FUN = function(x){sum(x>0)})) < prune.bfs
        # check that the pruning doesn't remove all basis functions
        while (all(prune.id)) {
          prune.bfs <- prune.bfs - 1
          prune.id <- do.call("apply", list(po.bf.matrix_pres, MARGIN = 2, FUN = function(x){sum(x>0)})) < prune.bfs
          warning(paste0("'prune.bfs' = ", prune.bfs + 1, " results in no valid basis functions. Trying 'prune.bfs' = ", prune.bfs, ". Otherwise fit a model without SRE."))
        }
        # prune from both basis function matrices - adjusting for the case when only 1 bf remains
        if (sum(!prune.id) == 1) {
          po.bf.matrix_pres <- as.matrix(data.frame(po.bf.matrix_pres[ , !prune.id]))
          po.bf.matrix_quad <- as.matrix(data.frame(po.bf.matrix_quad[ , !prune.id]))
          colnames(po.bf.matrix_pres) <- NULL
          colnames(po.bf.matrix_quad) <- NULL
          # re-adjust to sparse matrix if required
          if (bf.matrix.type == "sparse") {
            po.bf.matrix_pres <- methods::as(po.bf.matrix_pres, "sparseMatrix")
            po.bf.matrix_quad <- methods::as(po.bf.matrix_quad, "sparseMatrix")
          }
        } else {
          # adjusting for the case of a single presence or quad point:
          if (nrow(po.bf.matrix_pres) == 1) {
            po.bf.matrix_pres <- matrix(po.bf.matrix_pres[ , !prune.id], 1)
            if (bf.matrix.type == "sparse") {
              po.bf.matrix_pres <- methods::as(po.bf.matrix_pres, "sparseMatrix")
            }
          } else {
            po.bf.matrix_pres <- po.bf.matrix_pres[ , !prune.id]
          }
          if (nrow(po.bf.matrix_quad) == 1) {
            po.bf.matrix_quad <- matrix(po.bf.matrix_quad[ , !prune.id], 1)
            if (bf.matrix.type == "sparse") {
              po.bf.matrix_quad <- methods::as(po.bf.matrix_quad, "sparseMatrix")
            }
          } else {
            po.bf.matrix_quad <- po.bf.matrix_quad[ , !prune.id]
          }
        }
        # adjust the basis function information
        tmp <- bf.info
        bf.info <- tmp[!prune.id, ]
        attr(bf.info, "pruned") <- tmp[prune.id, ]
        # adjust the supplied basis functions depending on whether they are simple_basis or FRK package
        if (is(basis.functions, "bf.df")) {
          basis.functions <- basis.functions[!prune.id, ]
        } else {
          prune.idx <- (1:length(prune.id))[prune.id]
          basis.functions <- FRK::remove_basis(basis.functions, prund.idx)
        }
      }

    } else {
      # set a trivial example for not_sre models
      po.bf.matrix_pres <- matrix(rep(0, nrow(des.mat_pres)), ncol = 1)
      po.bf.matrix_quad <- matrix(rep(0, nrow(des.mat_quad)), ncol = 1)
      # pa.bf.matrix <- matrix(0, nrow = 1)
      bf.info <- cbind.data.frame(x = NA, y = NA, scale = NA, res = 1)
      basis.functions <- NULL
    }
    ############################################################################

    # TMB required data setup - depends on the type of presence-only bias handling
    dat.list <- switch(fixed.bias.type,
                       missing = list(
                         X_PO_pres = des.mat_pres,
                         X_PO_quad = des.mat_quad,
                         Z_PO_pres = po.bf.matrix_pres,
                         Z_PO_quad = po.bf.matrix_quad,
                         quad_size = data_quad[ , quad.weights.name],
                         bf_per_res = as.numeric(table(bf.info$res)),
                         approx_type = as.integer(which(approx.type == c("not_sre", "variational", "laplace")) - 1),
                         bf_type = as.integer(which(bf.matrix.type == c("sparse", "dense")) - 1),
                         fixed_bias_type = 0, # no accounting for biasing via fixed effects
                         random_bias_type = 0, # no accounting for biasing via random effects
                         model_type = 0 # PO data
                       ),
                       covariates = list(
                         X_PO_pres = des.mat_pres,
                         B_PO_pres = bias.des.mat_pres,
                         X_PO_quad = des.mat_quad,
                         B_PO_quad = bias.des.mat_quad,
                         Z_PO_pres = po.bf.matrix_pres,
                         Z_PO_quad = po.bf.matrix_quad,
                         quad_size = data_quad[ , quad.weights.name],
                         bf_per_res = as.numeric(table(bf.info$res)),
                         approx_type = as.integer(which(approx.type == c("not_sre", "variational", "laplace")) - 1),
                         bf_type = as.integer(which(bf.matrix.type == c("sparse", "dense")) - 1),
                         fixed_bias_type = 1, # accounting for biasing with provided fixed effects
                         random_bias_type = 0, # no accounting for biasing via random effects
                         model_type = 0 # PO data
                       )
    )
    ############################################################################

    ## Parameters ##############################################################

    # create the appropriate start parameters for the variance component w.r.t. approx. type
    var.starts <- switch(approx.type,
                         not_sre = rep(0, ncol(dat.list$Z_PO_pres)),
                         variational = rep(0, ncol(dat.list$Z_PO_pres)),
                         laplace = rep(0, length(dat.list$bf_per_res))
    )
    # initialise starting parameters at 0
    start.pars <- switch(fixed.bias.type,
                         missing = list(fixed = rep(0, ncol(dat.list$X_PO_pres)),
                                     random = rep(0, ncol(dat.list$Z_PO_pres)),
                                     log_variance_component = var.starts
                         ),
                         covariates = list(fixed = rep(0, ncol(dat.list$X_PO_pres)),
                                           bias = rep(0, ncol(dat.list$B_PO_pres)),
                                           random = rep(0, ncol(dat.list$Z_PO_pres)),
                                           log_variance_component = var.starts
                         )
    )

    # update to the warm starting parameters if provided
    if (!missing(starting.pars)) {
      start.pars <- update.starting.parameters(starting.pars, start.pars, target.approx.type = approx.type)
    }

    ############################################################################

    # collect the fixed effect names
    fixed.names <- colnames(des.mat_pres)
    # collect the biasing term names
    bias.names <- switch(fixed.bias.type,
                         missing = NULL,
                         covariates = colnames(bias.des.mat_pres)
    )
    random.bias.names <- NULL # cannot be present in the "PO" model case

    # get the random coefficient numbers (with <resolution level>.<coefficient number>)
    random.nos <- NULL
    if (length(dat.list$bf_per_res) == 1L) {
      random.nos <- 1L:dat.list$bf_per_res
    } else {
      for (lvl in 1L:length(dat.list$bf_per_res)) {
        random.nos <- c(random.nos, paste(lvl, 1L:dat.list$bf_per_res[lvl], sep = "."))
      }
    }
    # collect the random effect term names
    if (approx.type != "not_sre") {
      random.names <- paste0("u", random.nos)
    } else {
      random.names <- NULL
    }

    # add bias types to arg list
    arg.info$fixed.bias.type <- fixed.bias.type
    arg.info$random.bias.type <- "none"

    # collate info to be returned
    return.info <- list(tmb.data = dat.list, tmb.pars = start.pars, pt.quad.id = pt.quad.id, row.id = order(c(pres.rows, quad.rows)), fixed.names = fixed.names, bias.names = bias.names, random.names = random.names, random.bias.names = random.bias.names, bf.info = bf.info, basis.functions = basis.functions, args = arg.info)


  } else if (model.type == "PA") { # Presence/absence data

    # get the binary response
    Y = as.numeric(data[ , all.vars(formula[[2]])] > 0) # corrects in the case of abundance

    ## Fixed Effects ###########################################################

    # get the fixed effect design matrix
    des.mat <- get.design.matrix(formula, data)

    # get the offset term if present
    offset.vec <- get.offset(formula, data)

    ############################################################################

    ## Random Effects ##########################################################

    # Determine the basis functions to be used
    if (approx.type != "not_sre") {
      # when no basis functions are provided use a simple basis default
      if (missing(basis.functions)) {
        # get a rough guide for the number of basis functions (to be 50% of the presence points)
        sqrt_number_bfs <- sqrt(length(Y)*0.5)
        # set the basis function
        basis.functions <- simple_basis(sqrt_number_bfs, data, coord.names = coord.names)
      }
      # calculate the basis function matrix
      pa.bf.matrix <- get.bf.matrix(basis.functions, point.locations = data[ , coord.names], bf.matrix.type = bf.matrix.type)
      # store the basis function information
      bf.info <- attr(pa.bf.matrix, "bf.df")
    } else {
      # set a trivial example for not_sre models
      pa.bf.matrix <- matrix(rep(0, nrow(des.mat)), ncol = 1)
      # pa.bf.matrix <- matrix(0, nrow = 1)
      bf.info <- cbind.data.frame(x = NA, y = NA, scale = NA, res = 1)
      basis.functions <- NULL
    }


    ############################################################################

    # TMB required data setup
    dat.list <- list(
      X_PA = as.matrix(des.mat),
      Z_PA = pa.bf.matrix,
      Y = Y,
      OFFSET = as.vector(offset.vec),
      bf_per_res = as.numeric(table(bf.info$res)),
      approx_type = as.integer(which(approx.type == c("not_sre", "variational", "laplace")) - 1),
      bf_type = as.integer(which(bf.matrix.type == c("sparse", "dense")) - 1),
      model_type = 1,
      pa_offset = if (attr(offset.vec, "check")) {1} else {0}
    )

    ## Parameters ##############################################################

    # create the appropriate start parameters for the variance component (only laplace approach so no switch)
    var.starts <- rep(0, length(dat.list$bf_per_res))
    # initialise starting parameters at 0
    start.pars <- list(fixed = rep(0, ncol(dat.list$X_PA)),
                       random = rep(0, ncol(dat.list$Z_PA)),
                       log_variance_component = var.starts
    )

    # update to the warm starting parameters if provided
    if (!missing(starting.pars)) {
      start.pars <- update.starting.parameters(starting.pars, start.pars, target.approx.type = approx.type)
    }

    ############################################################################

    # collect the fixed effect names
    fixed.names <- colnames(des.mat)
    # collect the biasing term names
    bias.names <- NULL
    random.bias.names <- NULL

    # get the random coefficient numbers (with <resolution level>.<coefficient number>)
    random.nos <- NULL
    if (length(dat.list$bf_per_res) == 1L) {
      random.nos <- 1L:dat.list$bf_per_res
    } else {
      for (lvl in 1L:length(dat.list$bf_per_res)) {
        random.nos <- c(random.nos, paste(lvl, 1L:dat.list$bf_per_res[lvl], sep = "."))
      }
    }
    # collect the random effect term names
    if (approx.type != "not_sre") {
      random.names <- paste0("u", random.nos)
    } else {
      random.names <- NULL
    }

    # adjust the bias formula if mistakenly provided:
    if (!missing(bias.formula)) {
      # warning(paste0("'bias.formula' will be ignored since it is not compatible with model.type = 'PA'."))
      arg.info$bias.formula <- NULL
    }

    # add bias types to arg list
    arg.info$fixed.bias.type <- "missing"
    arg.info$random.bias.type <- "none"

    # collate info to be returned
    return.info <- list(tmb.data = dat.list, tmb.pars = start.pars, pt.quad.id = NA, row.id = NA, fixed.names = fixed.names, bias.names = bias.names, random.names = random.names, random.bias.names = random.bias.names, bf.info = bf.info, basis.functions = basis.functions, args = arg.info)


  } else { # Integrated data

    # get the presence point/ quadrature point identifier
    pt.quad.id <- as.numeric(data[ , all.vars(formula[[2]])])

    # track row number splits on presence and quadrature points
    row.id <- 1:nrow(data)
    pres.rows <- row.id[pt.quad.id == 1]
    quad.rows <- row.id[pt.quad.id == 0]

    # split the data into presence points and quadrature
    data_pres <- data[pt.quad.id == 1, ]
    data_quad <- data[pt.quad.id == 0, ]

    # get the binary response
    Y = as.numeric(IDM.presence.absence.df[ , all.vars(formula[[2]])] > 0) # corrects in the case of abundance

    ## Fixed Effects ###########################################################

    # get the fixed effect design matrix on the presence-only data
    des.mat <- get.design.matrix(formula, data)
    # need to adjust for a single column design matrix
    if (ncol(des.mat) == 1) {
      po.des.mat_pres <- as.matrix(data.frame(des.mat[pt.quad.id == 1, ]))
      po.des.mat_quad <- as.matrix(data.frame(des.mat[pt.quad.id == 0, ]))
      colnames(po.des.mat_pres) <- colnames(des.mat)
      colnames(po.des.mat_quad) <- colnames(des.mat)
    } else {
      # need to adjust for a single row design matrix
      if (sum(pt.quad.id == 1) == 1) {
        po.des.mat_pres <- as.matrix(as.data.frame(t(des.mat[pt.quad.id == 1, ])))
      } else {
        po.des.mat_pres <- des.mat[pt.quad.id == 1, ]
      }
      if (sum(pt.quad.id == 0) == 1) {
        po.des.mat_quad <- as.matrix(as.data.frame(t(des.mat[pt.quad.id == 0, ])))
      } else {
        po.des.mat_quad <- des.mat[pt.quad.id == 0, ]
      }
    }
    pa.des.mat <- get.design.matrix(formula, IDM.presence.absence.df)

    # get the offset term if present
    offset.vec <- get.offset(formula, IDM.presence.absence.df)

    # get the bias predictor design matrix
    if (missing(bias.formula)) {
      # assign fixed.bias.type indicator
      fixed.bias.type <- "missing"
      warning("No 'bias.formula' provided. It is strongly recommended that users include an intercept for the PO biasing process with 'bias.formula = ~ 1'")

    } else if (is(bias.formula, "formula")) {
      bias.des.mat <- get.design.matrix(bias.formula, data)
      # need to adjust for a single column design matrix
      if (ncol(bias.des.mat) == 1) {
        bias.des.mat_pres <- as.matrix(data.frame(bias.des.mat[pt.quad.id == 1, ]))
        bias.des.mat_quad <- as.matrix(data.frame(bias.des.mat[pt.quad.id == 0, ]))
        colnames(bias.des.mat_pres) <- colnames(bias.des.mat)
        colnames(bias.des.mat_quad) <- colnames(bias.des.mat)
      } else {
        # need to adjust for a single row design matrix
        if (sum(pt.quad.id == 1) == 1) {
          bias.des.mat_pres <- as.matrix(as.data.frame(t(bias.des.mat[pt.quad.id == 1, ])))
        } else {
          bias.des.mat_pres <- bias.des.mat[pt.quad.id == 1, ]
        }
        if (sum(pt.quad.id == 0) == 1) {
          bias.des.mat_quad <- as.matrix(as.data.frame(t(bias.des.mat[pt.quad.id == 0, ])))
        } else {
          bias.des.mat_quad <- bias.des.mat[pt.quad.id == 0, ]
        }
      }
      # Adjust the Intercept names if required
      if (any(grepl("(Intercept)", colnames(bias.des.mat_pres), fixed = T))) {
        colnames(bias.des.mat_pres)[grepl("(Intercept)", colnames(bias.des.mat_pres), fixed = T)] <- "(Bias Intercept)"
      }
      if (any(grepl("(Intercept)", colnames(bias.des.mat_quad), fixed = T))) {
        colnames(bias.des.mat_quad)[grepl("(Intercept)", colnames(bias.des.mat_quad), fixed = T)] <- "(Bias Intercept)"
      }
      # assign fixed.bias.type indicator
      fixed.bias.type <- "covariates"

    } else {
      stop(paste0("'bias.formula' provided is not of class formula"))
    }

    # since the model is an IDM, check whether it has spatial random effects then, account for bias using an additional latent field if needed
    if (latent.po.biasing) {
      if (missing(po.biasing.basis.functions) & !missing(basis.functions) & approx.type != "not_sre") {
        random.bias.type <- "field1"
      } else {
        random.bias.type <- "field2"
      }
    } else {
      random.bias.type <- "none"
    }
    ############################################################################

    ## Random Effects ##########################################################

    # Determine the basis functions to be used
    if (approx.type != "not_sre") {
      # when no basis functions are provided use a simple basis default
      if (missing(basis.functions)) {
        # get a rough guide for the number of basis functions (to be 25% of the presence points)
        sqrt_number_bfs <- sqrt(sum(pt.quad.id)*0.25)
        # set the basis function
        basis.functions <- simple_basis(sqrt_number_bfs, data, coord.names = coord.names)
      }
      # calculate the basis function matrices
      po.bf.matrix_pres <- get.bf.matrix(basis.functions, point.locations = data_pres[ , coord.names], bf.matrix.type = bf.matrix.type)
      po.bf.matrix_quad <- get.bf.matrix(basis.functions, point.locations = data_quad[ , coord.names], bf.matrix.type = bf.matrix.type)
      pa.bf.matrix <- get.bf.matrix(basis.functions, point.locations = IDM.presence.absence.df[ , coord.names], bf.matrix.type = bf.matrix.type)
      # store the basis function information
      bf.info <- attr(po.bf.matrix_pres, "bf.df")

      # prune the basis function if required
      if (FALSE) { # THIS IS NOT NEEDED SINCE PA DATA APPEARS TO STABLISE THE FIT
      # if (prune.bfs != 0) {
        # determine basis functions that do not intersect any presence points
        prune.id <- do.call("apply", list(po.bf.matrix_pres, MARGIN = 2, FUN = function(x){sum(x>0)})) < prune.bfs
        # check that the pruning doesn't remove all basis functions
        while (all(prune.id)) {
          prune.bfs <- prune.bfs - 1
          prune.id <- do.call("apply", list(po.bf.matrix_pres, MARGIN = 2, FUN = function(x){sum(x>0)})) < prune.bfs
          warning(paste0("'prune.bfs' = ", prune.bfs + 1, " results in no valid basis functions. Trying 'prune.bfs' = ", prune.bfs, ". Otherwise fit a model without SRE."))
        }
        # prune from both basis function matrices - adjusting for the case when only 1 bf remains
        if (sum(!prune.id) == 1) {
          po.bf.matrix_pres <- as.matrix(data.frame(po.bf.matrix_pres[ , !prune.id]))
          po.bf.matrix_quad <- as.matrix(data.frame(po.bf.matrix_quad[ , !prune.id]))
          pa.bf.matrix <- as.matrix(data.frame(pa.bf.matrix[ , !prune.id]))
          colnames(po.bf.matrix_pres) <- NULL
          colnames(po.bf.matrix_quad) <- NULL
          colnames(pa.bf.matrix) <- NULL
        } else {
          po.bf.matrix_pres <- po.bf.matrix_pres[ , !prune.id]
          po.bf.matrix_quad <- po.bf.matrix_quad[ , !prune.id]
          pa.bf.matrix <- pa.bf.matrix[ , !prune.id]
        }
        # adjust the basis function information
        tmp <- bf.info
        bf.info <- tmp[!prune.id, ]
        attr(bf.info, "pruned") <- tmp[prune.id, ]
        # adjust the supplied basis functions depending on whether they are simple_basis or FRK package
        if (is(basis.functions, "bf.df")) {
          basis.functions <- basis.functions[!prune.id, ]
        } else {
          prune.idx <- (1:length(prune.id))[prune.id]
          basis.functions <- FRK::remove_basis(basis.functions, prund.idx)
        }
      }
    } else {
      # set a trivial example for not_sre models
      po.bf.matrix_pres <- matrix(rep(0, nrow(po.des.mat_pres)), ncol = 1)
      po.bf.matrix_quad <- matrix(rep(0, nrow(po.des.mat_quad)), ncol = 1)
      pa.bf.matrix <- matrix(rep(0, nrow(pa.des.mat)), ncol = 1)
      # adjust for sparsity if needed
      if (bf.matrix.type == "sparse") {
        po.bf.matrix_pres <- methods::as(po.bf.matrix_pres, "sparseMatrix")
        po.bf.matrix_quad <- methods::as(po.bf.matrix_quad, "sparseMatrix")
        pa.bf.matrix  <- methods::as(pa.bf.matrix , "sparseMatrix")
      }
      bf.info <- cbind.data.frame(x = NA, y = NA, scale = NA, res = 1)
      basis.functions <- NULL
    }

    # Determine the additional basis functions to be used for presence-only biasing
    if (random.bias.type == "field2") {
      # when no basis functions are provided use a simple basis default
      if (missing(po.biasing.basis.functions)) {
        # get a rough guide for the number of basis functions (to be 25% of the presence points)
        sqrt_number_bfs <- sqrt(sum(pt.quad.id)*0.25)
        # set the basis functions
        po.biasing.basis.functions <- simple_basis(sqrt_number_bfs, data, coord.names = coord.names)
      }
      # calculate the basis function matrices
      po.bias.bf.matrix_pres <- get.bf.matrix(po.biasing.basis.functions, point.locations = data_pres[ , coord.names], bf.matrix.type = bf.matrix.type)
      po.bias.bf.matrix_quad <- get.bf.matrix(po.biasing.basis.functions, point.locations = data_quad[ , coord.names], bf.matrix.type = bf.matrix.type)
      # store the basis function information
      bias.bf.info <- attr(po.bias.bf.matrix_pres, "bf.df")

      # prune the basis functions if required
      if (prune.bfs != 0) {
        # determine basis functions that do not intersect any presence points
        bias.prune.id <- do.call("apply", list(po.bias.bf.matrix_pres, MARGIN = 2, FUN = function(x){sum(x>0)})) < prune.bfs
        # check that the pruning doesn't remove all basis functions
        while (all(bias.prune.id)) {
          prune.bfs <- prune.bfs - 1
          bias.prune.id <- do.call("apply", list(po.bias.bf.matrix_pres, MARGIN = 2, FUN = function(x){sum(x>0)})) < prune.bfs
          warning(paste0("'prune.bfs' = ", prune.bfs + 1, " results in no valid basis functions. Trying 'prune.bfs' = ", prune.bfs, ". Otherwise fit a model without SRE."))
        }
        # prune from both basis function matrices - adjusting for the case when only 1 bf remains
        if (sum(!bias.prune.id) == 1) {
          po.bias.bf.matrix_pres <- as.matrix(data.frame(po.bias.bf.matrix_pres[ , !bias.prune.id]))
          po.bias.bf.matrix_quad <- as.matrix(data.frame(po.bias.bf.matrix_quad[ , !bias.prune.id]))
          colnames(po.bias.bf.matrix_pres) <- NULL
          colnames(po.bias.bf.matrix_quad) <- NULL
          # re-adjust to sparse matrix if required
          if (bf.matrix.type == "sparse") {
            po.bias.bf.matrix_pres <- methods::as(po.bias.bf.matrix_pres, "sparseMatrix")
            po.bias.bf.matrix_quad <- methods::as(po.bias.bf.matrix_quad, "sparseMatrix")
          }
        } else {
          po.bias.bf.matrix_pres <- po.bias.bf.matrix_pres[ , !bias.prune.id]
          po.bias.bf.matrix_quad <- po.bias.bf.matrix_quad[ , !bias.prune.id]
        }
        # adjust the basis function information
        tmp <- bias.bf.info
        bias.bf.info <- tmp[!bias.prune.id, ]
        attr(bias.bf.info, "pruned") <- tmp[bias.prune.id, ]
        # adjust the supplied basis functions depending on whether they are simple_basis or FRK package
        if (is(po.biasing.basis.functions, "bf.df")) {
          po.biasing.basis.functions <- po.biasing.basis.functions[!bias.prune.id, ]
        } else {
          bias.prune.idx <- (1:length(bias.prune.id))[bias.prune.id]
          po.biasing.basis.functions <- FRK::remove_basis(po.biasing.basis.functions, bias.prund.idx)
        }
      }
      # if the model is not a SRE, adjust the approx.type. This will fit the trivial basis functions for the shared field
      if (approx.type == "not_sre") {
        approx.type <- "laplace"
      }
    } else {
      bias.bf.info <- NULL
      po.biasing.basis.functions <- NULL
    }
    ############################################################################

    # TMB required data setup
    dat.list <- switch(fixed.bias.type,
                       missing = switch(random.bias.type,
                                        none = list(
                                          X_PO_pres = po.des.mat_pres,
                                          X_PO_quad = po.des.mat_quad,
                                          X_PA = pa.des.mat,
                                          Z_PO_pres = po.bf.matrix_pres,
                                          Z_PO_quad = po.bf.matrix_quad,
                                          Z_PA = pa.bf.matrix,
                                          quad_size = data_quad[ , quad.weights.name],
                                          Y = Y,
                                          OFFSET = as.vector(offset.vec),
                                          bf_per_res = as.numeric(table(bf.info$res)),
                                          approx_type = as.integer(which(approx.type == c("not_sre", "variational", "laplace")) - 1),
                                          bf_type = as.integer(which(bf.matrix.type == c("sparse", "dense")) - 1),
                                          fixed_bias_type = 0, # no accounting for biasing via fixed effects
                                          random_bias_type = 0, # no accounting for biasing via random effects
                                          model_type = 2, # integrated data
                                          pa_offset = if (attr(offset.vec, "check")) {1} else {0}
                                        ),
                                        field1 = list(
                                          X_PO_pres = po.des.mat_pres,
                                          X_PO_quad = po.des.mat_quad,
                                          X_PA = pa.des.mat,
                                          Z_PO_pres = po.bf.matrix_pres,
                                          Z_PO_quad = po.bf.matrix_quad,
                                          Z_PA = pa.bf.matrix,
                                          quad_size = data_quad[ , quad.weights.name],
                                          Y = Y,
                                          OFFSET = as.vector(offset.vec),
                                          bf_per_res = as.numeric(table(bf.info$res)),
                                          approx_type = as.integer(which(approx.type == c("not_sre", "variational", "laplace")) - 1),
                                          bf_type = as.integer(which(bf.matrix.type == c("sparse", "dense")) - 1),
                                          fixed_bias_type = 0, # no accounting for biasing via fixed effects
                                          random_bias_type = 1, # accounting for biasing using additional latent field on existing basis functions
                                          model_type = 2, # integrated data
                                          pa_offset = if (attr(offset.vec, "check")) {1} else {0}
                                        ),
                                        field2 = list(
                                          X_PO_pres = po.des.mat_pres,
                                          X_PO_quad = po.des.mat_quad,
                                          X_PA = pa.des.mat,
                                          Z_PO_pres = po.bf.matrix_pres,
                                          Z_PO_quad = po.bf.matrix_quad,
                                          Z2_PO_pres = po.bias.bf.matrix_pres,
                                          Z2_PO_quad = po.bias.bf.matrix_quad,
                                          Z_PA = pa.bf.matrix,
                                          quad_size = data_quad[ , quad.weights.name],
                                          Y = Y,
                                          OFFSET = as.vector(offset.vec),
                                          bf_per_res = as.numeric(table(bf.info$res)),
                                          bias_bf_per_res = as.numeric(table(bias.bf.info$res)),
                                          approx_type = as.integer(which(approx.type == c("not_sre", "variational", "laplace")) - 1),
                                          bf_type = as.integer(which(bf.matrix.type == c("sparse", "dense")) - 1),
                                          fixed_bias_type = 0, # no accounting for biasing via fixed effects
                                          random_bias_type = 2, # accounting for biasing using additional latent field on additional basis functions
                                          model_type = 2, # integrated data
                                          pa_offset = if (attr(offset.vec, "check")) {1} else {0}
                                        )
                       ),
                       covariates = switch(random.bias.type,
                                           none = list(
                                             X_PO_pres = po.des.mat_pres,
                                             B_PO_pres = bias.des.mat_pres,
                                             X_PO_quad = po.des.mat_quad,
                                             B_PO_quad = bias.des.mat_quad,
                                             X_PA = pa.des.mat,
                                             Z_PO_pres = po.bf.matrix_pres,
                                             Z_PO_quad = po.bf.matrix_quad,
                                             Z_PA = pa.bf.matrix,
                                             quad_size = data_quad[ , quad.weights.name],
                                             Y = Y,
                                             OFFSET = as.vector(offset.vec),
                                             bf_per_res = as.numeric(table(bf.info$res)),
                                             approx_type = as.integer(which(approx.type == c("not_sre", "variational", "laplace")) - 1),
                                             bf_type = as.integer(which(bf.matrix.type == c("sparse", "dense")) - 1),
                                             fixed_bias_type = 1, # accounting for biasing with fixed effects
                                             random_bias_type = 0, # no accounting for biasing via random effects
                                             model_type = 2, # integrated data
                                             pa_offset = if (attr(offset.vec, "check")) {1} else {0}
                                           ),
                                           field1 = list(
                                             X_PO_pres = po.des.mat_pres,
                                             B_PO_pres = bias.des.mat_pres,
                                             X_PO_quad = po.des.mat_quad,
                                             B_PO_quad = bias.des.mat_quad,
                                             X_PA = pa.des.mat,
                                             Z_PO_pres = po.bf.matrix_pres,
                                             Z_PO_quad = po.bf.matrix_quad,
                                             Z_PA = pa.bf.matrix,
                                             quad_size = data_quad[ , quad.weights.name],
                                             Y = Y,
                                             OFFSET = as.vector(offset.vec),
                                             bf_per_res = as.numeric(table(bf.info$res)),
                                             approx_type = as.integer(which(approx.type == c("not_sre", "variational", "laplace")) - 1),
                                             bf_type = as.integer(which(bf.matrix.type == c("sparse", "dense")) - 1),
                                             fixed_bias_type = 1, # accounting for biasing with fixed effects
                                             random_bias_type = 1, # accounting for biasing using additional latent field on existing basis functions
                                             model_type = 2, # integrated data
                                             pa_offset = if (attr(offset.vec, "check")) {1} else {0}
                                           ),
                                           field2 = list(
                                             X_PO_pres = po.des.mat_pres,
                                             B_PO_pres = bias.des.mat_pres,
                                             X_PO_quad = po.des.mat_quad,
                                             B_PO_quad = bias.des.mat_quad,
                                             X_PA = pa.des.mat,
                                             Z_PO_pres = po.bf.matrix_pres,
                                             Z_PO_quad = po.bf.matrix_quad,
                                             Z2_PO_pres = po.bias.bf.matrix_pres,
                                             Z2_PO_quad = po.bias.bf.matrix_quad,
                                             Z_PA = pa.bf.matrix,
                                             quad_size = data_quad[ , quad.weights.name],
                                             Y = Y,
                                             OFFSET = as.vector(offset.vec),
                                             bf_per_res = as.numeric(table(bf.info$res)),
                                             bias_bf_per_res = as.numeric(table(bias.bf.info$res)),
                                             approx_type = as.integer(which(approx.type == c("not_sre", "variational", "laplace")) - 1),
                                             bf_type = as.integer(which(bf.matrix.type == c("sparse", "dense")) - 1),
                                             fixed_bias_type = 1, # accounting for biasing with fixed effects
                                             random_bias_type = 2, # accounting for biasing using additional latent field on additional basis functions
                                             model_type = 2, # integrated data
                                             pa_offset = if (attr(offset.vec, "check")) {1} else {0}
                                           )
                       )
    )

    ## Parameters ##############################################################

    # initialise starting parameters at 0 (can only use Laplace approach so no var.starts switch needed)
    start.pars <- switch(fixed.bias.type,
                         missing = switch(random.bias.type,
                                          none = list(fixed = rep(0, ncol(dat.list$X_PO_pres)),
                                                      random = rep(0, ncol(dat.list$Z_PO_pres)),
                                                      log_variance_component = rep(0, length(dat.list$bf_per_res))
                                          ),
                                          field1 = list(fixed = rep(0, ncol(dat.list$X_PO_pres)),
                                                        random = rep(0, ncol(dat.list$Z_PO_pres)),
                                                        random_bias = rep(0, sum(dat.list$bf_per_res)),
                                                        log_variance_component = rep(0, length(dat.list$bf_per_res)),
                                                        log_variance_component_bias = rep(0, length(dat.list$bf_per_res))
                                          ),
                                          field2 = list(fixed = rep(0, ncol(dat.list$X_PO_pres)),
                                                        random = rep(0, ncol(dat.list$Z_PO_pres)),
                                                        random_bias = rep(0, sum(dat.list$bias_bf_per_res)),
                                                        log_variance_component = rep(0, length(dat.list$bf_per_res)),
                                                        log_variance_component_bias = rep(0, length(dat.list$bias_bf_per_res))
                                          )
                         ),
                         covariates = switch(random.bias.type,
                                             none = list(fixed = rep(0, ncol(dat.list$X_PO_pres)),
                                                         bias = rep(0, ncol(dat.list$B_PO_pres)),
                                                         random = rep(0, ncol(dat.list$Z_PO_pres)),
                                                         log_variance_component = rep(0, length(dat.list$bf_per_res))
                                             ),
                                             field1 = list(fixed = rep(0, ncol(dat.list$X_PO_pres)),
                                                           bias = rep(0, ncol(dat.list$B_PO_pres)),
                                                           random = rep(0, ncol(dat.list$Z_PO_pres)),
                                                           random_bias = rep(0, sum(dat.list$bf_per_res)),
                                                           log_variance_component = rep(0, length(dat.list$bf_per_res)),
                                                           log_variance_component_bias = rep(0, length(dat.list$bf_per_res))
                                             ),
                                             field2 = list(fixed = rep(0, ncol(dat.list$X_PO_pres)),
                                                           bias = rep(0, ncol(dat.list$B_PO_pres)),
                                                           random = rep(0, ncol(dat.list$Z_PO_pres)),
                                                           random_bias = rep(0, sum(dat.list$bias_bf_per_res)),
                                                           log_variance_component = rep(0, length(dat.list$bf_per_res)),
                                                           log_variance_component_bias = rep(0, length(dat.list$bias_bf_per_res))
                                             )
                         )
    )

    # update to the warm starting parameters if provided
    if (!missing(starting.pars)) {
      start.pars <- update.starting.parameters(starting.pars, start.pars, target.approx.type = approx.type)
    }

    ############################################################################

    # collect the fixed effect names
    fixed.names <- colnames(po.des.mat_pres)

    # get the random coefficient numbers (with <resolution level>.<coefficient number>)
    random.nos <- NULL
    if (length(dat.list$bf_per_res) == 1L) {
      random.nos <- 1L:dat.list$bf_per_res
    } else {
      for (lvl in 1L:length(dat.list$bf_per_res)) {
        random.nos <- c(random.nos, paste(lvl, 1L:dat.list$bf_per_res[lvl], sep = "."))
      }
    }
    # check if there is a bias field in an IDM model
    bias.random.nos <- NULL
    if (random.bias.type %in% c("field1", "field2")) {
      if (!is.null(dat.list$bias_bf_per_res)) {
        if (length(dat.list$bias_bf_per_res) == 1L) {
          bias.random.nos <- 1L:dat.list$bias_bf_per_res
        } else {
          for (lvl in 1L:length(dat.list$bias_bf_per_res)) {
            bias.random.nos <- c(bias.random.nos, paste(lvl, 1L:dat.list$bias_bf_per_res[lvl], sep = "."))
          }
        }
      } else {
        bias.random.nos <- random.nos
      }
    }
    # collect the fixed effect biasing term names
    bias.names <- switch(fixed.bias.type,
                               missing = NULL,
                               covariates = colnames(bias.des.mat_pres)
    )
    # collect the random effect biasing term names
    if(random.bias.type %in% c("field1", "field2")) {
      random.bias.names <- paste0("tau", bias.random.nos)
    } else {
      random.bias.names <- NULL
    }

    # collect the random effect term names
    if (approx.type != "not_sre") {
      random.names <- paste0("u", random.nos)
    } else {
      random.names <- NULL
    }

    # add bias types to arg list
    arg.info$fixed.bias.type <- fixed.bias.type
    arg.info$random.bias.type <- random.bias.type

    # adjust the approx.type within the arg.info list (in case this has been changed in the instance where there are no shared SRE but PO biasing SRE)
    arg.info$approx.type <- approx.type

    # collate info to be returned
    return.info <- list(tmb.data = dat.list, tmb.pars = start.pars, pt.quad.id = pt.quad.id, row.id = order(c(pres.rows, quad.rows)), fixed.names = fixed.names, bias.names = bias.names, random.names = random.names, random.bias.names = random.bias.names, bf.info = bf.info, bias.bf.info = bias.bf.info, basis.functions = basis.functions, po.biasing.basis.functions = po.biasing.basis.functions, args = arg.info)

  }
}
