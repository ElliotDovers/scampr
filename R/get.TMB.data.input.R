#' Internal scampr function that creates a list of data and starting parameters for scampr models
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the fixed effects of the model to be fitted. The 'response' must be a binary that indicates whether a datum is a presence or: quadrature point (for point process models)/ absence (for binary models). See GLM function for further formula details.
#' @param data a data frame containing response and predictors within \code{formula}.
#' @param bias.formula an object of class "formula" (or one that can be coerced to that class) OR the character string "latent". In the formula case, this is a symbolic description of the predictors included to account for bias in the presence-only data (no response term is needed). In the case of fitting an integrated data model, \code{bias.formula = "latent"} will fit an approximate latent Gaussian field to account for the bias.
#' @param IDM.presence.absence.df an optional data frame. When fitting an integrated data model use this to pass in the presence/absence data.
#' @param coord.names a vector of character strings describing the column names of the coordinates in both data frames.
#' @param quad.weights.name a charater string of the column name of quadrature weights in the data.
#' @param model.type a character string indicating the type of model to be used. May be one of 'laplace' or 'variational' for Cox Processes involving spatially correlated errors or 'ipp' for a model that follows an inhomgeneous Poisson process.
#' @param data.type a character string indicating the type of data to be used. May be one of 'PO' (for a presence-only PPM) or 'PA' (for a presence/absence Binary GLM) or 'IDM' (for an integrated data model).
#' @param basis.functions an optional object of class 'Basis' created by \code{FRK::auto_basis()} or 'bf.df' created by \code{scampr::simple_basis()}. Either object describes a set of basis functions for approximating the latent Gaussian field. If NULL the model will use default \code{FRK::auto_basis()} with \code{max_basis = 0.25 * # of points}.
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to use sparse or dense matrix computations for the basis functions created.
#' @param starting.pars an optional named list or scampr model object that gives warm starting values for the parameters of the model.
#' @param po.biasing.basis.functions an optional extra set of basis functions that can be used when \code{bias.formula = "latent"}, otherwise \code{basis.functions} are used.
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
#' tmb.input <- scampr:::get.TMB.data.input(pres ~ MNT + D.Main, sp1 ~ MNT, po.data = dat_po, pa.data = dat_pa, model.type = "ipp")
#' str(tmp.input)
get.TMB.data.input <- function(formula, data, bias.formula, IDM.presence.absence.df, coord.names = c("x", "y"), quad.weights.name = "quad.size", model.type = c("variational", "laplace", "ipp"), data.type = c("PO", "PA", "IDM"), basis.functions, bf.matrix.type = c("sparse", "dense"), starting.pars, po.biasing.basis.functions) {

  # checks for parameters of restricted strings
  data.type <- match.arg(data.type)
  model.type <- match.arg(model.type)
  bf.matrix.type <- match.arg(bf.matrix.type)

  # store the arguments
  arg.info <- as.list(environment())

  # approach based on the data type #

  if (data.type == "PO") { # Presence-only data

    # get the presence point/ quadrature point identifier
    pt.quad.id <- as.numeric(data[ , all.vars(formula[[2]])])

    # track row number splits on presence and quadrature points
    row.id <- 1:nrow(data)
    pres.rows <- row.id[pt.quad.id == 1]
    quad.rows <- row.id[pt.quad.id == 0]

    ## Fixed Effects ###########################################################

    # get the fixed effect design matrix
    des.mat <- get.design.matrix(formula, data)

    # get the bias predictor design matrix
    if (missing(bias.formula)) {
      bias.type <- "none"
    } else if (is(bias.formula, "formula")) {
      bias.des.mat <- get.design.matrix(bias.formula, data)
      # Adjust the Intercept name if required
      if (any(grepl("(Intercept)", colnames(bias.des.mat), fixed = T))) {
        colnames(bias.des.mat)[grepl("(Intercept)", colnames(bias.des.mat), fixed = T)] <- "(Bias Intercept)"
      }
      bias.type <- "covariates"
    } else {
      stop(paste0("The type of 'bias.formula' provided is not compatible with 'data.type' = ", data.type))
    }
    ############################################################################

    ## Random Effects ##########################################################

    # Determine the basis functions to be used
    if (model.type != "ipp") {
      # when no basis functions are provided use a simple basis default
      if (missing(basis.functions)) {
        # get a rough guide for the number of basis functions (to be 50% of the presence points)
        sqrt_number_bfs <- sqrt(sum(pt.quad.id)*0.5)
        # set the basis function
        basis.funtions <- simple_basis(sqrt_number_bfs, data, coord.names = coord.names)
      }
      # calculate the basis function matrix
      po.bf.matrix <- get.bf.matrix(basis.functions, point.locations = data[ , coord.names], bf.matrix.type = bf.matrix.type)
      # store the basis function information
      bf.info <- attr(po.bf.matrix, "bf.df")
    } else {
      # set a trivial example for IPP models
      po.bf.matrix <- matrix(rep(0, nrow(des.mat)), ncol = 1)
      # pa.bf.matrix <- matrix(0, nrow = 1)
      bf.info <- cbind.data.frame(x = cNA, y = NA, scale = NA, res = 1)
      basis.functions <- NULL
    }
    ############################################################################

    # TMB required data setup - depends on the type of presence-only bias handling
    dat.list <- switch(bias.type,
                       none = list(
                         X_PO_pres = as.matrix(des.mat[pt.quad.id == 1, ]),
                         X_PO_quad = as.matrix(des.mat[pt.quad.id == 0, ]),
                         Z_PO_pres = po.bf.matrix[pt.quad.id == 1, ],
                         Z_PO_quad = po.bf.matrix[pt.quad.id == 0, ],
                         quad_size = data[pt.quad.id == 0, quad.weights.name],
                         bf_per_res = as.numeric(table(bf.info$res)),
                         mod_type = as.integer(which(model.type == c("ipp", "variational", "laplace")) - 1),
                         bf_type = as.integer(which(bf.matrix.type == c("sparse", "dense")) - 1),
                         bias_type = 0, # no accounting for biasing (outside of fixed effects)
                         data_type = 0 # PO data
                       ),
                       covariates = list(
                         X_PO_pres = as.matrix(des.mat[pt.quad.id == 1, ]),
                         B_PO_pres = as.matrix(bias.des.mat[pt.quad.id == 1, ]),
                         X_PO_quad = as.matrix(des.mat[pt.quad.id == 0, ]),
                         B_PO_quad = as.matrix(bias.des.mat[pt.quad.id == 0, ]),
                         Z_PO_pres = po.bf.matrix[pt.quad.id == 1, ],
                         Z_PO_quad = po.bf.matrix[pt.quad.id == 0, ],
                         quad_size = data[pt.quad.id == 0, quad.weights.name],
                         bf_per_res = as.numeric(table(bf.info$res)),
                         mod_type = as.integer(which(model.type == c("ipp", "variational", "laplace")) - 1),
                         bf_type = as.integer(which(bf.matrix.type == c("sparse", "dense")) - 1),
                         bias_type = 1, # accounting for biasing with provided covariates
                         data_type = 0 # PO data
                       )
    )
    ############################################################################

    ## Parameters ##############################################################

    # create the appropriate start parameters for the variance component w.r.t. approx. type
    var.starts <- switch(model.type,
                         ipp = rep(0, ncol(dat.list$Z_PO_pres)),
                         variational = rep(0, ncol(dat.list$Z_PO_pres)),
                         laplace = rep(0, length(dat.list$bf_per_res))
    )
    # initialise starting parameters at 0
    start.pars <- switch(bias.type,
                       none = list(fixed = rep(0, ncol(dat.list$X_PO_pres)),
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
      start.pars <- update.starting.parameters(starting.pars, start.pars, target.model.type = model.type)
    }

    ############################################################################

    # collect the fixed effect names
    fixed.names <- colnames(des.mat)
    # collect the biasing term names
    bias.names <- NULL
    if (bias.type == "covariates") {
      bias.names <-colnames(bias.des.mat)
    }
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
    if (model.type != "ipp") {
      random.names <- paste0("u", random.nos)
    } else {
      random.names <- NULL
    }

    # add bias.type to arg list
    arg.info$bias.type <- bias.type

    # collate info to be returned
    return.info <- list(tmb.data = dat.list, tmb.pars = start.pars, pt.quad.id = pt.quad.id, row.id = order(c(pres.rows, quad.rows)), fixed.names = fixed.names, bias.names = bias.names, random.names = random.names, bf.info = bf.info, basis.functions = basis.functions, args = arg.info)


  } else if (data.type == "PA") { # Presence/absence data

    # get the binary response
    Y = as.numeric(data[ , all.vars(formula[[2]])] > 0) # corrects in the case of abundance

    ## Fixed Effects ###########################################################

    # get the fixed effect design matrix
    des.mat <- get.design.matrix(formula, data)

    ############################################################################

    ## Random Effects ##########################################################

    # Determine the basis functions to be used
    if (model.type != "ipp") {
      # when no basis functions are provided use a simple basis default
      if (missing(basis.functions)) {
        # get a rough guide for the number of basis functions (to be 50% of the presence points)
        sqrt_number_bfs <- sqrt(length(Y)*0.5)
        # set the basis function
        basis.funtions <- simple_basis(sqrt_number_bfs, data, coord.names = coord.names)
      }
      # calculate the basis function matrix
      pa.bf.matrix <- get.bf.matrix(basis.functions, point.locations = data[ , coord.names], bf.matrix.type = bf.matrix.type)
      # store the basis function information
      bf.info <- attr(pa.bf.matrix, "bf.df")
    } else {
      # set a trivial example for IPP models
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
      bf_per_res = as.numeric(table(bf.info$res)),
      mod_type = as.integer(which(model.type == c("ipp", "variational", "laplace")) - 1),
      bf_type = as.integer(which(bf.matrix.type == c("sparse", "dense")) - 1),
      data_type = 1
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
      start.pars <- update.starting.parameters(starting.pars, start.pars, target.model.type = model.type)
    }

    ############################################################################

    # collect the fixed effect names
    fixed.names <- colnames(des.mat)
    # collect the biasing term names
    bias.names <- NULL

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
    if (model.type != "ipp") {
      random.names <- paste0("u", random.nos)
    } else {
      random.names <- NULL
    }

    # add bias.type to arg list
    arg.info$bias.type <- "none"

    # collate info to be returned
    return.info <- list(tmb.data = dat.list, tmb.pars = start.pars, pt.quad.id = NA, row.id = NA, fixed.names = fixed.names, bias.names = bias.names, random.names = random.names, bf.info = bf.info, basis.functions = basis.functions, args = arg.info)


  } else { # Integrated data

    # get the presence point/ quadrature point identifier
    pt.quad.id <- as.numeric(data[ , all.vars(formula[[2]])])

    # track row number splits on presence and quadrature points
    row.id <- 1:nrow(data)
    pres.rows <- row.id[pt.quad.id == 1]
    quad.rows <- row.id[pt.quad.id == 0]

    # get the binary response
    Y = as.numeric(IDM.presence.absence.df[ , all.vars(formula[[2]])] > 0) # corrects in the case of abundance

    ## Fixed Effects ###########################################################

    # get the fixed effect design matrix on the presence-only data
    po.des.mat <- get.design.matrix(formula, data)
    pa.des.mat <- get.design.matrix(formula, IDM.presence.absence.df)

    # get the bias predictor design matrix
    if (missing(bias.formula)) {
      # assign bias.type indicator
      bias.type <- "none"
    } else if (is(bias.formula, "formula")) {
      bias.des.mat <- get.design.matrix(bias.formula, data)
      # Adjust the Intercept name if required
      if (any(grepl("(Intercept)", colnames(bias.des.mat), fixed = T))) {
        colnames(bias.des.mat)[grepl("(Intercept)", colnames(bias.des.mat), fixed = T)] <- "(Bias Intercept)"
      }
      # assign bias.type indicator
      bias.type <- "covariates"
    } else if (bias.formula == "latent") {
      if (missing(po.biasing.basis.functions)) {
        # assign bias.type indicator
        bias.type <- "latent"
      } else {
        # assign bias.type indicator
        bias.type <- "new_latent"
      }
    } else {
      stop(paste0("The type of 'bias.formula' provided is not compatible with 'data.type' = ", data.type))
    }
    ############################################################################

    ## Random Effects ##########################################################

    # Determine the basis functions to be used
    if (model.type != "ipp") {
      # when no basis functions are provided use a simple basis default
      if (missing(basis.functions)) {
        # get a rough guide for the number of basis functions (to be 50% of the presence points)
        sqrt_number_bfs <- sqrt(sum(pt.quad.id)*0.5)
        # set the basis function
        basis.funtions <- simple_basis(sqrt_number_bfs, data, coord.names = coord.names)
      }
      # calculate the basis function matrices
      po.bf.matrix <- get.bf.matrix(basis.functions, point.locations = data[ , coord.names], bf.matrix.type = bf.matrix.type)
      pa.bf.matrix <- get.bf.matrix(basis.functions, point.locations = IDM.presence.absence.df[ , coord.names], bf.matrix.type = bf.matrix.type)
      # store the basis function information
      bf.info <- attr(po.bf.matrix, "bf.df")
    } else {
      # set a trivial example for IPP models
      po.bf.matrix <- matrix(rep(0, nrow(po.des.mat)), ncol = 1)
      pa.bf.matrix <- matrix(rep(0, nrow(pa.des.mat)), ncol = 1)
      # pa.bf.matrix <- matrix(0, nrow = 1)
      bf.info <- cbind.data.frame(x = NA, y = NA, scale = NA, res = 1)
      basis.functions <- NULL
    }
    # Determine the additional basis functions to be used for presence
    if (bias.type == "new_latent") {
      # calculate the basis function matrices
      po.bias.bf.matrix <- get.bf.matrix(po.biasing.basis.functions, point.locations = data[ , coord.names], bf.matrix.type = bf.matrix.type)
      # store the basis function information
      bias.bf.info <- attr(po.bias.bf.matrix, "bf.df")
    } else {
      bias.bf.info <- NULL
    }
    ############################################################################

    # TMB required data setup
    dat.list <- switch(bias.type,
                       none = list(
                         X_PO_pres = as.matrix(po.des.mat[pt.quad.id == 1, ]),
                         X_PO_quad = as.matrix(po.des.mat[pt.quad.id == 0, ]),
                         X_PA = as.matrix(pa.des.mat),
                         Z_PO_pres = po.bf.matrix[pt.quad.id == 1, ],
                         Z_PO_quad = po.bf.matrix[pt.quad.id == 0, ],
                         Z_PA = pa.bf.matrix,
                         quad_size = data[pt.quad.id == 0, quad.weights.name],
                         Y = Y,
                         bf_per_res = as.numeric(table(bf.info$res)),
                         mod_type = as.integer(which(model.type == c("ipp", "variational", "laplace")) - 1),
                         bf_type = as.integer(which(bf.matrix.type == c("sparse", "dense")) - 1),
                         bias_type = 0, # no accounting for biasing (outside of fixed effects)
                         data_type = 2 # integrated data
                       ),
                       covariates = list(
                         X_PO_pres = as.matrix(po.des.mat[pt.quad.id == 1, ]),
                         B_PO_pres = as.matrix(bias.des.mat[pt.quad.id == 1, ]),
                         X_PO_quad = as.matrix(po.des.mat[pt.quad.id == 0, ]),
                         B_PO_quad = as.matrix(bias.des.mat[pt.quad.id == 0, ]),
                         X_PA = as.matrix(pa.des.mat),
                         Z_PO_pres = po.bf.matrix[pt.quad.id == 1, ],
                         Z_PO_quad = po.bf.matrix[pt.quad.id == 0, ],
                         Z_PA = pa.bf.matrix,
                         quad_size = data[pt.quad.id == 0, quad.weights.name],
                         Y = Y,
                         bf_per_res = as.numeric(table(bf.info$res)),
                         mod_type = as.integer(which(model.type == c("ipp", "variational", "laplace")) - 1),
                         bf_type = as.integer(which(bf.matrix.type == c("sparse", "dense")) - 1),
                         bias_type = 1, # accounting for biasing using covariates
                         data_type = 2 # integrated data
                       ),
                       latent = list(
                         X_PO_pres = as.matrix(po.des.mat[pt.quad.id == 1, ]),
                         X_PO_quad = as.matrix(po.des.mat[pt.quad.id == 0, ]),
                         X_PA = as.matrix(pa.des.mat),
                         Z_PO_pres = po.bf.matrix[pt.quad.id == 1, ],
                         Z_PO_quad = po.bf.matrix[pt.quad.id == 0, ],
                         Z_PA = pa.bf.matrix,
                         quad_size = data[pt.quad.id == 0, quad.weights.name],
                         Y = Y,
                         bf_per_res = as.numeric(table(bf.info$res)),
                         mod_type = as.integer(which(model.type == c("ipp", "variational", "laplace")) - 1),
                         bf_type = as.integer(which(bf.matrix.type == c("sparse", "dense")) - 1),
                         bias_type = 2, # accounting for biasing using additional latent field on existing basis functions
                         data_type = 2 # integrated data
                       ),
                       new_latent = list(
                         X_PO_pres = as.matrix(po.des.mat[pt.quad.id == 1, ]),
                         X_PO_quad = as.matrix(po.des.mat[pt.quad.id == 0, ]),
                         X_PA = as.matrix(pa.des.mat),
                         Z_PO_pres = po.bf.matrix[pt.quad.id == 1, ],
                         Z_PO_quad = po.bf.matrix[pt.quad.id == 0, ],
                         Z2_PO_pres = po.bias.bf.matrix[pt.quad.id == 1, ],
                         Z2_PO_quad = po.bias.bf.matrix[pt.quad.id == 0, ],
                         Z_PA = pa.bf.matrix,
                         quad_size = data[pt.quad.id == 0, quad.weights.name],
                         Y = Y,
                         bf_per_res = as.numeric(table(bf.info$res)),
                         bias_bf_per_res = as.numeric(table(bias.bf.info$res)),
                         mod_type = as.integer(which(model.type == c("ipp", "variational", "laplace")) - 1),
                         bf_type = as.integer(which(bf.matrix.type == c("sparse", "dense")) - 1),
                         bias_type = 3, # accounting for biasing using additional latent field on additional basis functions
                         data_type = 2 # integrated data
                       )
    )

    ## Parameters ##############################################################

    # initialise starting parameters at 0 (can only use Laplace approach so no var.starts switch needed)
    start.pars <- switch(bias.type,
                         none = list(fixed = rep(0, ncol(dat.list$X_PO_pres)),
                                     random = rep(0, ncol(dat.list$Z_PO_pres)),
                                     log_variance_component = rep(0, length(dat.list$bf_per_res))
                         ),
                         covariates = list(fixed = rep(0, ncol(dat.list$X_PO_pres)),
                                           bias = rep(0, ncol(dat.list$B_PO_pres)),
                                           random = rep(0, ncol(dat.list$Z_PO_pres)),
                                           log_variance_component = rep(0, length(dat.list$bf_per_res))
                         ),
                         latent = list(fixed = rep(0, ncol(dat.list$X_PO_pres)),
                                     random = rep(0, ncol(dat.list$Z_PO_pres)),
                                     bias = rep(0, sum(dat.list$bf_per_res)),
                                     log_variance_component = rep(0, length(dat.list$bf_per_res)),
                                     log_variance_component_bias = rep(0, length(dat.list$bf_per_res))
                         ),
                         new_latent = list(fixed = rep(0, ncol(dat.list$X_PO_pres)),
                                       random = rep(0, ncol(dat.list$Z_PO_pres)),
                                       bias = rep(0, sum(dat.list$bias_bf_per_res)),
                                       log_variance_component = rep(0, length(dat.list$bf_per_res)),
                                       log_variance_component_bias = rep(0, length(dat.list$bias_bf_per_res))
                         )
    )

    # update to the warm starting parameters if provided
    if (!missing(starting.pars)) {
      start.pars <- update.starting.parameters(starting.pars, start.pars, target.model.type = model.type)
    }

    ############################################################################

    # collect the fixed effect names
    fixed.names <- colnames(po.des.mat)

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
    if (bias.type %in% c("latent", "new_latent")) {
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
    # collect the biasing term names
    bias.names <- switch(bias.type,
                         none = NULL,
                         covariates = colnames(bias.des.mat),
                         latent = paste0("tau", bias.random.nos),
                         new_latent = paste0("tau", bias.random.nos)
    )
    # collect the random effect term names
    if (model.type != "ipp") {
      random.names <- paste0("u", random.nos)
    } else {
      random.names <- NULL
    }

    # add bias.type to arg list
    arg.info$bias.type <- bias.type

    # collate info to be returned
    return.info <- list(tmb.data = dat.list, tmb.pars = start.pars, pt.quad.id = pt.quad.id, row.id = order(c(pres.rows, quad.rows)), fixed.names = fixed.names, bias.names = bias.names, random.names = random.names, bf.info = bf.info, bias.bf.info = bias.bf.info, basis.functions = basis.functions, args = arg.info)

  }
}
