#' Spatially correlated, approximate modelling of presences in R
#'
#' @description This is the main function for modelling presences within the \code{scampr} framework. The type of model will depend on the arguments provided. This can be used to model point patterns as a log-Gaussian Cox process (LGCP) or Inhomogeneous Poisson process (IPP), as well as, jointly fitting a model to presence-only and presence/absence data if both are provided. This function can also fit binary regression with a complimentary log-log link function (with optional spatial random effects) when \code{model.type = "PA"}.
#'
#' When \code{model.type = "PO"}, the function will fit either an IPP, or LGCP model to the point pattern (depending on argument \code{include.sre}). This uses numerical quadrature (provided with the data, see e.g. scampr::gorillas) to approximate the spatial integral. If fitting a LGCP, uses one of either a Laplace or variational approximation to marginalise over the latent field according to \code{sre.approx}.
#'
#' When \code{model.type = "PA"}, the function will fits a binary regression model to presence/absence data using a complimentary log-log link function. Can accommodate an approximate latent field as spatial random effects (depending on argument \code{include.sre}).
#'
#' When \code{model.type = "IDM"}, the function jointly fits a model to presence-only and presence/absence data as linked by response to environmental predictors provided in each formula. The presence-only formula must also contain biasing predictors to account for opportunistic collection. If argument \code{include.sre} equals \code{TRUE}, then the two data sources will additionally share a latent Gaussian random field. Observer bias in the presence-only data can be modelled using measured variables via the argument \code{bias.formula} or as unmeasured effects with a second latent Gaussian field if \code{latent.po.biasing == TRUE} (the default). Either formula accepts an offset term - an offset in \code{formula} will apply to the presence/absence data while one in \code{bias.formula} will apply to the presence-only data.
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the fixed effects of the model to be fitted. The 'response' must be a binary that indicates whether a datum is a presence or: quadrature point (for point process models)/ absence (for binary models). See GLM function for further formula details.
#' @param data a data frame containing response and predictors within \code{formula}. Usually defaults to the point pattern and quadrature but will attempt to cleverly determine this via \code{model.type}, e.g., this is assumed to be presence/absence data when \code{model.type='PA'}. Otherwise, use \code{po.data} and/or \code{pa.dtaa} to supply the appropriate data frame explicitly.
#' @param bias.formula an object of class "formula" (or one that can be coerced to that class). This is a symbolic description of the predictors included to account for bias in the presence-only data (no response term is needed).
#' @param pa.data an optional data frame. When fitting an integrated data model use this to pass in the presence/absence data.
#' @param po.data an optional data frame. Can be use to explicitly define the presence-only data (and quadrature) when fitting an integrated data model.
#' @param coord.names a vector of character strings describing the column names of the coordinates in both data frames.
#' @param quad.weights.name a character string of the column name of quadrature weights in the data.
#' @param include.sre a logical indicating whether to fit the model with spatial random effects (SRE).
#' @param sre.approx a character string indicating the type of approximation to be used to marginalise over the spatial random effects. May be one of 'laplace' or 'variational'.
#' @param model.type a character string indicating the type of data to be used. May be one of 'PO' (for a presence-only PPM) or 'PA' (for a presence/absence Binary GLM) or 'IDM' (for an integrated data model).
#' @param basis.functions an optional object of class 'Basis' created by \code{FRK::auto_basis()} or 'bf.df' created by \code{scampr::simple_basis()}. Either object describes a set of basis functions for approximating the latent Gaussian field. If NULL the model will use default \code{FRK::auto_basis()} with \code{max_basis = 0.25 * # of points}.
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to use sparse or dense matrix computations for the basis functions created.
#' @param starting.pars an optional named list or scampr model object that gives warm starting values for the parameters of the model.
#' @param latent.po.biasing a logical, applying only to IDM, indicating whether biasing in the presence-only data should be accounted for via an additional latent field. Default is true as this is the most flexible approach unless good measured PO biasing variables are available.
#' @param po.biasing.basis.functions an optional extra set of basis functions that can be used when \code{latent.po.biasing = TRUE}, otherwise \code{basis.functions} are used.
#' @param prune.bfs an integer indicating the number of presence-only records required within a basis function's radius for it NOT to be pruned. Applies to the PO and IDM model (additionally, within the presence-only biasing basis functions in the IDM case) to assist with stability in model convergence. Default is zero, i.e. no pruning.
#' @param se a logical indicating whether standard errors should be calculated.
#' @param starting.pars an optional named list or previously fit scampr model object that gives warm starting values for the parameters of the model.
#' @param subset an optional vector describing a subset of the data to be used. Not applicable to integrated data models.
#' @param maxit a numeric indicating the maximum number of iterations for the optimizer. Default is 100 as the optimizer uses a gradient based approach.
#'
#' @return a scampr model object
#' @export
#'
#' @importFrom methods as
#' @importFrom stats optim terms
#' @importFrom TMB sdreport
#'
#' @examples
#' # Get the flora data for one of the species
#' dat_po <- flora$po$sp1
#' dat_pa <- flora$pa
#'
#' # obtain a sample of 10,000 quadrature points for the point process model
#' set.seed(1)
#' quad.pts <- flora$quad[sample(1:nrow(flora$quad), 10000, replace = F), ]
#' set.seed(NULL)
#'
#' # Attach the quadrature points to the presence-only data
#' dat_po <- rbind.data.frame(dat_po, quad.pts)
#'
#' # Ensure the "response" variable in each data set shares the same name
#' dat_po$presence <- dat_po$pres
#' dat_pa$presence <- dat_pa$sp1
#'
#' # Fit models without a latent effects
#' # Point Process Model
#' m.ipp <- scampr(presence ~ MNT + D.Main, dat_po, include.sre = F)
#' # Binary Regression
#' m.bin <- scampr(presence ~ MNT, dat_pa, include.sre = F, model.type = "PA")
#' # Integrated Data Model
#' m.comb <- scampr(presence ~ MNT, dat_po, bias.formula = ~ D.Main,
#' pa.data = dat_pa, include.sre = F, model.type = "IDM", latent.po.biasing = F)
#'
#' # Set up a simple 2D grid of basis functions to fit a LGCP model to the data
#' bfs <- simple_basis(nodes.on.long.edge = 9, data = dat_po)
#'
#' \dontrun{
#' # Fit with a shared latent field (LGCP) #
#' # Point Process Model
#' m.lgcp <- scampr(presence ~ MNT + D.Main, dat_po, basis.functions = bfs)
#' # Binary Regression with spatial random effects
#' m.bin_w_sre <- scampr(presence ~ MNT, dat_pa, basis.functions = bfs, model.type = "PA")
#' # Integrated Data Model with spatial random effects
#' m.comb_w_sre <- scampr(presence ~ MNT, dat_po, ~ D.Main,
#' dat_pa, basis.functions = bfs, model.type = "IDM")
#' }
scampr <- function(formula, data, bias.formula, pa.data, po.data, coord.names = c("x", "y"), quad.weights.name = "quad.size", include.sre = TRUE, sre.approx = c("variational", "laplace"), model.type = c("PO", "PA", "IDM"), basis.functions, bf.matrix.type = c("sparse", "dense"), latent.po.biasing = TRUE, po.biasing.basis.functions, prune.bfs = 4, se = TRUE, starting.pars, subset, maxit = 100, ...) {

  ## checks ####################################################################

  # checks for parameters of restricted strings
  model.type <- match.arg(model.type)
  sre.approx <- match.arg(sre.approx)
  bf.matrix.type <- match.arg(bf.matrix.type)

  # determine the data frames supplied (TODO: this is probably a mess of logic gates but I have done this in a rush)
  if (missing(data) & missing(po.data) & missing(pa.data)) {
    if (model.type == "IDM") {
      stop("At least two of 'data', 'pa.data' and 'po.data' must be supplied for an integrated data model")
    } else {
      stop("At least one of 'data', 'pa.data' and 'po.data' must be supplied")
    }
  }
  # individual model cases:
  if (model.type == "PO") {
    if (missing(data) & missing(po.data) & !missing(pa.data)) {
      warning("'pa.data' supplied will be assumed to be presence-only data for model.type = 'PO'...")
      data <- pa.data
    } else if (!missing(data) & !missing(po.data)) {
      warning("Both 'data' and 'po.data' supplied for model.type = 'PO'\nModel will be fitted to the 'po.data'...")
      data <- po.data
    } else if (missing(data) & !missing(po.data)) {
      data <- po.data
    }
  }
  if (model.type == "PA") {
    if (missing(data) & !missing(po.data) & missing(pa.data)) {
      warning("'po.data' supplied will be assumed to be presence/absence data for model.type = 'PA'...")
      data <- po.data
    } else if (!missing(data) & !missing(pa.data)) {
      warning("Both 'data' and 'pa.data' supplied for model.type = 'PA'\nModel will be fitted to the 'pa.data'...")
      data <- pa.data
    } else if (missing(data) & !missing(pa.data)) {
      data <- pa.data
    }
  }
  if (model.type == "IDM") {
    if (!missing(data) & !missing(po.data) & !missing(pa.data)) {
      warning("Each of 'data', 'po.data' and 'pa.data' supplied for model.type = 'IDM'.\n'data' will be ignored...")
      data <- po.data
    } else if (!missing(data) & !missing(po.data) & missing(pa.data)) {
      pa.data <- data
      data <- po.data
    } else if (missing(data) & !missing(po.data) & !missing(pa.data)) {
      data <- po.data
    }
  }


  # get the variables for checking
  resp <- all.vars(formula[[2]])
  # pred <- labels(stats::terms(formula)) # THIS BREAKS FOR I(X^2) FOR EXAMPLE
  pred <- all.vars(formula[[3]])

  if (length(resp) != 1) {
    stop("'formula' can only take a single response")
  }
  if (!all(c(resp, pred) %in% colnames(data))) {
    stop("'data' does not contain all of the terms in 'formula'")
  }
  if ((!all(coord.names %in% colnames(data))) & (model.type != "PA")) {
    stop(paste0("'coord.names', ", coord.names, ", not found 'data'"))
  }
  if ((!quad.weights.name %in% colnames(data)) & (model.type != "PA")) {
    stop(paste0("'quad.weights.name', ", quad.weights.name, ", not found in 'data'"))
  }
  if (!missing(bias.formula)) {
    if (is(bias.formula, "formula")) {
      if (!all(labels(terms(bias.formula)) %in% colnames(data))) {
        stop("'data' does not contain the bias.formula terms")
      }
    } else {
      stop("'bias.formula' must be of class formula")
    }
  }
  if (!is.logical(se)) {
    stop(paste0("'se' must be a logcial indicating whether or not to calculate standard errors"))
  }
  if (!is.logical(include.sre)) {
    stop(paste0("'include.sre' must be a logcial indicating whether or not to fit a spatial random effects model (SRE)"))
  }
  if (model.type %in% c("PA", "IDM") & sre.approx == "variational" & include.sre) {
    # warning(paste0("'model.type' = ", model.type, " is not compatible with a variational approx. Will instead use a Laplace approx."))
    sre.approx <- "laplace"
  }
  # if (model.type == "PA" & latent.po.biasing) {
  #   warning(paste0("'model.type' = ", model.type, " is not compatible with 'latent.po.biasing' - will be ignored"))
  #   latent.po.biasing <- FALSE
  # }
  # if (model.type == "PO" & latent.po.biasing) {
  #   warning(paste0("'latent.po.biasing' is not identifiable for 'model.type' = ", model.type, " using spatial random effects - will be ignored"))
  #   latent.po.biasing <- FALSE
  # }
  # if (model.type == "PO" & latent.po.biasing) {
  #   if (include.sre) {
  #     warning(paste0("'latent.po.biasing' is not identifiable for 'model.type' = ", model.type, " using spatial random effects - will be ignored"))
  #   } else {
  #     stop(paste0("Please use 'include.sre' = TRUE instead of 'latent.po.biasing' for an Inhomogeneous Poisson Process model"))
  #   }
  #   latent.po.biasing <- FALSE
  # }
  if (model.type == "IDM") {
    if (missing(pa.data)) {
      stop("Please provide presence/absence data as 'pa.data' for the Integrated Data Model")
    } else {
      if (!all(c(resp, pred) %in% colnames(pa.data))) {
        stop("'pa.data' does not contain all of the terms in 'formula'")
      }
    }
  }
  if (!latent.po.biasing & !missing(po.biasing.basis.functions)) {
    warning("'po.biasing.basis.functions' provided without indicating a need to account for presence-only biasing using spatial random effects - will be ignored")
  }

  # create 3 level indicator for model type and marginal approximation type
  if (include.sre) {
    approx.type <- sre.approx
  } else {
    approx.type <- "not_sre"
  }

  ##############################################################################

  ## NA action #################################################################

  # default na.action is to remove any data rows with na (for terms involved in the model)
  if (model.type != "PA") {
    rm.rows <- attr(na.omit(data[ , c(coord.names, quad.weights.name, resp, pred)]), "na.action")
  } else {
    rm.rows <- attr(na.omit(data[ , c(resp, pred)]), "na.action")
  }
  if (!is.null(rm.rows)) {
    data <- data[-rm.rows, ]
  }

  if (!missing(pa.data)) {
    rm.rows.df2 <- attr(na.omit(pa.data[ , c(resp, pred)]), "na.action")
    if (!is.null(rm.rows.df2)) {
      pa.data <- pa.data[-rm.rows.df2, ]
    }
  }
  ##############################################################################

  ## Sub-setting (not available for IDM) #######################################

  # apply sub-setting to data if supplied (with checks)
  if (model.type != "IDM") {
    if (!missing(subset)) {
      if (!is.vector(subset)) {
        stop("subset must be a vector")
      } else {
        if (is(subset, "logical")) {
          if (length(subset) != nrow(data)) {
            stop("Logical subset must be of same dimension as data provided")
          } else {
            data <- data[subset, ]
          }
        } else if (is(subset, "integer") | is(subset, "numeric")) {
          if (!all(subset %in% 1:nrow(data))) {
            stop("numeric subset must include row numbers of the data provided")
          }
          data <- data[subset, ]
        } else {
          stop("subset has an incorrect format")
        }
      }
    }
  } else {
    if (!missing(subset)) {
      warning("'subset' is not currently applicable for Integrated Data Models and will be ignored")
    }
  }
  ##############################################################################

  ## Create the TMB data and parameter inputs ##################################

  mc <- match.call() # gets the arguments (must be updated to the alterations above)

  call.list <- as.list(mc)
  # alter the call according to alterations made above
  call.list$bf.matrix.type <- bf.matrix.type
  call.list$model.type <- model.type
  call.list$sre.approx <- sre.approx
  call.list$data <- data
  if (!missing(pa.data)) {
    call.list$pa.data <- pa.data
  }
  mod.call <- as.call(call.list)
  # remove unused terms from the call
  call.list <- call.list[!names(call.list) %in% c("subset", "se", "include.sre", "sre.approx", "maxit")]
  # add approx.type to the call for get.TMB.data.input()
  call.list$approx.type <- approx.type

  # remove the call function
  call.list[[1]] <- NULL

  # quick fix for formulae - re-assign here so they are found in the TMB input call
  call.list$formula <- as.formula(formula)
  if (!missing(bias.formula)) {
    call.list$formula <- as.formula(bias.formula)
  }

  # get the TMB inputs
  inputs <- do.call("get.TMB.data.input", call.list)

  ##############################################################################

  ## Create the TMB Objective Function and Optimise ############################

  obj <- make.objective.function(inputs, maxit)
  try(assign("fit.time", system.time(assign("res", do.call("optim", obj)))), silent = T)
  if (!exists("res")) {
    stop(paste0("Fit failed: initial objective function values: ", obj$fn(), ". Try warm starting parameters or a simpler model"))
  }

  res$cpu <- c(opt = as.numeric(fit.time[3]), se = NA)

  ##############################################################################

  ## Handle Standard Errors  ###################################################

  # set up a flag for failure of standard error calculation
  res$se.flag <- 0

  # get standard errors if required
  if (se) {
    try(assign("se.time", system.time(assign("tmp.estimates", summary(TMB::sdreport(obj))))), silent = T)
    if (!exists("tmp.estimates")) {
      stop("Calculating Std. Errors failed. Try 'se' = FALSE")
    }
    # add the timing
    res$cpu["se"] <- as.numeric(se.time[3])
    # check if any standard errors where infinite, NaN, or NA and flag the model
    if (any(is.na(tmp.estimates[,2]) | any(is.infinite(tmp.estimates[,2])) | any(is.nan(tmp.estimates[,2])))) {
      res$se.flag <- 1
    }
  } else { # Extract the point estimates only but place in the same format as for se = TRUE

    # add in the posterior means of the random coefficients and their corresponding posterior/prior variance components
    if (inputs$args$approx.type == "laplace") {
      # posterior means of the random coefficients
      re.posterior.means <- obj$env$parList()$random
      names(re.posterior.means) <- rep("random", length(re.posterior.means))
      # posterior variances of the random coefficients
      var.comp_post <- NULL # don't exist in the Laplace approximation
      # prior variances of the random coefficients
      var.comp_prior <- exp(obj$env$parList()$log_variance_component)^2
      names(var.comp_prior) <- rep("PriorVar", length(var.comp_prior))

    } else if (inputs$args$approx.type == "variational") {
      # posterior means of the random coefficients
      re.posterior.means <- NULL # already included in res$par in this case
      # posterior variances of the random coefficients
      var.comp_post <- exp(obj$env$parList()$log_variance_component)
      names(var.comp_post) <- rep("PosteriorVar", length(var.comp_post))
      # prior variances of the random coefficients
      var.comp_prior <- sapply(split(var.comp_post + (obj$env$parList()$random^2), inputs$bf.info$res), mean) # calculate the mean of the sum of variances and squared posterior means, within each basis function resolution
      names(var.comp_prior) <- rep("PriorVar", length(var.comp_prior))
    } else {
      re.posterior.means <- NULL
      var.comp_post <- NULL
      var.comp_prior <- NULL
    }
    # TODO: add in the random_bias components
    if (inputs$args$random.bias.type %in% c("field1", "field2")) {
      # posterior means of the random bias coefficients
      re.posterior.means_bias <- obj$env$parList()$random_bias
      names(re.posterior.means_bias) <- rep("random_bias", length(re.posterior.means_bias))
      # prior variances of the random bias coefficients (included with the other prior variances)
      tmp.bias_var <- exp(obj$env$parList()$log_variance_component_bias)^2
      names(tmp.bias_var) <- rep("PriorVar_bias", length(tmp.bias_var))
      var.comp_prior <- c(tmp.bias_var, var.comp_prior)
    } else {
      re.posterior.means_bias <- NULL
    }

    # combine into the estimates data frame equivalent produced when se = T
    tmp.estimates <- cbind(Estimate = c(res$par, re.posterior.means, var.comp_post, re.posterior.means_bias, var.comp_prior),
                           `Std. Error` = rep(NA, sum(length(res$par) + length(re.posterior.means) + length(var.comp_post) + length(re.posterior.means_bias) + length(var.comp_prior))))
  }

  ##############################################################################

  ## Additional Information for the model object  ##############################

  # add coefficient information to the results list (useful for S3 methods coeff, etc.)
  res$coefficients <- res$par
  # fixed effect names
  names(res$coefficients)[names(res$par) == "fixed"] <- inputs$fixed.names
  # biasing variable effect names
  if (any(names(res$par) == "bias")) {
    names(res$coefficients)[names(res$par) == "bias"] <- inputs$bias.names
  }
  # random effect coefficient means (in the case of a variational model)
  if (any(names(res$par) == "random")) {
    names(res$coefficients)[names(res$par) == "random"] <- switch(inputs$args$approx.type,
                                                                  not_sre = stop("random effects found in model without spatial random effects!"),
                                                                  variational = paste0("VA Posterior Mean (", inputs$random.names, ")"),
                                                                  laplace = stop("random effect means found in a Laplace-based model")
    )
  }
  # random effect variance components
  if (any(names(res$par) == "log_variance_component")) {
    names(res$coefficients)[names(res$par) == "log_variance_component"] <- switch(inputs$args$approx.type,
                                                                                  not_sre = stop("random effects found in model without spatial random effects!"),
                                                                                  variational = paste0("VA Posterior log variance (", inputs$random.names, ")"),
                                                                                  laplace = paste0("Prior log sd(u) (res. ", 1:length(inputs$tmb.data$bf_per_res), ")")
    )
  }
  # biasing random effect variance components
  if (any(names(res$par) == "log_variance_component_bias")) {
    if (is.null(inputs$tmb.data$bias_bf_per_res)) {
      K <- length(inputs$tmb.data$bf_per_res)
    } else {
      K <- length(inputs$tmb.data$bias_bf_per_res)
    }
    names(res$coefficients)[names(res$par) == "log_variance_component_bias"] <- switch(inputs$args$approx.type,
                                                                                not_sre = stop("random effects found in model without spatial random effects!"),
                                                                                variational = stop("biasing random effects found in a VA-based model"),
                                                                                laplace = paste0("Prior log sd(tau) (res. ", 1:K, ")")
    )
  }

  # add the fixed effects data frame
  if (sum(rownames(tmp.estimates) == "fixed") == 1) { # check for a single fixed effect to adjust the resulting data frame
    res$fixed.effects <- data.frame(t(tmp.estimates[rownames(tmp.estimates) == "fixed", ]))
    colnames(res$fixed.effects)[2] <- "Std. Error" # need to correct the white space in name
  } else {
    res$fixed.effects <- tmp.estimates[rownames(tmp.estimates) == "fixed", ]
  }
  rownames(res$fixed.effects) <- inputs$fixed.names

  # add the biasing effects if present
  if (any(rownames(tmp.estimates) == "bias")) {
    if (sum(rownames(tmp.estimates) == "bias") == 1) { # check for a single fixed effect to adjust the resulting data frame
      res$fixed.bias.effects <- data.frame(t(tmp.estimates[rownames(tmp.estimates) == "bias", ]))
      colnames(res$fixed.bias.effects)[2] <- "Std. Error" # need to correct the white space in columnn name
      rownames(res$fixed.bias.effects) <- inputs$bias.names
    } else {
      res$fixed.bias.effects <- tmp.estimates[rownames(tmp.estimates) == "bias", ]
      rownames(res$fixed.bias.effects) <- inputs$bias.names
    }
  } else {
    res$fixed.bias.effects <- NULL
  }

  # add the random effects if present
  if (any(row.names(tmp.estimates) == "random")) {
    if (sum(rownames(tmp.estimates) == "random") == 1) { # check for a single random effect to adjust the resulting data frame
      res$random.effects <- cbind(data.frame(t(tmp.estimates[rownames(tmp.estimates) == "random", ])), inputs$bf.info)
      colnames(res$random.effects)[2] <- "Std. Error" # need to correct the white space in name
    } else {
      res$random.effects <- cbind(tmp.estimates[rownames(tmp.estimates) == "random", ], inputs$bf.info)
      rownames(res$random.effects) <- inputs$random.names
    }
  } else {
    res$random.effects <- NULL
  }

  # add the random biasing effects if present
  if (any(row.names(tmp.estimates) == "random_bias")) {
    if (sum(rownames(tmp.estimates) == "random_bias") == 1) { # check for a single random bias effect to adjust the resulting data frame
      res$random.bias.effects <- data.frame(t(tmp.estimates[rownames(tmp.estimates) == "random_bias", ]))
      colnames(res$random.bias.effects)[2] <- "Std. Error" # need to correct the white space in name
    } else {
      res$random.bias.effects <- tmp.estimates[rownames(tmp.estimates) == "random_bias", ]
      # check and attach the appropriate basis function information
      if (!is.null(inputs$bias.bf.info)) {
        res$random.bias.effects <- cbind(res$random.bias.effects, inputs$bias.bf.info)
      } else {
        res$random.bias.effects <- cbind(res$random.bias.effects, inputs$bf.info)
      }
      rownames(res$random.bias.effects) <- inputs$random.bias.names
    }
  } else {
    res$random.bias.effects <- NULL
  }

  # add the posterior variance/sd estimates if present
  if (any(grepl("PosteriorVar", row.names(tmp.estimates), fixed = T))) {
    if (sum(grepl("PosteriorVar", row.names(tmp.estimates), fixed = T)) == 1) { # check for a single posterior variance to adjust the resulting data frame
      res$post.variances <- data.frame(t(tmp.estimates[grepl("PosteriorVar", row.names(tmp.estimates), fixed = T), ]))
      colnames(res$post.variances)[2] <- "Std. Error" # need to correct the white space in name
      rownames(res$post.variances) <- row.names(tmp.estimates)[grepl("PosteriorVar", row.names(tmp.estimates), fixed = T)]
    } else {
      res$post.variances <- tmp.estimates[grepl("PosteriorVar", row.names(tmp.estimates), fixed = T), ]
    }
  } else {
    res$post.variances <- NULL
  }

  # add the prior variance/sd estimates if present
  if (any(grepl("PriorVar", row.names(tmp.estimates), fixed = T))) {
    if (sum(grepl("PriorVar", row.names(tmp.estimates), fixed = T)) == 1) { # check for a single fixed effect to adjust the resulting data frame
      res$prior.variances <- data.frame(t(tmp.estimates[grepl("PriorVar", row.names(tmp.estimates), fixed = T), ]))
      colnames(res$prior.variances)[2] <- "Std. Error" # need to correct the white space in name
      rownames(res$prior.variances) <- row.names(tmp.estimates)[grepl("PriorVar", row.names(tmp.estimates), fixed = T)]
    } else {
      res$prior.variances <- tmp.estimates[grepl("PriorVar", row.names(tmp.estimates), fixed = T), ]
    }
  } else {
    res$prior.variances <- NULL
  }
  ##############################################################################

  # add in the fitted values ###################################################

  # convoluted if statements to get all of the potential components of the linear predictor
  if (inputs$args$model.type == "PA") {
    # calculate the fixed effect component
    Xbeta_pa <- as.vector(inputs$tmb.data$X_PA %*% res$fixed.effects[, 1L])

    # get the random effect component if the model has spatial random effects
    if (inputs$args$approx.type != "not_sre") {
      Zmu_pa <- as.vector(inputs$tmb.data$Z_PA %*% res$random.effects[, 1L])
    } else {
      Zmu_pa <- NULL
    }

    # zero off the remaining components
    Xbeta <- NULL
    Btau <- NULL
    Zmu <- NULL
    Z2mu2 <- NULL
    # res$fitted.values <- as.vector(inputs$tmb.data$X_PA %*% res$fixed.effects[, 1L]) + as.vector(inputs$tmb.data$Z_PA %*% res$random.effects[, 1L])

    } else { # in either a PO or IDM case

      # need to put the presence-only data back together (pres and quad) and return to the ordering of the original data provided
      if (ncol(inputs$tmb.data$X_PO_pres) == 1) { # when there is only one variable
        oldX <- matrix(c(inputs$tmb.data$X_PO_pres, inputs$tmb.data$X_PO_quad)[inputs$row.id], ncol = 1)
      } else {
        oldX <- rbind(inputs$tmb.data$X_PO_pres, inputs$tmb.data$X_PO_quad)[inputs$row.id, ]
      }
      # calculate the fixed effect component (shared by all)
      Xbeta <- as.vector(oldX %*% res$fixed.effects[, 1L])

      # get the random effect component if the model has spatial random effects
      if (inputs$args$approx.type != "not_sre") {
        # if (ncol(inputs$tmb.data$Z_PO_pres) == 1) { # when there is only one variable
        #   oldZ <- matrix(c(inputs$tmb.data$Z_PO_pres, inputs$tmb.data$Z_PO_quad)[inputs$row.id], ncol = 1)
        # } else {
        #   oldZ <- rbind(inputs$tmb.data$Z_PO_pres, inputs$tmb.data$Z_PO_quad)[inputs$row.id, ]
        # } # THIS SHOULD BE REDUNDANT SINCE WE ALREADY MAKE THE INPUT A SINGLE COLUMN MATRIX!
        oldZ <- rbind(inputs$tmb.data$Z_PO_pres, inputs$tmb.data$Z_PO_quad)[inputs$row.id, ]
        if (ncol(inputs$tmb.data$Z_PO_pres) == 1) { # when there is only one variable
          Zmu <- as.vector(oldZ * res$random.effects[, 1L])
        } else {
          Zmu <- as.vector(oldZ %*% res$random.effects[, 1L])
        }
      } else {
        Zmu <- NULL
      }

      # get the biasing fixed effects component if included in the model
      if (inputs$args$fixed.bias.type != "missing") {
        if (ncol(inputs$tmb.data$B_PO_pres) == 1) {
          oldB <- matrix(c(inputs$tmb.data$B_PO_pres, inputs$tmb.data$B_PO_quad)[inputs$row.id], ncol = 1)
        } else {
          oldB <- rbind(inputs$tmb.data$B_PO_pres, inputs$tmb.data$B_PO_quad)[inputs$row.id, ]
        }
        Btau <- as.vector(oldB %*% res$fixed.bias.effects[, 1L])
      } else {
        Btau <- NULL
      }

      # get the biasing random effects component if included in the model
      if (inputs$args$random.bias.type != "none") {
        if (inputs$args$random.bias.type == "field1") {
          # if (ncol(inputs$tmb.data$Z_PO_pres) == 1) { # when there is only one variable
          #   oldZ <- matrix(c(inputs$tmb.data$Z_PO_pres, inputs$tmb.data$Z_PO_quad)[inputs$row.id], ncol = 1)
          # } else {
          #   oldZ <- rbind(inputs$tmb.data$Z_PO_pres, inputs$tmb.data$Z_PO_quad)[inputs$row.id, ]
          # } # THIS SHOULD BE REDUNDANT SINCE WE ALREADY MAKE THE INPUT A SINGLE COLUMN MATRIX!
          oldZ <- rbind(inputs$tmb.data$Z_PO_pres, inputs$tmb.data$Z_PO_quad)[inputs$row.id, ]
          if (ncol(inputs$tmb.data$Z_PO_pres) == 1) { # when there is only one variable
            Z2mu2 <- as.vector(oldZ * res$random.bias.effects[, 1L])
          } else {
            Z2mu2 <- as.vector(oldZ %*% res$random.bias.effects[, 1L])
          }
        } else {
          # if (ncol(inputs$tmb.data$Z2_PO_pres) == 1) { # when there is only one variable
          #   oldZ2 <- matrix(c(inputs$tmb.data$Z2_PO_pres, inputs$tmb.data$Z2_PO_quad)[inputs$row.id], ncol = 1)
          # } else {
          #   oldZ2 <- rbind(inputs$tmb.data$Z2_PO_pres, inputs$tmb.data$Z2_PO_quad)[inputs$row.id, ]
          # } # THIS SHOULD BE REDUNDANT SINCE WE ALREADY MAKE THE INPUT A SINGLE COLUMN MATRIX!
          oldZ2 <- rbind(inputs$tmb.data$Z2_PO_pres, inputs$tmb.data$Z2_PO_quad)[inputs$row.id, ]
          if (ncol(inputs$tmb.data$Z2_PO_pres) == 1) { # when there is only one variable
            Z2mu2 <- as.vector(oldZ2 * res$random.bias.effects[, 1L])
          } else {
            Z2mu2 <- as.vector(oldZ2 %*% res$random.bias.effects[, 1L])
          }
        }
      } else {
        Z2mu2 <- NULL
      }

      # get the PA components for the IDM
      if (inputs$args$model.type == "IDM") {
        # calculate the fixed effect component
        Xbeta_pa <- as.vector(inputs$tmb.data$X_PA %*% res$fixed.effects[, 1L])

        # get the random effect component if the model has spatial random effects
        if (inputs$args$approx.type != "not_sre") {
          Zmu_pa <- as.vector(inputs$tmb.data$Z_PA %*% res$random.effects[, 1L])
        } else {
          Zmu_pa <- NULL
        }
      } else {
        Xbeta_pa <- NULL
        Zmu_pa <- NULL
      }
    }

  # calculate the fitted values
  if (inputs$args$model.type == "PA") {
    res$fitted.values <- rowSums(cbind(Xbeta_pa, Zmu_pa))
    attr(res$fitted.values, "Xbeta") <- Xbeta_pa
    attr(res$fitted.values, "Zmu") <- Zmu_pa
  } else if (inputs$args$model.type == "PO") {
    res$fitted.values <- rowSums(cbind(Xbeta, Btau, Zmu, Z2mu2))
    attr(res$fitted.values, "Xbeta") <- Xbeta
    attr(res$fitted.values, "Btau") <- Btau
    attr(res$fitted.values, "Zmu") <- Zmu
    attr(res$fitted.values, "Z2mu2") <- Z2mu2
  } else {
    res$fitted.values <- rowSums(cbind(Xbeta, Btau, Zmu, Z2mu2))
    attr(res$fitted.values, "PA") <- rowSums(cbind(Xbeta_pa, Zmu_pa))
    attr(res$fitted.values, "Xbeta") <- Xbeta
    attr(res$fitted.values, "Btau") <- Btau
    attr(res$fitted.values, "Zmu") <- Zmu
    attr(res$fitted.values, "Z2mu2") <- Z2mu2
  }

  ##############################################################################

  # final additions for interfacing with other functions
  res$basis.per.res <- inputs$tmb.data$bf_per_res
  attr(res$basis.per.res, "bias") <- inputs$tmb.data$bias_bf_per_res
  res$basis.functions <- inputs$basis.functions
  if (inputs$args$random.bias.type == "field1") {
    res$po.biasing.basis.functions <- inputs$basis.functions
  } else {
    res$po.biasing.basis.functions <- inputs$po.biasing.basis.functions
  }
  res$basis.fn.info <- inputs$bf.info
  # adjust the TMB call list approx type
  call.list$approx.type <- inputs$args$approx.type
  res$approx.type <- approx.type # set this as the original approx.type
  res$starting.pars <- inputs$start.pars
  res$data <- data
  if (inputs$args$model.type == "IDM") {
    attr(res$data, "PA") <- pa.data
  }
  res$formula <- formula
  if (!missing(bias.formula)) {
    attr(res$formula, "bias") <- bias.formula
  }
  res$coord.names <- coord.names
  res$quad.weights.name <- quad.weights.name
  res$pt.quad.id <- inputs$pt.quad.id
  res$model.type <- inputs$args$model.type
  res$bf.matrix.type <- inputs$args$bf.matrix.type
  res$fixed.bias.type <- inputs$args$fixed.bias.type
  res$random.bias.type <- inputs$args$random.bias.type
  res$tmb.data <- inputs$tmb.data
  res$ll.components <- obj$report()
  res$call <- mod.call
  res$tmb.call.list <- call.list
  class(res) <- "scampr"

  return(res)

}
