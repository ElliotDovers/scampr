#' Spatially correlated, approximate modelling of presences in R
#'
#' @description This is the main function for modelling presences within the \code{scampr} framework. The type of model will depend on the arguments provided. This can be used to model point patterns as a log-Gaussian Cox process (LGCP) or Inhomogeneous Poisson process (IPP), as well as, jointly fitting a model to presence-only and presence/absence data if both are provided. This function can also fit binary regression with a complimentary log-log link function (with optional spatial random effects) when \code{data.type = "PA"}.
#'
#' When \code{data.type = "PO"}, the function will fit either an IPP, or LGCP model to the point pattern (depending on argument \code{model.type}). This uses numerical quadrature (provided with the data, see e.g. scampr::gorillas) to approximate the spatial integral. If fitting a LGCP, uses one of either a Laplace or variational approximation to marginalise over the latent field.
#'
#' When \code{data.type = "PA"}, the function will fits a binary regression model to presence/absence data using a complimentary log-log link function. Can accomodate an approximate latent field as spatial random effects (depending on argument \code{model.type}).
#'
#' When \code{data.type = "IDM}, the function jointly fits a model to presence-only and presence/absence data as linked by response to environmental predictors provided in each formula. The presence-only formula must also contain biasing predictors to account for opportunistic collection. If argument \code{model.type} is not equal to "ipp", then the two data sources will additionally share a latent Gaussian random field.
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the fixed effects of the model to be fitted. The 'response' must be a binary that indicates whether a datum is a presence or: quadrature point (for point process models)/ absence (for binary models). See GLM function for further formula details.
#' @param data a data frame containing response and predictors within \code{formula}.
#' @param bias.formula an object of class "formula" (or one that can be coerced to that class) OR the character string "latent". In the formula case, this is a symbolic description of the predictors included to account for bias in the presence-only data (no response term is needed). In the case of fitting an integrated data model, \code{bias.formula = "latent"} will fit an approximate latent Gaussian field to account for the bias.
#' @param IDM.presence.absence.df an optional data frame. When fitting an integrated data model use this to pass in the presence/absence data.
#' @param coord.names a vector of character strings describing the column names of the coordinates in both data frames.
#' @param quad.weights.name a character string of the column name of quadrature weights in the data.
#' @param model.type a character string indicating the type of model to be used. May be one of 'laplace' or 'variational' for Cox Processes involving spatially correlated errors or 'ipp' for a model that follows an inhomgeneous Poisson process.
#' @param data.type a character string indicating the type of data to be used. May be one of 'PO' (for a presence-only PPM) or 'PA' (for a presence/absence Binary GLM) or 'IDM' (for an integrated data model).
#' @param basis.functions an optional object of class 'Basis' created by \code{FRK::auto_basis()} or 'bf.df' created by \code{scampr::simple_basis()}. Either object describes a set of basis functions for approximating the latent Gaussian field. If NULL the model will use default \code{FRK::auto_basis()} with \code{max_basis = 0.25 * # of points}.
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to use sparse or dense matrix computations for the basis functions created.
#' @param starting.pars an optional named list or scampr model object that gives warm starting values for the parameters of the model.
#' @param po.biasing.basis.functions an optional extra set of basis functions that can be used when \code{bias.formula = "latent"}, otherwise \code{basis.functions} are used.
#' @param se a logical indicating whether standard errors should be calculated.
#' @param starting.pars an optional named list or previously fit scampr model object that gives warm starting values for the parameters of the model.
#' @param subset an optional vector describing a subset of the data to be used. Not applicable to integrated data models.
#' @param maxit a numeric indicating the maximum number of iterations for the optimizer. Default is 1000.
#'
#' @return a scampr model object
#' @export
#'
#' @importFrom methods as
#' @importFrom stats optim
#' @importFrom TMB sdreport
#'
#' @examples
#' # Get the flora data for one of the species
#' dat_po <- flora$po$sp1
#' dat_pa <- flora$pa
#'
#' # Attach the quadrature to the PO data
#' dat_po <- rbind.data.frame(dat_po, flora$quad)
#'
#' # Fit models without a latent effects (IPP) #
#' # Point Process Model
#' m.ipp <- scampr(pres ~ MNT + D.Main, dat_po, model.type = "ipp")
#' # Binary Regression
#' m.bin <- scampr(pa.formula = sp1 ~ MNT, pa.data = dat_pa, model.type = "ipp")
#' # Combined Data Model
#' m.comb <- scampr(pres ~ MNT + D.Main, dat_po, sp1 ~ MNT,
#' dat_pa, model.type = "ipp")
#'
#' # Set up a simple 2D grid of basis functions to fit a LGCP model to the data
#' bfs <- simple_basis(nodes.on.long.edge = 9, data = dat_po)
#'
#' \dontrun{
#' # Fit with a shared latent field (LGCP) #
#' # Point Process Model
#' m.lgcp <- scampr(pres ~ MNT + D.Main, dat_po, basis.functions = bfs)
#' # Binary Regression with spatial random effects
#' m.bin_w_sre <- scampr(pa.formula = sp1 ~ MNT, pa.data = dat_pa, basis.functions = bfs)
#' # Combined Data Model with spatial random effects
#' m.comb_w_sre <- scampr(pres ~ MNT + D.Main, dat_po, sp1 ~ MNT,
#' dat_pa, basis.functions = bfs)
#' }
scampr <- function(formula, data, bias.formula, IDM.presence.absence.df, coord.names = c("x", "y"), quad.weights.name = "quad.size", model.type = c("variational", "laplace", "ipp"), data.type = c("PO", "PA", "IDM"), basis.functions, bf.matrix.type = c("sparse", "dense"), po.biasing.basis.functions, se = TRUE, starting.pars, subset, maxit = 1000) {

  ## checks ####################################################################

  # checks for parameters of restricted strings
  data.type <- match.arg(data.type)
  model.type <- match.arg(model.type)
  bf.matrix.type <- match.arg(bf.matrix.type)

  # get the variables for checking
  resp <- all.vars(formula[[2]])
  pred <- labels(stats::terms(formula))

  if (length(resp) != 1) {
    stop("'formula' can only take a single response")
  }
  if (!all(c(resp, pred) %in% colnames(data))) {
    stop("'data' does not contain all of the terms in 'formula'")
  }
  if ((!all(coord.names %in% colnames(data))) & (data.type != "PA")) {
    stop(paste0("'coord.names', ", coord.names, ", not found 'data'"))
  }
  if ((!quad.weights.name %in% colnames(data)) & (data.type != "PA")) {
    stop(paste0("'quad.weights.name', ", quad.weights.name, ", not found in 'data'"))
  }
  if (!missing(bias.formula)) {
    if (bias.formula != "latent") {
      if (!all(labels(terms(bias.formula)) %in% colnames(data))) {
        stop("'data' does not contain the formula terms")
      }
    }
  }
  if (!is.logical(se)) {
    stop(paste0("'se' must be a logcial indicating whether or not to calculate standard errors"))
  }
  if (data.type %in% c("PA", "IDM") & model.type == "variational") {
    warning(paste0("'data.type' = ", data.type, " is not compatible with a variational approx. Will instead use a Laplace approx."))
    model.type <- "laplace"
  }
  if (data.type == "IDM") {
    if (missing(IDM.presence.absence.df)) {
      stop("Please provide presence/absence data as 'IDM.presence.absence.df' for the Integrated Data Model")
    } else {
      if (!all(c(resp, pred) %in% colnames(IDM.presence.absence.df))) {
        stop("'IDM.presence.absence.df' does not contain all of the terms in 'formula'")
      }
    }
  }
  ##############################################################################

  ## NA action #################################################################

  # default na.action is to remove any data rows with na (for terms involved in the model)
  if (data.type != "PA") {
    rm.rows <- attr(na.omit(data[ , c(coord.names, quad.weights.name, resp, pred)]), "na.action")
  } else {
    rm.rows <- attr(na.omit(data[ , c(resp, pred)]), "na.action")
  }
  if (!is.null(rm.rows)) {
    data <- data[-rm.rows, ]
  }

  if (!missing(IDM.presence.absence.df)) {
    rm.rows.df2 <- attr(na.omit(IDM.presence.absence.df[ , c(resp, pred)]), "na.action")
    if (!is.null(rm.rows.df2)) {
      IDM.presence.absence.df <- IDM.presence.absence.df[-rm.rows.df2, ]
    }
  }
  ##############################################################################

  ## Sub-setting (not available for IDM) #######################################

  # apply sub-setting to data if supplied (with checks)
  if (data.type != "IDM") {
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
  }
  ##############################################################################

  ## Create the TMB data and parameter inputs ##################################

  mc <- match.call() # gets the arguments (after being altered above)
  mc[[1]] <- quote(scampr::get.TMB.data.input) # tell the call the function to be evaluated
  inputs <- eval(mc, envir = parent.frame())
  ##############################################################################

  ## Create the TMB Objective Function and Optimise ############################

  obj <- make.objective.function(inputs, maxit)
  try(assign("fit.time", system.time(assign("res", do.call("stats::optim", obj)))), silent = T)
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
  } else {
    tmp.estimates <- cbind(Estimate = res$par, `Std. Error` = rep(NA, length(res$par)))
  }

  ##############################################################################

  ## Additional Information for the model object  ##############################

  # add coefficient information to the results list (useful for S3 methods coeff, etc.)
  res$coefficients <- res$par
  # fixed effect names
  names(res$coefficients)[names(res$par) == "fixed"] <- inputs$fixed.names
  # biasing variable effect names
  if (!is.null(inputs$bias.names)) {
    names(res$coefficients)[names(res$par) == "bias"] <- inputs$bias.names
  }
  # random effect coefficient means (in the case of a variational model)
  if (any(names(res$par) == "random")) {
    names(res$coefficients)[names(res$par) == "random"] <- switch(inputs$args$model.type,
                                                                  ipp = stop("random effects found in non-random model"),
                                                                  variational = paste0("VA Posterior Mean (", inputs$random.names, ")"),
                                                                  laplace = stop("random effect means found in a Laplace-based model")
    )
  }
  # random effect variance components
  if (any(names(res$par) == "log_variance_component")) {
    names(res$coefficients)[names(res$par) == "log_variance_component"] <- switch(inputs$args$model.type,
                                                                                  ipp = stop("random effects found in non-random model"),
                                                                                  variational = paste0("VA Posterior log variance (", inputs$random.names, ")"),
                                                                                  laplace = paste0("Prior log sd(u) (res. ", 1:length(inputs$tmb.data$bf_per_res), ")")
    )
  }
  # biasing random effect variance components
  if (any(names(res$par) == "log_variance_component_bias")) {
    names(res$coefficients)[names(res$par) == "log_variance_component_bias"] <- switch(inputs$args$model.type,
                                                                                ipp = stop("random effects found in non-random model"),
                                                                                variational = stop("biasing random effects found in a VA-based model"),
                                                                                laplace = paste0("Prior log sd(tau) (res. ", 1:length(inputs$tmb.data$bf_per_res), ")")
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
  if (!is.null(inputs$bias.names)) {
    if (sum(rownames(tmp.estimates) == "bias") == 1) { # check for a single fixed effect to adjust the resulting data frame
      res$bias.effects <- data.frame(t(tmp.estimates[rownames(tmp.estimates) == "bias", ]))
      colnames(res$bias.effects)[2] <- "Std. Error" # need to correct the white space in columnn name
      rownames(res$bias.effects) <- inputs$bias.names
    } else {
      res$bias.effects <- tmp.estimates[rownames(tmp.estimates) == "bias", ]
      if (inputs$args$bias.type %in% c("latent", "new_latent")) {
        if (!is.null(inputs$bias.bf.info)) {
          res$bias.effects <- cbind(res$bias.effects, inputs$bias.bf.info)
        } else {
          res$bias.effects <- cbind(res$bias.effects, inputs$bf.info)
        }
      }
      rownames(res$bias.effects) <- inputs$bias.names
    }
  } else {
    res$bias.effects <- NULL
  }

  # add the random effects if present
  if (any(row.names(tmp.estimates) == "random")) {
    if (sum(rownames(tmp.estimates) == "random") == 1) { # check for a single fixed effect to adjust the resulting data frame
      res$random.effects <- data.frame(t(tmp.estimates[rownames(tmp.estimates) == "random", ]))
      colnames(res$random.effects)[2] <- "Std. Error" # need to correct the white space in name
    } else {
      res$random.effects <- cbind(tmp.estimates[rownames(tmp.estimates) == "random", ], inputs$bf.info)
      rownames(res$random.effects) <- inputs$random.names
    }
  } else {
    res$random.effects <- NULL
  }

  # add the variance/sd estimates if present
  if (any(grepl("PriorVar", row.names(tmp.estimates), fixed = T))) {
    if (sum(grepl("PriorVar", row.names(tmp.estimates), fixed = T)) == 1) { # check for a single fixed effect to adjust the resulting data frame
      res$variances <- data.frame(t(tmp.estimates[grepl("PriorVar", row.names(tmp.estimates), fixed = T), ]))
      colnames(res$variances)[2] <- "Std. Error" # need to correct the white space in name
      rownames(res$variances) <- row.names(tmp.estimates)[grepl("PriorVar", row.names(tmp.estimates), fixed = T)]
    } else {
      res$variances <- tmp.estimates[grepl("PriorVar", row.names(tmp.estimates), fixed = T), ]
    }
  } else {
    res$variances <- NULL
  }
  ##############################################################################

  # add in the fitted values ###################################################
  if (inputs$args$model.type != "ipp") {
    if (inputs$args$data.type == "PA") {
      res$fitted.values <- as.vector(inputs$tmb.data$X_PA %*% res$fixed.effects[, 1L]) + as.vector(inputs$tmb.data$Z_PA %*% res$random.effects[, 1L])

    } else { # in either a PO or IDM case

      # need to put the presence-only data back together (pres and quad) and return to the ordering of the original data provided
      oldX <- rbind(inputs$tmb.data$X_PO_pres, inputs$tmb.data$X_PO_quad)[inputs$po.info[ , "row.id"], ]
      oldZ <- rbind(inputs$tmb.data$Z_PO_pres, inputs$tmb.data$Z_PO_quad)[inputs$po.info[ , "row.id"], ]

      if (inputs$args$bias.type == "none") { # when no biasing covariates are present
        res$fitted.values <- as.vector(oldX %*% res$fixed.effects[, 1L]) + as.vector(oldZ %*% res$random.effects[, 1L])
      } else if (inputs$args$bias.type == "covariates") {

        oldB <- rbind(inputs$tmb.data$B_PO_pres, inputs$tmb.data$B_PO_quad)[inputs$po.info[ , "row.id"], ]
        res$fitted.values <- as.vector(oldX %*% res$fixed.effects[, 1L]) + as.vector(oldB %*% res$bias.effects[, 1L]) + as.vector(oldZ %*% res$random.effects[, 1L])

      } else if (inputs$args$bias.type == "latent") {

        res$fitted.values <- as.vector(oldX %*% res$fixed.effects[, 1L]) + as.vector(oldZ %*% res$bias.effects[, 1L]) + as.vector(oldZ %*% res$random.effects[, 1L])

      } else if (inputs$args$bias.type == "new_latent") {

        oldZ2 <- rbind(inputs$tmb.data$Z2_PO_pres, inputs$tmb.data$Z2_PO_quad)[inputs$po.info[ , "row.id"], ]
        res$fitted.values <- as.vector(oldX %*% res$fixed.effects[, 1L]) + as.vector(oldZ2 %*% res$bias.effects[, 1L]) + as.vector(oldZ %*% res$random.effects[, 1L])

      }
      # add in the fitted values on the presence/absence data if the model is an IDM
      if (inputs$args$data.type == "IDM") {
        attr(res$fitted.values, "PA") <- as.vector(inputs$tmb.data$X_PA %*% res$fixed.effects[, 1L]) + as.vector(inputs$tmb.data$Z_PA %*% res$random.effects[, 1L])
      }
    }
  } else { # in the case of a model without latent effects
    res$fitted.values <- NA # need to update
  }
  ##############################################################################

  # final additions for interfacing with other functions
  res$basis.per.res <- inputs$tmb.data$bf_per_res
  res$basis.functions <- inputs$basis.functions
  res$basis.fn.info <- inputs$bf.info
  res$approx.type <- inputs$args$model.type
  res$starting.pars <- inputs$args$starting.pars
  res$data <- inputs$args$data
  if (inputs$args$data.type == "IDM") {
    attr(res$data, "PA") <- inputs$args$IDM.presence.absence.df
  }
  res$formula <- inputs$args$formula
  if (inputs$args$data.type == "IDM") {
    attr(res$formula, "bias") <- inputs$args$bias.formula
  }
  res$coord.names <- inputs$args$coord.names
  res$quad.weights.name <- inputs$args$quad.weights.name
  res$pt.quad.id <- inputs$po.info[ , "pt.quad.id"]
  res$data.model.type <- inputs$args$data.type
  res$bf.matrix.type <- inputs$args$bf.matrix.type
  res$bias.type <- inputs$args$bias.type
  class(res) <- "scampr"

  return(res)

}
