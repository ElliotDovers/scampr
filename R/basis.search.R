#' Search algorithm for simple 2D basis functions configurations on scampr models
#'
#' @description This function takes in a scampr model and calculates likelihoods and AIC for the list of basis functions supplied. If none are supplied then the algorithm fits increasingly dense regular grids of basis functions (of the type created by scampr::simple_basis). The algorithm starts with an IPP (i.e. zero basis functions) and increases to <= 'max.basis.functions'. If 'po.fold.id' and/or 'pa.fold.id' are supplied then the function will perform a k-fold hold-one-out cross validation to calculate out-of-sample likelihoods (conditional on the latent field).
#'
#' @param object a scampr model: object of class 'scampr' that provides the framework for the search algorithm. Recommended that an IPP model of the appropriate type is used.
#' @param po.fold.id Optional for cross-validation. An integer or factor vector the same length as the po.data that describes the CV fold that each location falls into.
#' @param pa.fold.id Optional for cross-validation. An integer or factor vector the same length as the pa.data that describes the CV fold that each location falls into.
#' @param basis.functions.list Optional. A list of basis function configuration to trial. See \code{?scampr} for details on accepted types.
#' @param max.basis.functions Optional. An integer describing a rough upper limit to the number of basis functions to search. Defaults to half the number of presences in the data sets.
#' @param radius.type a character string describing the type of radius length to use. One of 'diag' = diagonal dist. between nodes or 'limiting' = sqrt(Domain Area)/log(k).
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to use sparse or dense matrix computations for the basis functions created.
#' @param domain.data Optional. A data frame of columns 'coord.names' that contains at least the extremities of the domain of interest. Useful to ensure the same basis function configurations are created by 'simple_basis' if comparing to various searches.
#' @param trunc.pa.prob Optional. A small positive number by which the predicted probability of presence is truncated. This can be used to ensure infinite values are avoiding within the cross-validation.
#'
#' @return a data.frame with columns including- 'nodes.on.long.edge': number used in scampr::simple_basis to create basis configuration. 'bf': the number of basis functions. 'loglik': the fitting marginal log-likelihood. 'aic': the corresponding AIC. Optionally, 'predicted_cll_po': the conditional (on the latent field) Presence-only likelihood. 'predicted_cll_pa': the conditional (on the latent field) Presence/Absence likelihood. 'roc_auc': Area under the ROC curve on the Presence/Absence data. Optional columns are the results from a cross-validation described by 'po.fold.id' and/or 'pa.fold.id'. (_va or _lp subscript for approx. type if both are calculated).
#' @export
#'
#' @examples
#' #' # Get the gorilla nesting data
#' dat <- gorillas
#'
#' # Standardise the elevation covariate
#' dat$elev.std <- scale(dat$elevation)
#'
#' # Fit an IPP model to the point pattern
#' m.ipp <- scampr(pres ~ elev.std, data = dat, model.type = "ipp")
#'  \dontrun{
#' # Search through an increasingly dense regular grid of basis functions
#' res <- simple_basis_search(m.ipp)
#' }
basis.search <- function(object, po.fold.id, pa.fold.id, basis.functions.list, max.basis.functions, radius.type = c("diag", "limiting"), bf.matrix.type = c("sparse", "dense"), domain.data, trunc.pa.prob = 1e-7) {

  # checks not covered by model fitting
  radius.type <- match.arg(radius.type)
  bf.matrix.type <- match.arg(bf.matrix.type)

  # Use provided model data as domain.data if missing and check coords are present
  if (missing(domain.data)) {
    domain.data <- rbind(object$data[ , object$coord.names], attr(object$data, "PA")[ , object$coord.names])
  } else {
    if (!all(object$coord.names %in% colnames(domain.data))) {
      stop(paste0("Model coord.names, ", object$coord.names, ", not found in 'domain.data' frame provided"))
    }
  }

  ##############################################################################
  # Decide on the basis functions to be trialed #
  ##############################################################################

  if (missing(basis.functions.list)) { # if a list of basis configurations is not provided
    if (missing(max.basis.functions)) { # if a maximum number of basis functions to try is not provided
      # set the max number of basis functions to half the number of presences
      max.basis.functions <- 0.5 * max(c(sum(attr(object$data, "PA")[ , all.vars(object$formula[[2]])] > 0), sum(object$data[, all.vars(object$formula[[2]])])))
    }
    # initialise objects
    tmp.nodes <- NULL
    tmp.k <- NULL
    basis.functions.list <- list()
    # first iteration
    tmp.nodes <- c(tmp.nodes, 2)
    # simple basis function configuration
    tmp.bfs <- simple_basis(tmp.nodes[1], data = domain.data, radius.type = radius.type)
    basis.functions.list[[1]] <- tmp.bfs
    tmp.nodes <- c(tmp.nodes, 3)
    tmp.k <- c(tmp.k, nrow(tmp.bfs))
    counter <- 2
    while (nrow(tmp.bfs) < max.basis.functions) {
      # create new basis config.
      tmp.bfs <- simple_basis(tmp.nodes[counter], data = domain.data, radius.type = radius.type)
      # store in basis function list
      basis.functions.list[[counter]] <- tmp.bfs
      # record the number of basis functions
      tmp.k <- c(tmp.k, nrow(tmp.bfs))

      # set the next number of nodes to use
      tmp.nodes <- c(tmp.nodes, tmp.nodes[[counter]] + 1)
      # increase the counter (loop control variable)
      counter <- counter + 1
    }
  }

  ##############################################################################
  # Loop through basis function list and calculate required metrics #
  ##############################################################################

  # initialise the result objects
  fitted.ll <- NULL
  logloss <- NULL
  auc <- NULL
  predicted.ll <- NULL
  ks <- 0
  k_bias <- NULL
  timing <- NULL
  flag_fit <- NULL

  ## FOR AN IDM
  if (object$model.type == "IDM") {

    # If the model accounts for bias using a second latent field then there is an extra search loop
    if (object$random.bias.type == "field2") {
      # get the model's call
      call.list <- as.list(object$call)
      # adjust the call to not include SRE
      call.list$include.sre <- quote(FALSE)
      # fit the model
      call.list[[1]] <- NULL
      base.model <- do.call("scampr", call.list)

      # if spatial folds are provided, perform cross validation
      if (!missing(po.fold.id) & !missing(pa.fold.id)) {
        tmp.cv <- spatial.kfold.cv(base.model, po.fold.id, pa.fold.id)
        logloss <- c(logloss, tmp.cv$predicted.cll.pa)
        predicted.ll <- c(predicted.ll, tmp.cv$predicted.cll.po)
        auc <- c(auc, tmp.cv$auc)
      }

      # attach the bias field k
      k_bias <- c(k_bias, 0)
      # attach the fitted loglikelihood
      fitted.ll <- c(fitted.ll, logLik(base.model))
      # attach the timing
      timing <- c(timing, sum(base.model$cpu))
      # attach the convergence flag
      flag_fit <- c(flag_fit, base.model$se.flag != 0 | base.model$convergence != 0)

      # reset the 'include.sre' flag
      call.list$include.sre <- quote(TRUE)

      # loop through the basis configurations to get the model fits
      for (config in 1:length(basis.functions.list)) {
        for (config2 in 1:length(basis.functions.list)) {
          # adjust the basis functions
          call.list$basis.functions <- basis.functions.list[[config]]
          call.list$po.biasing.basis.functions <- basis.functions.list[[config2]]
          # fit the model
          tmp.mod <- do.call("scampr", call.list)
          # record the fitted loglikelihood
          fitted.ll <- c(fitted.ll, logLik(tmp.mod))
          # attach the timing
          timing <- c(timing, sum(tmp.mod$cpu))
          # attach the convergence flag
          flag_fit <- c(flag_fit, tmp.mod$se.flag != 0 | tmp.mod$convergence != 0)
          # record the number of basis functions
          ks <- c(ks, nrow(basis.functions.list[[config]]))
          # attach the bias field k
          k_bias <- c(k_bias, nrow(basis.functions.list[[config2]]))
          # if spatial folds are provided, perform cross validation
          if (!missing(po.fold.id) & !missing(pa.fold.id)) {
            tmp.cv <- spatial.kfold.cv(tmp.mod, po.fold.id, pa.fold.id)
            logloss <- c(logloss, tmp.cv$predicted.cll.pa)
            predicted.ll <- c(predicted.ll, tmp.cv$predicted.cll.po)
            auc <- c(auc, tmp.cv$auc)
          }
        }
      }
    } else {
      # get the model's call
      call.list <- as.list(object$call)
      # adjust the call to not include SRE
      call.list$include.sre <- quote(FALSE)
      # fit the model
      call.list[[1]] <- NULL
      base.model <- do.call("scampr", call.list)

      # if spatial folds are provided, perform cross validation
      if (!missing(po.fold.id) & !missing(pa.fold.id)) {
        tmp.cv <- spatial.kfold.cv(base.model, po.fold.id, pa.fold.id)
        logloss <- c(logloss, tmp.cv$predicted.cll.pa)
        predicted.ll <- c(predicted.ll, tmp.cv$predicted.cll.po)
        auc <- c(auc, tmp.cv$auc)
      }

      # attach the fitted loglikelihood
      fitted.ll <- c(fitted.ll, logLik(base.model))
      # attach the timing
      timing <- c(timing, sum(base.model$cpu))
      # attach the convergence flag
      flag_fit <- c(flag_fit, base.model$se.flag != 0 | base.model$convergence != 0)

      # reset the 'include.sre' flag
      call.list$include.sre <- quote(TRUE)

      # loop through the basis configurations to get the model fits
      for (config in 1:length(basis.functions.list)) {
        # adjust the basis functions
        call.list$basis.functions <- basis.functions.list[[config]]
        # fit the model
        tmp.mod <- do.call("scampr", call.list)
        # record the fitted loglikelihood
        fitted.ll <- c(fitted.ll, logLik(tmp.mod))
        # attach the timing
        timing <- c(timing, sum(tmp.mod$cpu))
        # attach the convergence flag
        flag_fit <- c(flag_fit, tmp.mod$se.flag != 0 | tmp.mod$convergence != 0)
        # record the number of basis functions
        ks <- c(ks, nrow(basis.functions.list[[config]]))
        # if spatial folds are provided, perform cross validation
        if (!missing(po.fold.id) & !missing(pa.fold.id)) {
          tmp.cv <- spatial.kfold.cv(tmp.mod, po.fold.id, pa.fold.id)
          logloss <- c(logloss, tmp.cv$predicted.cll.pa)
          predicted.ll <- c(predicted.ll, tmp.cv$predicted.cll.po)
          auc <- c(auc, tmp.cv$auc)
        }
      }
    }

  ## FOR A PPM
  } else  if (object$model.type == "PO") {

    # get the model's call
    call.list <- as.list(object$call)
    # adjust the call to not include SRE
    call.list$include.sre <- quote(FALSE)
    # fit the model
    call.list[[1]] <- NULL
    base.model <- do.call("scampr", call.list)

    # if spatial fold is provided, perform cross validation
    if (!missing(po.fold.id)) {
      tmp.cv <- spatial.kfold.cv(base.model, po.fold.id)
      # logloss <- c(logloss, tmp.cv$predicted.cll.pa)
      predicted.ll <- c(predicted.ll, tmp.cv$predicted.cll.po)
      # auc <- c(auc, tmp.cv$auc)
    }

    # attach the fitted loglikelihood
    fitted.ll <- c(fitted.ll, logLik(base.model))
    # attach the timing
    timing <- c(timing, sum(base.model$cpu))
    # attach the convergence flag
    flag_fit <- c(flag_fit, base.model$se.flag != 0 | base.model$convergence != 0)

    # reset the 'include.sre' flag
    call.list$include.sre <- quote(TRUE)

    # loop through the basis configurations to get the model fits
    for (config in 1:length(basis.functions.list)) {
      # adjust the basis functions
      call.list$basis.functions <- basis.functions.list[[config]]
      # fit the model
      tmp.mod <- do.call("scampr", call.list)
      # record the fitted loglikelihood
      fitted.ll <- c(fitted.ll, logLik(tmp.mod))
      # attach the timing
      timing <- c(timing, sum(tmp.mod$cpu))
      # attach the convergence flag
      flag_fit <- c(flag_fit, tmp.mod$se.flag != 0 | tmp.mod$convergence != 0)
      # record the number of basis functions
      ks <- c(ks, nrow(basis.functions.list[[config]]))
      # if spatial fold is provided, perform cross validation
      if (!missing(po.fold.id)) {
        tmp.cv <- spatial.kfold.cv(tmp.mod, po.fold.id)
        # logloss <- c(logloss, tmp.cv$predicted.cll.pa)
        predicted.ll <- c(predicted.ll, tmp.cv$predicted.cll.po)
        # auc <- c(auc, tmp.cv$predicted.cll.po)
      }
    }

  ## FOR A CLOGLOG BINOMIAL
  } else if (object$model.type == "PA") {

    # get the model's call
    call.list <- as.list(object$call)
    # adjust the call to not include SRE
    call.list$include.sre <- quote(FALSE)
    # fit the model
    call.list[[1]] <- NULL
    base.model <- do.call("scampr", call.list)

    # if spatial folds are provided, perform cross validation
    if (!missing(pa.fold.id)) {
      tmp.cv <- spatial.kfold.cv(base.model, pa.fold.id = pa.fold.id)
      logloss <- c(logloss, tmp.cv$predicted.cll.pa)
      # predicted.ll <- c(predicted.ll, tmp.cv$predicted.cll.po)
      auc <- c(auc, tmp.cv$auc)
    }

    # attach the fitted loglikelihood
    fitted.ll <- c(fitted.ll, logLik(base.model))
    # attach the timing
    timing <- c(timing, sum(base.model$cpu))
    # attach the convergence flag
    flag_fit <- c(flag_fit, base.model$se.flag != 0 | base.model$convergence != 0)

    # reset the 'include.sre' flag
    call.list$include.sre <- quote(TRUE)

    # loop through the basis configurations to get the model fits
    for (config in 1:length(basis.functions.list)) {
      # adjust the basis functions
      call.list$basis.functions <- basis.functions.list[[config]]
      # fit the model
      tmp.mod <- do.call("scampr", call.list)
      # record the fitted loglikelihood
      fitted.ll <- c(fitted.ll, logLik(tmp.mod))
      # attach the timing
      timing <- c(timing, sum(tmp.mod$cpu))
      # attach the convergence flag
      flag_fit <- c(flag_fit, tmp.mod$se.flag != 0 | tmp.mod$convergence != 0)
      # record the number of basis functions
      ks <- c(ks, nrow(basis.functions.list[[config]]))
      # if spatial folds are provided, perform cross validation
      if (!missing(pa.fold.id)) {
        tmp.cv <- spatial.kfold.cv(tmp.mod, po.fold.id, pa.fold.id)
        logloss <- c(logloss, tmp.cv$predicted.cll.pa)
        # predicted.ll <- c(predicted.ll, tmp.cv$predicted.cll.po)
        auc <- c(auc, tmp.cv$predicted.cll.po)
      }
    }
  } else {
    stop(paste0("model type not reognised: ", object$model.type))
  }

  ##############################################################################
  # Collate the results #
  ##############################################################################

  if (object$model.type == "IDM" & object$random.bias.type == "field2") {
    res <- data.frame(cbind(k = ks[-1], k_bias = k_bias[-1], fitted.ll = fitted.ll[-1], predicted.ll_po = predicted.ll[-1], predicted.ll_pa = logloss[-1], AUC = auc[-1], cpu = timing[-1], convergence = flag_fit[-1]))
    attr(res, "basis.functions.list") <- basis.functions.list
    attr(res, "baseline") <- data.frame(cbind(k = ks[1], k_bias = k_bias[1], fitted.ll = fitted.ll[1], predicted.ll_po = predicted.ll[1], predicted.ll_pa = logloss[1], AUC = auc[1], cpu = timing[-1], convergence = flag_fit[-1]))
  } else {
    res <- data.frame(cbind(k = ks, k_bias = k_bias, fitted.ll = fitted.ll, predicted.ll_po = predicted.ll, predicted.ll_pa = logloss, AUC = auc, cpu = timing, convergence = flag_fit))
    attr(res, "basis.functions.list") <- basis.functions.list
  }

  # add a secondary class for default plotting, etc.
  class(res) <- c(class(res), "basis.search")

  return(res)
}
