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
#' @param in.parallel a logical indicating whether to calculate the spatial cross validation on parallel cores (only relevant if po.fold.id or/and pa.fold.id are supplied)
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
basis.search <- function(object, po.fold.id, pa.fold.id, basis.functions.list, max.basis.functions, radius.type = c("diag", "limiting"), bf.matrix.type = c("sparse", "dense"), domain.data, trunc.pa.prob = 1e-7, in.parallel = FALSE) {

  # checks not covered by model fitting
  radius.type <- match.arg(radius.type)
  bf.matrix.type <- match.arg(bf.matrix.type)
  # set switch for parallel computing
  if (in.parallel) {
    para.switch <- "yes"
  } else {
    para.switch <- "no"
  }

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
  fitted.ll <- c()
  fitted.ll_po <- c()
  fitted.ll_pa <- c()
  logloss <- c()
  auc <- c()
  predicted.ll <- c()
  ks <- c()
  k_bias <- c()
  timing <- c()
  flag_fit <- c()

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
        tmp.cv <- switch(para.switch,
                         no = spatial.kfold.cv(base.model, po.fold.id, pa.fold.id),
                         yes = spatial.kfold.cv_parallel(base.model, po.fold.id, pa.fold.id)
        )
        logloss[1] <- tmp.cv$predicted.cll.pa
        predicted.ll[1] <- tmp.cv$predicted.cll.po
        auc[1] <- tmp.cv$auc
      }
      # record the current k
      ks[1] <- 0
      # record the bias field k
      k_bias[1] <- 0
      # record the fitted loglikelihood
      fitted.ll[1] <- logLik(base.model)
      fitted.ll_po[1] <- base.model$ll.components$LL_PO_pres + base.model$ll.components$LL_PO_quad + base.model$ll.components$LL_random
      fitted.ll_pa[1] <- base.model$ll.components$LL_PA + base.model$ll.components$LL_random
      # record the timing
      timing[1] <- sum(base.model$cpu)
      # record the convergence flag
      flag_fit[1] <- base.model$se.flag != 0 | base.model$convergence != 0

      # reset the 'include.sre' flag
      call.list$include.sre <- quote(TRUE)

      # loop through the basis configurations to get the model fits
      for (config in 1:length(basis.functions.list)) {
        for (config2 in 1:length(basis.functions.list)) {

          # set the storage index for the nested loop
          lcv <- config2 + ((config - 1) * length(basis.functions.list))

          # skip the model fit if the configurations are the same
          if (config == config2) {

            # adjust the basis functions
            call.list$basis.functions <- basis.functions.list[[config]]
            call.list$po.biasing.basis.functions <- basis.functions.list[[config2]]
            # fit the model
            tmp.mod <- do.call("scampr", call.list)
            # record the fitted loglikelihood
            fitted.ll[lcv + 1] <- NA
            fitted.ll_po[lcv + 1] <- NA
            fitted.ll_pa[lcv + 1] <- NA
            # record the timing
            timing[lcv + 1] <- 0
            # record the convergence flag
            flag_fit[lcv + 1] <- 0
            # record the number of basis functions
            ks[lcv + 1] <- nrow(basis.functions.list[[config]])
            # record the bias field k
            k_bias[lcv + 1] <- nrow(basis.functions.list[[config2]])
            # if spatial folds are provided, perform cross validation
            if (!missing(po.fold.id) & !missing(pa.fold.id)) {
              logloss[lcv + 1] <- NA
              predicted.ll[lcv + 1] <- NA
              auc[lcv + 1] <- NA
            }
          } else {

            # adjust the basis functions
            call.list$basis.functions <- basis.functions.list[[config]]
            call.list$po.biasing.basis.functions <- basis.functions.list[[config2]]
            # fit the model
            tmp.mod <- do.call("scampr", call.list)
            # record the fitted loglikelihood
            fitted.ll[lcv + 1] <- logLik(tmp.mod)
            fitted.ll_po[lcv + 1] <- tmp.mod$ll.components$LL_PO_pres + tmp.mod$ll.components$LL_PO_quad + tmp.mod$ll.components$LL_random
            fitted.ll_pa[lcv + 1] <- tmp.mod$ll.components$LL_PA + tmp.mod$ll.components$LL_random
            # record the timing
            timing[lcv + 1] <- sum(tmp.mod$cpu)
            # record the convergence flag
            flag_fit[lcv + 1] <- tmp.mod$se.flag != 0 | tmp.mod$convergence != 0
            # record the number of basis functions
            ks[lcv + 1] <- nrow(basis.functions.list[[config]])
            # record the bias field k
            k_bias[lcv + 1] <- nrow(basis.functions.list[[config2]])
            # if spatial folds are provided, perform cross validation
            if (!missing(po.fold.id) & !missing(pa.fold.id)) {
              tmp.cv <- switch(para.switch,
                               no = spatial.kfold.cv(tmp.mod, po.fold.id, pa.fold.id),
                               yes = spatial.kfold.cv_parallel(tmp.mod, po.fold.id, pa.fold.id)
              )
              logloss[lcv + 1] <- tmp.cv$predicted.cll.pa
              predicted.ll[lcv + 1] <- tmp.cv$predicted.cll.po
              auc[lcv + 1] <- tmp.cv$auc
            }
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
        tmp.cv <- switch(para.switch,
                         no = spatial.kfold.cv(base.model, po.fold.id, pa.fold.id),
                         yes = spatial.kfold.cv_parallel(base.model, po.fold.id, pa.fold.id)
        )
        logloss[1] <- tmp.cv$predicted.cll.pa
        predicted.ll[1] <- tmp.cv$predicted.cll.po
        auc[1] <- tmp.cv$auc
      }
      # record the current k
      ks[1] <- 0
      # record the fitted loglikelihood
      fitted.ll[1] <- logLik(base.model)
      fitted.ll_po[1] <- base.model$ll.components$LL_PO_pres + base.model$ll.components$LL_PO_quad + base.model$ll.components$LL_random
      fitted.ll_pa[1] <- base.model$ll.components$LL_PA + base.model$ll.components$LL_random
      # record the timing
      timing[1] <- sum(base.model$cpu)
      # record the convergence flag
      flag_fit[1] <- base.model$se.flag != 0 | base.model$convergence != 0

      # reset the 'include.sre' flag
      call.list$include.sre <- quote(TRUE)

      # loop through the basis configurations to get the model fits
      for (config in 1:length(basis.functions.list)) {
        # adjust the basis functions
        call.list$basis.functions <- basis.functions.list[[config]]
        # fit the model
        tmp.mod <- do.call("scampr", call.list)
        # record the fitted loglikelihood
        fitted.ll[config + 1] <- logLik(tmp.mod)
        fitted.ll_po[config + 1] <- tmp.mod$ll.components$LL_PO_pres + tmp.mod$ll.components$LL_PO_quad + tmp.mod$ll.components$LL_random
        fitted.ll_pa[config + 1] <- tmp.mod$ll.components$LL_PA + tmp.mod$ll.components$LL_random
        # record the timing
        timing[config + 1] <- sum(tmp.mod$cpu)
        # record the convergence flag
        flag_fit[config + 1] <- tmp.mod$se.flag != 0 | tmp.mod$convergence != 0
        # record the number of basis functions
        ks[config + 1] <- nrow(basis.functions.list[[config]])
        # if spatial folds are provided, perform cross validation
        if (!missing(po.fold.id) & !missing(pa.fold.id)) {
          tmp.cv <- switch(para.switch,
                           no = spatial.kfold.cv(tmp.mod, po.fold.id, pa.fold.id),
                           yes = spatial.kfold.cv_parallel(tmp.mod, po.fold.id, pa.fold.id)
          )
          logloss[config + 1] <- tmp.cv$predicted.cll.pa
          predicted.ll[config + 1] <- tmp.cv$predicted.cll.po
          auc[config + 1] <- tmp.cv$auc
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
      tmp.cv <- switch(para.switch,
                       no = spatial.kfold.cv(base.model, po.fold.id),
                       yes = spatial.kfold.cv_parallel(base.model, po.fold.id)
      )
      predicted.ll[1] <- tmp.cv$predicted.cll.po
    }
    # record the current k
    ks[1] <- 0
    # record the fitted loglikelihood
    fitted.ll[1] <- logLik(base.model)
    fitted.ll_po[1] <- base.model$ll.components$LL_PO_pres + base.model$ll.components$LL_PO_quad + base.model$ll.components$LL_random
    fitted.ll_pa[1] <- base.model$ll.components$LL_PA + base.model$ll.components$LL_random
    # record the timing
    timing[1] <- sum(base.model$cpu)
    # record the convergence flag
    flag_fit[1] <- base.model$se.flag != 0 | base.model$convergence != 0

    # reset the 'include.sre' flag
    call.list$include.sre <- quote(TRUE)

    # loop through the basis configurations to get the model fits
    for (config in 1:length(basis.functions.list)) {
      # adjust the basis functions
      call.list$basis.functions <- basis.functions.list[[config]]
      # fit the model
      tmp.mod <- do.call("scampr", call.list)
      # record the fitted loglikelihood
      fitted.ll[config + 1] <- logLik(tmp.mod)
      fitted.ll_po[config + 1] <- tmp.mod$ll.components$LL_PO_pres + tmp.mod$ll.components$LL_PO_quad + tmp.mod$ll.components$LL_random
      fitted.ll_pa[config + 1] <- tmp.mod$ll.components$LL_PA + tmp.mod$ll.components$LL_random
      # record the timing
      timing[config + 1] <- sum(tmp.mod$cpu)
      # record the convergence flag
      flag_fit[config + 1] <- tmp.mod$se.flag != 0 | tmp.mod$convergence != 0
      # record the number of basis functions
      ks[config + 1] <- nrow(basis.functions.list[[config]])
      # if spatial fold is provided, perform cross validation
      if (!missing(po.fold.id)) {
        tmp.cv <- switch(para.switch,
                         no = spatial.kfold.cv(tmp.mod, po.fold.id),
                         yes = spatial.kfold.cv(tmp.mod, po.fold.id)
        )
        predicted.ll[config + 1] <- tmp.cv$predicted.cll.po
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
      tmp.cv <- switch(para.switch,
                       no = spatial.kfold.cv(base.model, pa.fold.id = pa.fold.id),
                       yes = spatial.kfold.cv_parallel(base.model, pa.fold.id = pa.fold.id)
      )
      logloss[1] <- tmp.cv$predicted.cll.pa
      # predicted.ll <- c(predicted.ll, tmp.cv$predicted.cll.po)
      auc[1] <- tmp.cv$auc
    }
    # record the current k
    ks[1] <- 0
    # record the fitted loglikelihood
    fitted.ll[1] <- logLik(base.model)
    fitted.ll_po[1] <- base.model$ll.components$LL_PO_pres + base.model$ll.components$LL_PO_quad + base.model$ll.components$LL_random
    fitted.ll_pa[1] <- base.model$ll.components$LL_PA + base.model$ll.components$LL_random
    # record the timing
    timing[1] <- sum(base.model$cpu)
    # record the convergence flag
    flag_fit[1] <- base.model$se.flag != 0 | base.model$convergence != 0

    # reset the 'include.sre' flag
    call.list$include.sre <- quote(TRUE)

    # loop through the basis configurations to get the model fits
    for (config in 1:length(basis.functions.list)) {
      # adjust the basis functions
      call.list$basis.functions <- basis.functions.list[[config]]
      # fit the model
      tmp.mod <- do.call("scampr", call.list)
      # record the fitted loglikelihood
      fitted.ll[config + 1] <- logLik(tmp.mod)
      fitted.ll_po[config + 1] <- tmp.mod$ll.components$LL_PO_pres + tmp.mod$ll.components$LL_PO_quad + tmp.mod$ll.components$LL_random
      fitted.ll_pa[config + 1] <- tmp.mod$ll.components$LL_PA + tmp.mod$ll.components$LL_random
      # record the timing
      timing[config + 1] <- sum(tmp.mod$cpu)
      # record the convergence flag
      flag_fit[config + 1] <- tmp.mod$se.flag != 0 | tmp.mod$convergence != 0
      # record the number of basis functions
      ks[config + 1] <- nrow(basis.functions.list[[config]])
      # if spatial folds are provided, perform cross validation
      if (!missing(pa.fold.id)) {
        tmp.cv <- switch(para.switch,
                         no = spatial.kfold.cv(tmp.mod, po.fold.id, pa.fold.id),
                         yes=  spatial.kfold.cv_parallel(tmp.mod, po.fold.id, pa.fold.id)
        )
        logloss[config + 1] <- tmp.cv$predicted.cll.pa
        auc[config + 1] <- tmp.cv$auc
      }
    }
  } else {
    stop(paste0("model type not reognised: ", object$model.type))
  }

  ##############################################################################
  # Collate the results #
  ##############################################################################

  if (object$model.type == "IDM" & object$random.bias.type == "field2") {
    res <- data.frame(cbind(k = ks[-1], k_bias = k_bias[-1], fitted.ll = fitted.ll[-1], fitted.ll_pa = fitted.ll_pa[-1], fitted.ll_po = fitted.ll_po[-1], predicted.ll_po = predicted.ll[-1], predicted.ll_pa = logloss[-1], AUC = auc[-1], cpu = timing[-1], convergence = flag_fit[-1]))
    # determine the maximum fitted log-likelihood (excluding fits with poor convergence)

    # with dense k and coarse k_bias (1)
    tmp.res1 <- res[res$k_bias < res$k, ]
    tmp.res1$fitted.ll[tmp.res1$convergence == 1] <- NA
    k1 <- tmp.res1$k[which.max(tmp.res1$fitted.ll)]
    k1_bias <- tmp.res1$k_bias[which.max(tmp.res1$fitted.ll)]
    best.bfs1 <- list(k = basis.functions.list[[match(k1, lapply(basis.functions.list, nrow))]],
                      k_bias = basis.functions.list[[match(k1_bias, lapply(basis.functions.list, nrow))]]
    )
    attr(best.bfs1, "fitted logLik") <- tmp.res1$fitted.ll[which.max(tmp.res1$fitted.ll)]
    attr(res, "basis.config: k dense k_bias coarse") <- best.bfs1

    # with coarse k and dense k_bias (2)
    tmp.res2 <- res[res$k_bias > res$k, ]
    tmp.res2$fitted.ll[tmp.res2$convergence == 1] <- NA
    k2 <- tmp.res2$k[which.max(tmp.res2$fitted.ll)]
    k2_bias <- tmp.res2$k_bias[which.max(tmp.res2$fitted.ll)]
    best.bfs2 <- list(k = basis.functions.list[[match(k2, lapply(basis.functions.list, nrow))]],
                      k_bias = basis.functions.list[[match(k2_bias, lapply(basis.functions.list, nrow))]]
    )
    attr(best.bfs2, "fitted logLik") <- tmp.res2$fitted.ll[which.max(tmp.res2$fitted.ll)]
    attr(res, "basis.config: k coarse k_bias dense") <- best.bfs2

    # add on the basis function list searched
    attr(res, "basis.functions.list") <- basis.functions.list
    # add in the baseline (no SRE) scenario
    attr(res, "baseline") <- data.frame(cbind(k = ks[1], k_bias = k_bias[1], fitted.ll = fitted.ll[1], fitted.ll_pa = fitted.ll_pa[1], fitted.ll_po = fitted.ll_po[1], predicted.ll_po = predicted.ll[1], predicted.ll_pa = logloss[1], AUC = auc[1], cpu = timing[1], convergence = flag_fit[1]))
  } else {
    res <- data.frame(cbind(k = ks, k_bias = k_bias, fitted.ll = fitted.ll, fitted.ll_pa = fitted.ll_pa, fitted.ll_po = fitted.ll_po, predicted.ll_po = predicted.ll, predicted.ll_pa = logloss, AUC = auc, cpu = timing, convergence = flag_fit))
    # determine the maximum fitted log-likelihood (excluding fits with poor convergence)
    tmp.res <- res
    tmp.res$fitted.ll[tmp.res$convergence == 1] <- NA
    best.bfs <- basis.functions.list[[which.max(tmp.res$fitted.ll) - 1]] # minus one since the search results include the NULL model in this case
    attr(best.bfs, "fitted logLik") <- tmp.res$fitted.ll[which.max(tmp.res$fitted.ll)]
    attr(res, "basis.config") <- best.bfs
    # add on the basis function list searched
    attr(res, "basis.functions.list") <- basis.functions.list
  }

  # add a secondary class for default plotting, etc.
  class(res) <- c(class(res), "basis.search")

  return(res)
}
