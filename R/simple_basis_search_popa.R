#' Internal Basis Function Search Algorithm for POPA models
#'
#' @param po.formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the Presence-only data model to be fitted. The 'response' must be a binary that indicates whether a datum is a presence or quadrature point. See GLM function for further formula details.
#' @param pa.formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the Presence/Absence data model to be fitted. All predictor terms must also be included in po.formula. The 'response' must be a must be either the site abundance or a binary indicating whether there is a presence or not. See GLM function for further formula details.
#' @param po.data a data frame containing predictors at both presence-records and quadrature as well as the po.formula 'response'.
#' @param pa.data a data frame containing predictors and response for the pa.formula.
#' @param po.fold.id Optional for cross-validation. An integer or factor vector the same length as the po.data that describes the CV fold that each location falls into.
#' @param pa.fold.id Optional for cross-validation. An integer or factor vector the same length as the pa.data that describes the CV fold that each location falls into.
#' @param max.basis.functions Optional. An integer describing a rough upper limit to the number of basis functions to search. Defaults to half the number of presences in the data sets.
#' @param coord.names a vector of character strings describing the column names of the coordinates in both data frames.
#' @param quad.weights.name a charater string of the column name of quadrature weights in the po.data.
#' @param radius.type a character string describing the type of radius length to use. One of 'diag' = diagonal dist. between nodes or 'limiting' = sqrt(Domain Area)/log(k).
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to use sparse or dense matrix computations for the basis functions created.
#' @param domain.data Optional. A data frame of columns 'coord.names' that contains at least the extremities of the domain of interest. Useful to ensure the same basis function configurations are created by 'simple_basis' if comparing to various searches.
#' @param approx.with a charater string indicating the type of approximation to use for the intractible marginalisation. One of 'variational' or 'laplace'.
#' @param trunc.pa.prob Optional. A small positve number by which the predicted probability of presence is truncated. This can be used to ensure infinite values are avoiding within the cross-validation.
#'
#' @return a data.frame with columns including- 'nodes.on.long.edge': number used in scampr::simple_basis to create basis configuration. 'bf': the number of basis functions. 'loglik': the fitting marginal log-likelihood. 'aic': the corresponding AIC. Optionally, 'predicted_cll_po': the conditional (on the latent field) Presence-only likelihood. 'predicted_cll_pa': the conditional (on the latent field) Presence/Absence likelihood. 'roc_auc': Area under the ROC curve on the Presence/Absence data. Optional columns are the results from a cross-validation described by 'po.fold.id' and/or 'pa.fold.id'. (_va or _lp subscript for approx. type if both are calculated)
#' @noRd
#'
#' @importFrom pROC roc auc
#'
#' @examples
#' #' # Get the Eucalypt data
#' dat_po <- eucalypt[["po"]]
#' dat_pa <- eucalypt[["pa"]]
#'
#' \dontrun{
#' # Fit with a shared latent field
#' res <- simple_basis_search_popa(pres ~ TMP_MIN + D_MAIN_RDS, Y ~ TMP_MIN,
#' po.data = dat_po, pa.data = dat_pa)
#' }
simple_basis_search_popa <- function(po.formula, pa.formula, po.data, pa.data, po.fold.id, pa.fold.id, max.basis.functions, coord.names = c("x", "y"), quad.weights.name = "quad.size", radius.type = c("diag", "limiting"), bf.matrix.type = c("sparse", "dense"), domain.data, approx.with = c("laplace", "variational"), trunc.pa.prob = 1e-7) {

  # Use provided model data as domain.data if missing and check coords are present
  if (missing(domain.data)) {
    if (!all(coord.names %in% colnames(pa.data)) & !all(coord.names %in% colnames(po.data))) {
      stop(paste0("coord.names, ", coord.names, ", not found in both data sets provided"))
    }
    domain.data <- rbind(pa.data[ , coord.names], po.data[ , coord.names])
  } else {
    if (!all(coord.names %in% colnames(domain.data))) {
      stop(paste0("coord.names, ", coord.names, ", not found in 'domain.data' frame provided"))
    }
  }
  # checks not covered by model fitting
  radius.type <- match.arg(radius.type)
  bf.matrix.type <- match.arg(bf.matrix.type)
  approx.with <- match.arg(approx.with)
  if (approx.with == "variational") {
    stop("variational approx. not available for POPA models")
  }
  # Assign the maximum number of basis functions to searched, if missing
  if(missing(max.basis.functions)) {
    #  set at half the maximum number of presences in either data set
    max.basis.functions <- 0.5 * max(c(sum(pa.data[ , all.vars(pa.formula[[2]])] > 0), sum(po.data[, all.vars(po.formula[[2]])])))
  }

  # If either CV folds vector is missing just perform in-sample likelihood search
  if (missing(po.fold.id) | missing(pa.fold.id)) {

    # initialise the storage objects
    nbf <- NULL
    loglik <- NULL
    aic <- NULL
    pred_loglik_po <- NULL # these will remain NULL in this case
    pred_loglik_pa <- NULL # these will remain NULL in this case
    roc_auc <- NULL # these will remain NULL in this case
    counter.counter <- NULL
    counter <- 1
    counter.counter <- c(counter.counter, counter)

    # First iteration hard coded for IPP model, i.e. zero basis functions #

    m_ipp <- popa(po.formula, pa.formula, po.data, pa.data, model.type = "ipp", coord.names = coord.names, quad.weights.name = quad.weights.name, bf.matrix.type = bf.matrix.type)

    # Store Results
    nbf[counter] <- 0
    loglik[counter] <- logLik.scampr(m_ipp)
    aic[counter] <- AIC.scampr(m_ipp)
    print(paste0("Completed fit with 0 basis functions (IPP)"))
    counter <- 2
    counter.counter <- c(counter.counter, counter)
    bfs <- cbind(NA, NA)

    # Start looping until we hit the max number of basis functions
    while (nrow(bfs) <= max.basis.functions) {
      # Set the current simple basis function configuration
      bfs <- simple_basis(counter, data = domain.data, radius.type = radius.type)
      # try to fit the full model
      m <- NULL
      try(assign("m", popa(po.formula, pa.formula, po.data, pa.data, simple.basis = bfs, starting.pars = m_ipp, coord.names = coord.names, quad.weights.name = quad.weights.name, bf.matrix.type = bf.matrix.type)))

      if (is.null(m)) {
        # Store missing results in this case
        loglik[counter] <- NA
        aic[counter] <- NA
      } else {
        # Store results
        loglik[counter] <-  logLik.scampr(m)
        aic[counter] <- AIC.scampr(m)
        m <- NULL
      }
      # Store common results
      nbf[counter] <- nrow(bfs)
      print(paste0("Completed fit with ", nrow(bfs), " basis functions"))
      counter <- counter + 1
      counter.counter <- c(counter.counter, counter)
    }
  } else { # IF CV BOTH FOLDS ARE PRESENT

    # Convert CV folds to numerics in case factor is provided
    po.fold.id <- as.numeric(po.fold.id)
    pa.fold.id <- as.numeric(pa.fold.id)

    # check the CV folds are matching numbers
    if (!identical(sort(unique(po.fold.id)), sort(unique(pa.fold.id)))) {
      stop("Both fold.id vectors must contain the same fold numbers")
    }

    # initialise the storage objects
    nbf <- NULL
    loglik <- NULL
    aic <- NULL
    pred_loglik_po <- NULL
    pred_loglik_pa <- NULL
    roc_auc <- NULL
    counter.counter <- NULL
    counter <- 1
    counter.counter <- c(counter.counter, counter)

    # First iteration hard coded for IPP model, i.e. zero basis functions #

    m_ipp <- popa(po.formula, pa.formula, po.data, pa.data, model.type = "ipp", coord.names = coord.names, quad.weights.name = quad.weights.name, bf.matrix.type = bf.matrix.type)

    # Initialise Objects #
    train_mods_ipp <- list() # training models
    test.abund <- NULL # predicted abundance rate
    test.inten <- NULL # predicted intensity rate
    pa.resp <- all.vars(pa.formula[[2]]) # PA response
    po.resp <- all.vars(po.formula[[2]]) # PO response
    test.pa.resp <- NULL # test PA response
    test.po.resp <- NULL # test PO response
    test.po.quad.sizes <- NULL
    # Loop through CV folds
    for (i in sort(unique(po.fold.id))) {
      train_mods_ipp[[i]] <- popa(po.formula, pa.formula, po.data[po.fold.id != i, ], pa.data[pa.fold.id != i, ], model.type = "ipp", coord.names = coord.names, quad.weights.name = quad.weights.name, bf.matrix.type = bf.matrix.type)
      test.abund <- c(test.abund, predict.scampr(train_mods_ipp[[i]], newdata = pa.data[pa.fold.id == i, ], process = "abundance"))
      test.inten <- c(test.inten, predict.scampr(train_mods_ipp[[i]], newdata = po.data[po.fold.id == i, ]))
      test.pa.resp <- c(test.pa.resp, as.numeric(pa.data[pa.fold.id == i, pa.resp] > 0))
      test.po.resp <- c(test.po.resp, po.data[po.fold.id == i, po.resp])
      test.po.quad.sizes <- c(test.po.quad.sizes, po.data[po.fold.id == i, quad.weights.name])
    }
    # predicted presence probability
    pres_prob <- 1 - exp(-exp(test.abund))
    # NEED TO TRUNCATE PRESENCE PROBABILITY BY 'trunc.pa.prob'
    pres_prob[pres_prob <= trunc.pa.prob] <- trunc.pa.prob
    pres_prob[1 - pres_prob <= trunc.pa.prob] <- 1 - trunc.pa.prob
    # predicted conditional log-likelihood on the PA data
    pll_pa <- sum((test.pa.resp * log(pres_prob)) + ((1 - test.pa.resp) * log(1 - pres_prob)))
    # predicted conditional log-likelihood on the PO data
    pll_po <- sum(test.inten[test.po.resp == 1]) - sum(test.po.quad.sizes[test.po.resp == 0] * exp(test.inten[test.po.resp == 0]))
    # area under ROC curve
    roccurve <- pROC::roc(test.pa.resp, as.vector(pres_prob), quiet = T)
    auc_val <- as.numeric(pROC::auc(roccurve))

    # Store Results
    nbf[counter] <- 0
    loglik[counter] <- logLik.scampr(m_ipp)
    aic[counter] <- AIC.scampr(m_ipp)
    pred_loglik_po[counter] <- pll_po
    pred_loglik_pa[counter] <- pll_pa
    roc_auc[counter] <- auc_val
    print(paste0("Completed fits with 0 basis functions (IPP)"))
    counter <- 2
    counter.counter <- c(counter.counter, counter)
    bfs <- cbind(NA, NA)
    rm(pres_prob, pll_pa, pll_po, roccurve, auc_val)

    # Start looping until we hit the max number of basis functions
    while (nrow(bfs) <= max.basis.functions) {
      # Set the current simple basis function configuration
      bfs <- simple_basis(counter, data = domain.data, radius.type = radius.type)
      # try to fit the full model
      m <- NULL
      try(assign("m", popa(po.formula, pa.formula, po.data, pa.data, simple.basis = bfs, starting.pars = m_ipp, coord.names = coord.names, quad.weights.name = quad.weights.name, bf.matrix.type = bf.matrix.type)))
      # Initialise Objects #
      train_mods <- list() # training models
      test.abund <- NULL # predicted abundance rate
      test.inten <- NULL # predicted intensity rate
      tmp.m <- NULL
      # Loop through CV folds
      for (i in sort(unique(po.fold.id))) {
        try(assign("tmp.m", popa(po.formula, pa.formula, po.data[po.fold.id != i, ], pa.data[pa.fold.id != i, ], simple.basis = bfs, starting.pars =  train_mods_ipp[[i]], coord.names = coord.names, quad.weights.name = quad.weights.name, bf.matrix.type = bf.matrix.type)))
        # Only get predictions if the current model fits
        if (!is.null(tmp.m)) {
          train_mods[[i]] <- tmp.m
          test.abund <- c(test.abund, predict.scampr(train_mods[[i]], newdata = pa.data[pa.fold.id == i, ], process = "abundance"))
          test.inten <- c(test.inten, predict.scampr(train_mods[[i]], newdata = po.data[po.fold.id == i, ]))
          tmp.m <- NULL
        }
      }

      # Can only get the test results if ALL the training models fit
      if (length(train_mods) != length(unique(po.fold.id))) {
        # Store missing results in this case
        pred_loglik_po[counter] <- NA
        pred_loglik_pa[counter] <- NA
        roc_auc[counter] <- NA
      } else {
        # predicted presence probability
        pres_prob <- 1 - exp(-exp(test.abund))
        # NEED TO TRUNCATE PRESENCE PROBABILITY BY 'trunc.pa.prob'
        pres_prob[pres_prob <= trunc.pa.prob] <- trunc.pa.prob
        pres_prob[1 - pres_prob <= trunc.pa.prob] <- 1 - trunc.pa.prob
        # predicted conditional log-likelihood on the PA data
        pll_pa <- sum((test.pa.resp * log(pres_prob)) + ((1 - test.pa.resp) * log(1 - pres_prob)))
        # predicted conditional log-likelihood on the PO data
        pll_po <- sum(test.inten[test.po.resp == 1]) - sum(test.po.quad.sizes[test.po.resp == 0] * exp(test.inten[test.po.resp == 0]))
        # area under ROC curve
        roccurve <- pROC::roc(test.pa.resp, as.vector(pres_prob), quiet = T)
        auc_val <- as.numeric(pROC::auc(roccurve))

        # Store results
        pred_loglik_po[counter] <- pll_po
        pred_loglik_pa[counter] <- pll_pa
        roc_auc[counter] <- auc_val
        rm(pres_prob, pll_pa, pll_po, roccurve, auc_val)
      }

      if (is.null(m)) {
        # Store missing results in this case
        loglik[counter] <- NA
        aic[counter] <- NA
      } else {
        # Store results
        loglik[counter] <-  logLik.scampr(m)
        aic[counter] <- AIC.scampr(m)
        m <- NULL
      }
      # Store common results
      nbf[counter] <- nrow(bfs)
      print(paste0("Completed fits with ", nrow(bfs), " basis functions"))
      counter <- counter + 1
      counter.counter <- c(counter.counter, counter)
    }
  }

  # Adjust names of return object accoding to approx. type
  ret.frame <- as.data.frame(cbind(nodes.on.long.edge = counter.counter[1:(length(counter.counter) - 1)], bf = nbf, loglik = loglik, aic = aic, predicted_cll_po = pred_loglik_po, predicted_cll_pa = pred_loglik_pa, roc_auc = roc_auc))
  if (approx.with == "variational") {
    attr(ret.frame, "approx") <- "variational"
  } else {
    attr(ret.frame, "approx") <- "laplace"
  }
  return(ret.frame)
}
