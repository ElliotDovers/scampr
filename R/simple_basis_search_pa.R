simple_basis_search_pa <- function(pa.formula, pa.data, pa.fold.id, max.basis.functions, coord.names = c("x", "y"), radius.type = c("diag", "limiting"), bf.matrix.type = c("sparse", "dense"), domain.data, approx.with = c("laplace", "variational")) {

  # Use provided model data as domain.data if missing and check coords are present
  if (missing(domain.data)) {
    if (!all(coord.names %in% colnames(pa.data))) {
      stop(paste0("coord.names, ", coord.names, ", not found in data set provided"))
    }
    domain.data <- pa.data[ , coord.names]
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
    stop("variational approx. not available for PA models")
  }

  # Assign the maximum number of basis functions to searched, if missing
  if(missing(max.basis.functions)) {
    #  set at half the maximum number of presences in the data set
    max.basis.functions <- 0.5 * max(sum(pa.data[ , all.vars(pa.formula[[2]])] > 0))
  }

  # If CV folds vector is missing just perform in-sample likelihood search
  if (missing(pa.fold.id)) {

    # initialise the storage objects
    nbf <- NULL
    loglik <- NULL
    aic <- NULL
    pred_loglik_po <- NULL # legacy object to remain NULL
    pred_loglik_pa <- NULL # these will remain NULL in this case
    roc_auc <- NULL # these will remain NULL in this case
    counter.counter <- NULL
    counter <- 1
    counter.counter <- c(counter.counter, counter)

    # First iteration hard coded for IPP model, i.e. zero basis functions #

    m_ipp <- scampr::pa(pa.formula, pa.data, model.type = "ipp", coord.names = coord.names, bf.matrix.type = bf.matrix.type)

    # Store Results
    nbf[counter] <- 0
    loglik[counter] <- logLik(m_ipp)
    aic[counter] <- AIC(m_ipp)
    print(paste0("Completed fit with 0 basis functions (IPP)"))
    counter <- 2
    counter.counter <- c(counter.counter, counter)
    bfs <- cbind(NA, NA)

    # Start looping until we hit the max number of basis functions
    while (nrow(bfs) <= max.basis.functions) {
      # Set the current simple basis function configuration
      bfs <- scampr::simple_basis(counter, data = domain.data, radius.type = radius.type)
      # try to fit the full model
      m <- NULL
      try(assign("m", pa(pa.formula, pa.data, simple.basis = bfs, starting.pars = m_ipp, coord.names = coord.names, bf.matrix.type = bf.matrix.type)))

      if (is.null(m)) {
        # Store missing results in this case
        loglik[counter] <- NA
        aic[counter] <- NA
      } else {
        # Store results
        loglik[counter] <-  logLik(m)
        aic[counter] <- AIC(m)
        m <- NULL
      }
      # Store common results
      nbf[counter] <- nrow(bfs)
      print(paste0("Completed fit with ", nrow(bfs), " basis functions"))
      counter <- counter + 1
      counter.counter <- c(counter.counter, counter)
    }
  } else { # IF CV FOLDS ARE PRESENT

    # Convert CV folds to numerics in case factor is provided
    pa.fold.id <- as.numeric(pa.fold.id)

    # initialise the storage objects
    nbf <- NULL
    loglik <- NULL
    aic <- NULL
    pred_loglik_po <- NULL # legacy object to remain NULL
    pred_loglik_pa <- NULL
    roc_auc <- NULL
    counter.counter <- NULL
    counter <- 1
    counter.counter <- c(counter.counter, counter)

    # First iteration hard coded for IPP model, i.e. zero basis functions #

    m_ipp <- scampr::pa(pa.formula, pa.data, model.type = "ipp", coord.names = coord.names, bf.matrix.type = bf.matrix.type)

    # Initialise Objects #
    train_mods_ipp <- list() # training models
    test.abund <- NULL # predicted abundance rate
    pa.resp <- all.vars(pa.formula[[2]]) # PA response
    test.pa.resp <- NULL # test PA response
    # Loop through CV folds
    for (i in sort(unique(pa.fold.id))) {
      train_mods_ipp[[i]] <- pa(pa.formula, pa.data[pa.fold.id != i, ], model.type = "ipp", coord.names = coord.names, bf.matrix.type = bf.matrix.type)
      test.abund <- c(test.abund, predict(train_mods_ipp[[i]], newdata = pa.data[pa.fold.id == i, ], process = "abundance"))
      test.pa.resp <- c(test.pa.resp, as.numeric(pa.data[pa.fold.id == i, pa.resp] > 0))
    }
    # predicted presence probability
    pres_prob <- 1 - exp(-exp(test.abund))
    # predicted conditional log-likelihood on the PA data
    pll_pa <- sum((test.pa.resp * log(pres_prob)) - ((1 - test.pa.resp) * exp(test.abund)))
    # area under ROC curve
    roccurve <- pROC::roc(test.pa.resp, as.vector(pres_prob))
    auc_val <- as.numeric(pROC::auc(roccurve))

    # Store Results
    nbf[counter] <- 0
    loglik[counter] <- logLik(m_ipp)
    aic[counter] <- AIC(m_ipp)
    pred_loglik_pa[counter] <- pll_pa
    roc_auc[counter] <- auc_val
    print(paste0("Completed fits with 0 basis functions (IPP)"))
    counter <- 2
    counter.counter <- c(counter.counter, counter)
    bfs <- cbind(NA, NA)
    rm(pres_prob, pll_pa, roccurve, auc_val)

    # Start looping until we hit the max number of basis functions
    while (nrow(bfs) <= max.basis.functions) {
      # Set the current simple basis function configuration
      bfs <- scampr::simple_basis(counter, data = domain.data, radius.type = radius.type)
      # try to fit the full model
      m <- NULL
      try(assign("m", pa(pa.formula, pa.data, simple.basis = bfs, starting.pars = m_ipp, coord.names = coord.names, bf.matrix.type = bf.matrix.type)))
      # Initialise Objects #
      train_mods <- list() # training models
      test.abund <- NULL # predicted abundance rate
      tmp.m <- NULL
      # Loop through CV folds
      for (i in sort(unique(pa.fold.id))) {
        try(assign("tmp.m", pa(pa.formula, pa.data[pa.fold.id != i, ], simple.basis = bfs, starting.pars =  train_mods_ipp[[i]], coord.names = coord.names, bf.matrix.type = bf.matrix.type)))
        # Only get predictions if the current model fits
        if (!is.null(tmp.m)) {
          train_mods[[i]] <- tmp.m
          test.abund <- c(test.abund, predict(train_mods[[i]], newdata = pa.data[pa.fold.id == i, ], process = "abundance"))
          tmp.m <- NULL
        }
      }

      # Can only get the test results if ALL the training models fit
      if (length(train_mods) != length(unique(pa.fold.id))) {
        # Store missing results in this case
        pred_loglik_pa[counter] <- NA
        roc_auc[counter] <- NA
      } else {
        # predicted presence probability
        pres_prob <- 1 - exp(-exp(test.abund))
        # predicted conditional log-likelihood on the PA data
        pll_pa <- sum((test.pa.resp * log(pres_prob)) - ((1 - test.pa.resp) * exp(test.abund)))
        # area under ROC curve
        roccurve <- pROC::roc(test.pa.resp, as.vector(pres_prob))
        auc_val <- as.numeric(pROC::auc(roccurve))

        # Store results
        pred_loglik_pa[counter] <- pll_pa
        roc_auc[counter] <- auc_val
        rm(pres_prob, pll_pa, roccurve, auc_val)
      }

      if (is.null(m)) {
        # Store missing results in this case
        loglik[counter] <- NA
        aic[counter] <- NA
      } else {
        # Store results
        loglik[counter] <-  logLik(m)
        aic[counter] <- AIC(m)
        m <- NULL
      }
      # Store common results
      nbf[counter] <- nrow(bfs)
      print(paste0("Completed fits with ", nrow(bfs), " basis functions"))
      counter <- counter + 1
      counter.counter <- c(counter.counter, counter)
    }
  }

  plot(nbf, loglik, xlab = "# Basis Functions", ylab = "log-Likelihood (in-sample)", type = "l")
  points(nbf, loglik)

  # Adjust names of return object accoding to approx. type
  ret.frame <- as.data.frame(cbind(counter = counter.counter[1:(length(counter.counter) - 1)], bf = nbf, loglik = loglik, aic = aic, predicted_cll_pa = pred_loglik_pa, roc_auc = roc_auc))
  if (approx.with == "variational") {
    attr(ret.frame, "approx") <- "variational"
  } else {
    attr(ret.frame, "approx") <- "laplace"
  }
  return(ret.frame)
}
