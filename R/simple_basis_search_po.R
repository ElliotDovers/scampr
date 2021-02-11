simple_basis_search_po <- function(po.formula, po.data, po.fold.id, max.basis.functions, approx.with = c("variational", "laplace"), coord.names = c("x", "y"), quad.weights.name = "quad.size", radius.type = c("diag", "limiting"), bf.matrix.type = c("sparse", "dense"), domain.data) {

  # Use provided model data as domain.data if missing and check coords are present
  if (missing(domain.data)) {
    if (!all(coord.names %in% colnames(po.data))) {
      stop(paste0("coord.names, ", coord.names, ", not found in data set provided"))
    }
    domain.data <- po.data[ , coord.names]
  } else {
    if (!all(coord.names %in% colnames(domain.data))) {
      stop(paste0("coord.names, ", coord.names, ", not found in 'domain.data' frame provided"))
    }
  }
  # checks not covered by model fitting
  radius.type <- match.arg(radius.type)
  bf.matrix.type <- match.arg(bf.matrix.type)
  approx.with <- match.arg(approx.with)
  # Assign the maximum number of basis functions to searched, if missing
  if(missing(max.basis.functions)) {
    #  set at half the maximum number of presences in either data set
    max.basis.functions <- 0.5 * sum(po.data[, all.vars(po.formula[[2]])])
  }

  # If CV folds vector is missing just perform in-sample likelihood search
  if (missing(po.fold.id)) {

    # initialise the storage objects
    nbf <- NULL
    loglik <- NULL
    aic <- NULL
    pred_loglik_po <- NULL # these will remain NULL in this case
    counter.counter <- NULL
    counter <- 1
    counter.counter <- c(counter.counter, counter)

    # First iteration hard coded for IPP model, i.e. zero basis functions #

    m_ipp <- scampr::po(po.formula, po.data, model.type = "ipp", coord.names = coord.names, quad.weights.name = quad.weights.name, bf.matrix.type = bf.matrix.type)

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
      try(assign("m", po(po.formula, po.data, simple.basis = bfs, starting.pars = m_ipp, coord.names = coord.names, quad.weights.name = quad.weights.name, model.type = approx.with, bf.matrix.type = bf.matrix.type)))

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
    po.fold.id <- as.numeric(po.fold.id)

    # initialise the storage objects
    nbf <- NULL
    loglik <- NULL
    aic <- NULL
    pred_loglik_po <- NULL
    counter.counter <- NULL
    counter <- 1
    counter.counter <- c(counter.counter, counter)

    # First iteration hard coded for IPP model, i.e. zero basis functions #

    m_ipp <- scampr::po(po.formula, po.data, model.type = "ipp", coord.names = coord.names, quad.weights.name = quad.weights.name, bf.matrix.type = bf.matrix.type)

    # Initialise Objects #
    train_mods_ipp <- list() # training models
    test.inten <- NULL # predicted intensity rate
    po.resp <- all.vars(po.formula[[2]]) # PO response
    test.po.resp <- NULL # test PO response
    test.po.quad.sizes <- NULL
    # Loop through CV folds
    for (i in sort(unique(po.fold.id))) {
      train_mods_ipp[[i]] <- po(po.formula, po.data[po.fold.id != i, ], model.type = "ipp", coord.names = coord.names, quad.weights.name = quad.weights.name, bf.matrix.type = bf.matrix.type)
      test.inten <- c(test.inten, predict(train_mods_ipp[[i]], newdata = po.data[po.fold.id == i, ]))
      test.po.resp <- c(test.po.resp, po.data[po.fold.id == i, po.resp])
      test.po.quad.sizes <- c(test.po.quad.sizes, po.data[po.fold.id == i, quad.weights.name])
    }
    # predicted conditional log-likelihood on the PO data
    pll_po <- sum(test.inten[test.po.resp == 1]) - sum(test.po.quad.sizes[test.po.resp == 0] * exp(test.inten[test.po.resp == 0]))

    # Store Results
    nbf[counter] <- 0
    loglik[counter] <- logLik(m_ipp)
    aic[counter] <- AIC(m_ipp)
    pred_loglik_po[counter] <- pll_po
    print(paste0("Completed fits with 0 basis functions (IPP)"))
    counter <- 2
    counter.counter <- c(counter.counter, counter)
    bfs <- cbind(NA, NA)
    rm(pll_po)

    # Start looping until we hit the max number of basis functions
    while (nrow(bfs) <= max.basis.functions) {
      # Set the current simple basis function configuration
      bfs <- scampr::simple_basis(counter, data = domain.data, radius.type = radius.type)
      # try to fit the full model
      m <- NULL
      try(assign("m", po(po.formula, po.data, simple.basis = bfs, starting.pars = m_ipp, coord.names = coord.names, quad.weights.name = quad.weights.name, model.type = approx.with, bf.matrix.type = bf.matrix.type)))
      # Initialise Objects #
      train_mods <- list() # training models
      test.inten <- NULL # predicted intensity rate
      tmp.m <- NULL
      # Loop through CV folds
      for (i in sort(unique(po.fold.id))) {
        try(assign("tmp.m", po(po.formula, po.data[po.fold.id != i, ], simple.basis = bfs, starting.pars =  train_mods_ipp[[i]], coord.names = coord.names, quad.weights.name = quad.weights.name, model.type = approx.with, bf.matrix.type = bf.matrix.type)))
        # Only get predictions if the current model fits
        if (!is.null(tmp.m)) {
          train_mods[[i]] <- tmp.m
          test.inten <- c(test.inten, predict(train_mods[[i]], newdata = po.data[po.fold.id == i, ]))
          tmp.m <- NULL
        }
      }

      # Can only get the test results if ALL the training models fit
      if (length(train_mods) != length(unique(po.fold.id))) {
        # Store missing results in this case
        pred_loglik_po[counter] <- NA
      } else {
        # predicted conditional log-likelihood on the PO data
        pll_po <- sum(test.inten[test.po.resp == 1]) - sum(test.po.quad.sizes[test.po.resp == 0] * exp(test.inten[test.po.resp == 0]))

        # Store results
        pred_loglik_po[counter] <- pll_po
        rm(pll_po)
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
  ret.frame <- as.data.frame(cbind(counter = counter.counter[1:(length(counter.counter) - 1)], bf = nbf, loglik = loglik, aic = aic, predicted_cll_po = pred_loglik_po))
  if (approx.with == "variational") {
    attr(ret.frame, "approx") <- "variational"
  } else {
    attr(ret.frame, "approx") <- "laplace"
  }
  return(ret.frame)
}
