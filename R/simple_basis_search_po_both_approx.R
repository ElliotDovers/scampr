#' Internal Basis Function Search Algorithm for PO models - stable version fitting both VA and Laplace models.
#'
#' @param po.formula object of class "formula" (or one that can be coerced to that class): a symbolic description of the Presence-only data model to be fitted. The 'response' must be a binary that indicates whether a datum is a presence or quadrature point. See GLM function for further formula details.
#' @param po.data a data frame containing predictors at both presence-records and quadrature as well as the po.formula 'response'.
#' @param po.fold.id Optional for cross-validation. An integer or factor vector the same length as the po.data that describes the CV fold that each location falls into.
#' @param max.basis.functions Optional. An integer describing a rough upper limit to the number of basis functions to search. Defaults to half the number of presences in the data sets.
#' @param coord.names a vector of character strings describing the column names of the coordinates in both data frames.
#' @param quad.weights.name a charater string of the column name of quadrature weights in the po.data.
#' @param radius.type a character string describing the type of radius length to use. One of 'diag' = diagonal dist. between nodes or 'limiting' = sqrt(Domain Area)/log(k).
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to use sparse or dense matrix computations for the basis functions created.
#' @param domain.data Optional. A data frame of columns 'coord.names' that contains at least the extremities of the domain of interest. Useful to ensure the same basis function configurations are created by 'simple_basis' if comparing to various searches.
#'
#' @return a data.frame with columns including- 'nodes.on.long.edge': number used in scampr::simple_basis to create basis configuration. 'bf': the number of basis functions. 'loglik': the fitting marginal log-likelihood. 'aic': the corresponding AIC. Optionally, 'predicted_cll_po': the conditional (on the latent field) Presence-only likelihood. Optional columns are the results from a cross-validation described by 'po.fold.id' and/or 'pa.fold.id'. (_va or _lp subscript for approx. type if both are calculated)
#' @noRd
#'
#' @examples
#' #' #' # Get the Eucalypt data
#' dat_po <- eucalypt[["po"]]
#' dat_pa <- eucalypt[["pa"]]
#'
#' \dontrun{
#' # Fit with a shared latent field
#' res <- simple_basis_search_po_both_approx(pres ~ TMP_MIN + D_MAIN_RDS,
#' po.data = dat_po)
#' }
simple_basis_search_po_both_approx <- function(po.formula, po.data, po.fold.id, max.basis.functions, coord.names = c("x", "y"), quad.weights.name = "quad.size", radius.type = c("diag", "limiting"), bf.matrix.type = c("sparse", "dense"), domain.data) {

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
  # Assign the maximum number of basis functions to searched, if missing
  if(missing(max.basis.functions)) {
    #  set at half the maximum number of presences in either data set
    max.basis.functions <- 0.5 * sum(po.data[, all.vars(po.formula[[2]])])
  }

  # If CV folds vector is missing just perform in-sample likelihood search
  if (missing(po.fold.id)) {

    # initialise the storage objects
    nbf <- NULL
    loglik_va <- NULL
    loglik_lp <- NULL
    aic_va <- NULL
    aic_lp <- NULL
    timing_va <- NULL
    timing_lp <- NULL
    pred_loglik_po_va <- NULL # these will remain NULL in this case
    pred_loglik_po_lp <- NULL # these will remain NULL in this case
    counter.counter <- NULL
    counter <- 1
    counter.counter <- c(counter.counter, counter)

    # First iteration hard coded for IPP model, i.e. zero basis functions #

    cpu.time <- system.time(assign("m_ipp", po(po.formula, po.data, model.type = "ipp", coord.names = coord.names, quad.weights.name = quad.weights.name, bf.matrix.type = bf.matrix.type)))

    # Store Results
    nbf[counter] <- 0
    loglik_va[counter] <- logLik.scampr(m_ipp)
    aic_va[counter] <- AIC.scampr(m_ipp)
    loglik_lp[counter] <- logLik.scampr(m_ipp)
    aic_lp[counter] <- AIC.scampr(m_ipp)
    timing_va[counter] <- cpu.time
    timing_lp[counter] <- cpu.time
    print(paste0("Completed fit with 0 basis functions (IPP)"))
    counter <- 2
    counter.counter <- c(counter.counter, counter)
    bfs <- cbind(NA, NA)

    # Start looping until we hit the max number of basis functions
    while (nrow(bfs) <= max.basis.functions) {
      # Set the current simple basis function configuration
      bfs <- simple_basis(counter, data = domain.data, radius.type = radius.type)
      # try to fit the full model
      m_va <- NULL
      m_lp <- NULL
      try(assign("cpu.time_va", system.time(assign("m_va", po(po.formula, po.data, simple.basis = bfs, starting.pars = m_ipp, coord.names = coord.names, quad.weights.name = quad.weights.name, model.type = "variational", bf.matrix.type = bf.matrix.type)))))
      if (!is.null(m_va)) {
        try(assign("cpu.time_lp", system.time(assign("m_lp", po(po.formula, po.data, simple.basis = bfs, starting.pars = m_va, coord.names = coord.names, quad.weights.name = quad.weights.name, model.type = "laplace", bf.matrix.type = bf.matrix.type)))))
      } else {
        try(assign("cpu.time_lp", system.time(assign("m_lp", po(po.formula, po.data, simple.basis = bfs, starting.pars = m_ipp, coord.names = coord.names, quad.weights.name = quad.weights.name, model.type = "laplace", bf.matrix.type = bf.matrix.type)))))
      }

      if (is.null(m_va)) {
        # Store missing results in this case
        loglik_va[counter] <- NA
        aic_va[counter] <- NA
        timing_va[counter] <- NA
      } else {
        # Store results
        loglik_va[counter] <-  logLik.scampr(m_va)
        aic_va[counter] <- AIC.scampr(m_va)
        timing_va[counter] <- cpu.time_va
        m_va <- NULL
      }
      if (is.null(m_lp)) {
        # Store missing results in this case
        loglik_lp[counter] <- NA
        aic_lp[counter] <- NA
        timing_lp[counter] <- NA
      } else {
        # Store results
        loglik_lp[counter] <-  logLik.scampr(m_lp)
        aic_lp[counter] <- AIC.scampr(m_lp)
        timing_lp[counter] <- cpu.time_lp
        m_lp <- NULL
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
    loglik_va <- NULL
    loglik_lp <- NULL
    aic_va <- NULL
    aic_lp <- NULL
    pred_loglik_po_va <- NULL
    pred_loglik_po_lp <- NULL
    counter.counter <- NULL
    counter <- 1
    counter.counter <- c(counter.counter, counter)

    # First iteration hard coded for IPP model, i.e. zero basis functions #

    m_ipp <- po(po.formula, po.data, model.type = "ipp", coord.names = coord.names, quad.weights.name = quad.weights.name, bf.matrix.type = bf.matrix.type)

    # Initialise Objects #
    train_mods_ipp <- list() # training models
    test.inten <- NULL # predicted intensity rate
    po.resp <- all.vars(po.formula[[2]]) # PO response
    test.po.resp <- NULL # test PO response
    test.po.quad.sizes <- NULL
    # Loop through CV folds
    for (i in sort(unique(po.fold.id))) {
      train_mods_ipp[[i]] <- po(po.formula, po.data[po.fold.id != i, ], model.type = "ipp", coord.names = coord.names, quad.weights.name = quad.weights.name, bf.matrix.type = bf.matrix.type)
      test.inten <- c(test.inten, predict.scampr(train_mods_ipp[[i]], newdata = po.data[po.fold.id == i, ]))
      test.po.resp <- c(test.po.resp, po.data[po.fold.id == i, po.resp])
      test.po.quad.sizes <- c(test.po.quad.sizes, po.data[po.fold.id == i, quad.weights.name])
    }
    # predicted conditional log-likelihood on the PO data
    pll_po <- sum(test.inten[test.po.resp == 1]) - sum(test.po.quad.sizes[test.po.resp == 0] * exp(test.inten[test.po.resp == 0]))

    # Store Results
    nbf[counter] <- 0
    loglik_va[counter] <- logLik.scampr(m_ipp)
    loglik_lp[counter] <- logLik.scampr(m_ipp)
    aic_va[counter] <- AIC.scampr(m_ipp)
    aic_lp[counter] <- AIC.scampr(m_ipp)
    pred_loglik_po_va[counter] <- pll_po
    pred_loglik_po_lp[counter] <- pll_po
    print(paste0("Completed fits with 0 basis functions (IPP)"))
    counter <- 2
    counter.counter <- c(counter.counter, counter)
    bfs <- cbind(NA, NA)
    rm(pll_po)

    # Start looping until we hit the max number of basis functions
    while (nrow(bfs) <= max.basis.functions) {
      # Set the current simple basis function configuration
      bfs <- simple_basis(counter, data = domain.data, radius.type = radius.type)
      # try to fit the full model
      m_va <- NULL
      m_lp <- NULL
      try(assign("m_va", po(po.formula, po.data, simple.basis = bfs, starting.pars = m_ipp, coord.names = coord.names, quad.weights.name = quad.weights.name, model.type = "variational", bf.matrix.type = bf.matrix.type)))
      if (!is.null(m_va)) {
        try(assign("m_lp", po(po.formula, po.data, simple.basis = bfs, starting.pars = m_va, coord.names = coord.names, quad.weights.name = quad.weights.name, model.type = "laplace", bf.matrix.type = bf.matrix.type)))
      } else {
        try(assign("m_lp", po(po.formula, po.data, simple.basis = bfs, starting.pars = m_ipp, coord.names = coord.names, quad.weights.name = quad.weights.name, model.type = "laplace", bf.matrix.type = bf.matrix.type)))
      }

      # Initialise Objects #
      train_mods_va <- list() # training models
      test.inten_va <- NULL # predicted intensity rate
      tmp.m_va <- NULL
      train_mods_lp <- list() # training models
      test.inten_lp <- NULL # predicted intensity rate
      tmp.m_lp <- NULL
      # Loop through CV folds
      for (i in sort(unique(po.fold.id))) {
        # Try VA model
        try(assign("tmp.m_va", po(po.formula, po.data[po.fold.id != i, ], simple.basis = bfs, starting.pars =  train_mods_ipp[[i]], coord.names = coord.names, quad.weights.name = quad.weights.name, model.type = "variational", bf.matrix.type = bf.matrix.type)))
        # Only get predictions if the current VA model fits
        if (!is.null(tmp.m_va)) {
          train_mods_va[[i]] <- tmp.m_va
          test.inten_va <- c(test.inten_va, predict.scampr(train_mods_va[[i]], newdata = po.data[po.fold.id == i, ]))
        }
        # If VA model exists try LP model with VA starts
        if (!is.null(tmp.m_va)) {
          try(assign("tmp.m_lp", po(po.formula, po.data[po.fold.id != i, ], simple.basis = bfs, starting.pars =  train_mods_va[[i]], coord.names = coord.names, quad.weights.name = quad.weights.name, model.type = "laplace", bf.matrix.type = bf.matrix.type)))
          # Only get predictions if the current LP model fits
          if (!is.null(tmp.m_lp)) {
            train_mods_lp[[i]] <- tmp.m_lp
            test.inten_lp <- c(test.inten_lp, predict.scampr(train_mods_lp[[i]], newdata = po.data[po.fold.id == i, ]))
          } else { #  LP model failed with VA starts, so try Lp with IPP starts
            try(assign("tmp.m_lp", po(po.formula, po.data[po.fold.id != i, ], simple.basis = bfs, starting.pars =  train_mods_ipp[[i]], coord.names = coord.names, quad.weights.name = quad.weights.name, model.type = "laplace", bf.matrix.type = bf.matrix.type)))
            # Only get predictions if the current LP model fits
            if (!is.null(tmp.m_lp)) {
              train_mods_lp[[i]] <- tmp.m_lp
              test.inten_lp <- c(test.inten_lp, predict.scampr(train_mods_lp[[i]], newdata = po.data[po.fold.id == i, ]))
            }
          }
        } else { # If VA failed, try (optimistically) LP model with IPP starts
          try(assign("tmp.m_lp", po(po.formula, po.data[po.fold.id != i, ], simple.basis = bfs, starting.pars =  train_mods_ipp[[i]], coord.names = coord.names, quad.weights.name = quad.weights.name, model.type = "laplace", bf.matrix.type = bf.matrix.type)))
          # Only get predictions if the current LP model fits
          if (!is.null(tmp.m_lp)) {
            train_mods_lp[[i]] <- tmp.m_lp
            test.inten_lp <- c(test.inten_lp, predict.scampr(train_mods_lp[[i]], newdata = po.data[po.fold.id == i, ]))
          }
        }
        # Reset the temp. models
        tmp.m_va <- NULL
        tmp.m_lp <- NULL
      }

      # Can only get the VA test results if ALL the VA training models fit
      if (length(train_mods_va) != length(unique(po.fold.id))) {
        # Store missing results in this case
        pred_loglik_po_va[counter] <- NA
      } else {
        # predicted conditional log-likelihood on the PO data
        pll_po_va <- sum(test.inten_va[test.po.resp == 1]) - sum(test.po.quad.sizes[test.po.resp == 0] * exp(test.inten_va[test.po.resp == 0]))

        # Store results
        pred_loglik_po_va[counter] <- pll_po_va
        rm(pll_po_va)
      }

      # Can only get the LP test results if ALL the LP training models fit
      if (length(train_mods_lp) != length(unique(po.fold.id))) {
        # Store missing results in this case
        pred_loglik_po_lp[counter] <- NA
      } else {
        # predicted conditional log-likelihood on the PO data
        pll_po_lp <- sum(test.inten_lp[test.po.resp == 1]) - sum(test.po.quad.sizes[test.po.resp == 0] * exp(test.inten_lp[test.po.resp == 0]))

        # Store results
        pred_loglik_po_lp[counter] <- pll_po_lp
        rm(pll_po_lp)
      }

      # Store VA in sample likelihood and aic if available
      if (is.null(m_va)) {
        # Store missing results in this case
        loglik_va[counter] <- NA
        aic_va[counter] <- NA
      } else {
        # Store results
        loglik_va[counter] <-  logLik.scampr(m_va)
        aic_va[counter] <- AIC.scampr(m_va)
        m_va <- NULL
      }
      # Store LP in sample likelihood and aic if available
      if (is.null(m_lp)) {
        # Store missing results in this case
        loglik_lp[counter] <- NA
        aic_lp[counter] <- NA
      } else {
        # Store results
        loglik_lp[counter] <-  logLik.scampr(m_lp)
        aic_lp[counter] <- AIC.scampr(m_lp)
        m_lp <- NULL
      }
      # Store common results
      nbf[counter] <- nrow(bfs)
      print(paste0("Completed fits with ", nrow(bfs), " basis functions"))
      counter <- counter + 1
      counter.counter <- c(counter.counter, counter)
    }
  }

  # Adjust names of return object accoding to approx. type
  ret.frame <- as.data.frame(cbind(nodes.on.long.edge = counter.counter[1:(length(counter.counter) - 1)], bf = nbf, loglik_va = loglik_va, loglik_lp = loglik_lp, aic_va = aic_va, aic_lp = aic_lp, timing_va = timing_va, timing_lp = timing_lp, predicted_cll_po_va = pred_loglik_po_va, predicted_cll_po_lp = pred_loglik_po_lp))
  attr(ret.frame, "approx") <- "both"
  return(ret.frame)
}
