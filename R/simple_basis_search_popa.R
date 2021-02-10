simple_basis_search_popa <- function(po.formula, pa.formula, po.data, pa.data, po.fold.id, pa.fold.id, quad.weights.name = "quad.size", max.basis.functions) {

  # checks
  if (!identical(sort(unique(po.fold.id)), sort(unique(pa.fold.id)))) {
    stop("Both fold.id vectors must contain the same fold numbers")
  }
  if(missing(max.basis.functions)) {
    max.basis.functions <- 0.5 * max(c(sum(pa.data[ , all.vars(pa.formula[[2]])] > 0), sum(po.data[, all.vars(po.formula[[2]])])))
  }

  rad_type <- "diag"

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

  m_ipp <- scampr::popa(po.formula, pa.formula, po.data, pa.data, model.type = "ipp")

  # Initialise Objects #
  train_mods_ipp <- list() # training models
  test.abund <- NULL # predicted abundance rate
  test.inten <- NULL # predicted intensity rate
  pa.resp <- all.vars(pa.formula[[2]]) # PA response
  po.resp <- all.vars(po.formula[[2]]) # PO response
  test.pa.resp <- NULL # test PA response
  test.po.resp <- NULL # test PO response
  test.po.quad.sizes <- NULL
  for (i in sort(unique(po.fold.id))) {
    train_mods_ipp[[i]] <- popa(po.formula, pa.formula, po.data[po.fold.id != i, ], pa.data[pa.fold.id != i, ], model.type = "ipp")
    test.abund <- c(test.abund, predict(train_mods_ipp[[i]], newdata = pa.data[pa.fold.id == i, ], process = "abundance"))
    test.inten <- c(test.inten, predict(train_mods_ipp[[i]], newdata = po.data[po.fold.id == i, ]))
    test.pa.resp <- c(test.pa.resp, as.numeric(pa.data[pa.fold.id == i, pa.resp] > 0))
    test.po.resp <- c(test.po.resp, po.data[po.fold.id == i, po.resp])
    test.po.quad.sizes <- c(test.po.quad.sizes, po.data[po.fold.id == i, quad.weights.name])
  }
  # predicted presence probability
  pres_prob <- 1 - exp(-exp(test.abund))
  # predicted conditional log-likelihood on the PA data
  pll_pa <- sum((test.pa.resp * log(pres_prob)) - ((1 - test.pa.resp) * exp(test.abund)))
  # predicted conditional log-likelihood on the PO data
  pll_po <- sum(test.inten[test.po.resp == 1]) - sum(test.po.quad.sizes[test.po.resp == 0] * exp(test.inten[test.po.resp == 0]))
  # area under ROC curve
  roccurve <- pROC::roc(test.pa.resp, as.vector(pres_prob))
  auc_val <- as.numeric(pROC::auc(roccurve))

  # Store Results
  nbf[counter] <- 0
  loglik[counter] <- logLik(m_ipp)
  aic[counter] <- AIC(m_ipp)
  pred_loglik_po[counter] <- pll_po
  pred_loglik_pa[counter] <- pll_pa
  roc_auc[counter] <- auc_val
  counter <- 2
  counter.counter <- c(counter.counter, counter)
  bfs <- cbind(NA, NA)
  rm(pres_prob, pll_pa, pll_po, roccurve, auc_val)

  # Start looping until we hit the max number of basis functions
  while ((nrow(bfs) <= max.basis.functions)) {

    bfs <- scampr::simple_basis(counter, data = po.data, radius.type = rad_type)
    # try to fit the full model
    m <- NULL
    try(assign("m", popa(po.formula, pa.formula, po.data, pa.data, simple.basis = bfs, starting.pars = m_ipp)))
    # Initialise Objects #
    train_mods <- list() # training models
    test.abund <- NULL # predicted abundance rate
    test.inten <- NULL # predicted intensity rate
    tmp.m <- NULL
    for (i in sort(unique(po.fold.id))) {
      try(assign("tmp.m", popa(po.formula, pa.formula, po.data[po.fold.id != i, ], pa.data[pa.fold.id != i, ], simple.basis = bfs, starting.pars =  train_mods_ipp[[i]])))
      if (!is.null(tmp.m)) {
        train_mods[[i]] <- tmp.m
        test.abund <- c(test.abund, predict(train_mods[[i]], newdata = pa.data[pa.fold.id == i, ], process = "abundance"))
        test.inten <- c(test.inten, predict(train_mods[[i]], newdata = po.data[po.fold.id == i, ]))
        tmp.m <- NULL
      }
    }

    # Can only get the appropriate results if ALL the training models fit
    if (length(train_mods) != length(unique(po.fold.id))) {
      # Store missing results in this case
      pred_loglik_po[counter] <- NA
      pred_loglik_pa[counter] <- NA
      roc_auc[counter] <- NA
    } else {
      # predicted presence probability
      pres_prob <- 1 - exp(-exp(test.abund))
      # predicted conditional log-likelihood on the PA data
      pll_pa <- sum((test.pa.resp * log(pres_prob)) - ((1 - test.pa.resp) * exp(test.abund)))
      # predicted conditional log-likelihood on the PO data
      pll_po <- sum(test.inten[test.po.resp == 1]) - sum(test.po.quad.sizes[test.po.resp == 0] * exp(test.inten[test.po.resp == 0]))
      # area under ROC curve
      roccurve <- pROC::roc(test.pa.resp, as.vector(pres_prob))
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
      loglik[counter] <-  logLik(m)
      aic[counter] <- AIC(m)
      m <- NULL
    }
    # Store common results
    nbf[counter] <- nrow(bfs)
    print(counter)
    counter <- counter + 1
    counter.counter <- c(counter.counter, counter)
  }

  plot(nbf, loglik, xlab = "# Basis Functions", ylab = "log-Likelihood", type = "l")
  points(nbf, loglik)
  # abline(v = nbf[best.mod.id.loglik], col = "darkgreen")

  return(as.data.frame(cbind(counter = counter.counter[1:(length(counter.counter) - 1)], bf = nbf, loglik = loglik, aic = aic, pred_loglik_po = pred_loglik_po, pred_loglik_pa = pred_loglik_pa, roc_auc = roc_auc)))
}
