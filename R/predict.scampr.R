#' Predict function for objects of class 'scampr'
#'
#' @description Functions the same as predict.glm with additional functionality. Can select whether predictions come from the realised latent field (dens = "posterior") or unrealised (dens = "prior") - depending on whether a user wants to make data specific predictions or be more broad. Particular to combined data models (popa) the user can specify whether to return the rate at which presence records occur (process = intensity) or the abundance rate (process = abundance).
#'
#' @param object a scampr model object
#' @param newdata a data frame of point locations to predict over as well as predictors involved in the model
#' @param type a character string , one of 'link' or 'response', indicating the type of prediction to be returned. Eitehr log-intensity or intensity respectively.
#' @param dens a character string, one of 'posterior' or 'prior', indicating the probability density of the random effects to take the expectation from.
#' @param process a character string, one of 'intensity' or 'abundance', indictating the process to be estimated. Only available for combined data models.
#' @param ... NA
#'
#' @return a numeric vector of length newdata (or length of fitted data) containing the predictions.
#' @exportS3Method stats::predict scampr
#'
#' @importFrom methods as
#' @importFrom stats na.omit
#'
#' @examples
#' # Get the flora data for one of the species
#' dat_po <- flora$po$sp1
#' dat_pa <- flora$pa
#'
#' # Attach the quadrature to the PO data
#' dat_po <- rbind.data.frame(dat_po, flora$quad)
#'
#' # Set a train and test set
#' train_po <- dat_po[dat_po$x <= mean(c(dat_po$x, dat_pa$x)), ]
#' test_po <- dat_po[dat_po$x > mean(c(dat_po$x, dat_pa$x)), ]
#' train_pa <- dat_pa[dat_pa$x <= mean(c(dat_po$x, dat_pa$x)), ]
#' test_pa <- dat_pa[dat_pa$x > mean(c(dat_po$x, dat_pa$x)), ]
#'
#' # Set up a simple 2D grid of basis functions to fit a LGCP model to the data
#' bfs <- simple_basis(nodes.on.long.edge = 9, data = dat_po)
#'
#' # Fit an IPP model to the point pattern
#' m.ipp <- scampr(pres ~ MNT + D.Main, train_po, model.type = "ipp")
#'
#' # Fit a combined data model
#' m.comb <- scampr(pres ~ MNT + D.Main, dat_po, sp1 ~ MNT,
#' dat_pa, model.type = "ipp")
#'
#' # Fit presence/absence model
#' m.pa <- scampr(pa.formula = sp1 ~ MNT, pa.data = dat_pa, model.type = "ipp")
#'
#' \dontrun{
#' # Fit a LGCP model to the point pattern
#' m.lgcp_va <- scampr(pres ~ MNT + D.Main, dat_po, simple.basis = bfs)
#'
#' predict(m.ipp, test_po)
#' predict(m.comb, test_po, process = "intensity")
#' predict(m.comb, test_pa, process = "abundance")
#' predict(m.pa, test_pa)
#' predict(m.lgcp_va, test_po)
#' predict(m.lgcp_va, test_po, dens = "prior")
#' }
predict.scampr <- function(object, ..., newdata, type = c("link", "response"), dens = c("posterior", "prior"), process = c("intensity", "abundance")) {

  ## checks ##
  type <- match.arg(type)
  dens <- match.arg(dens)
  process <- match.arg(process)
  # Adjust the calculation based on required prediction (abundance or intensity)
  if (object$data.model.type == "popa") {
    # obtain the pa data
    data.pa <- attr(object$data, "pa")
    # obtain the formulas
    form.po <- object$formula
    form.pa <- attr(object$formula, "pa")
    # obtain bias covariates
    pa.pred <- all.vars(form.pa[[3]])
    po.pred <- all.vars(form.po[[3]])
    bias.preds <- po.pred[!po.pred %in% pa.pred]
    # adjust the model object to reflect abundance if required
    if (process == "abundance") {
      object$formula <- form.pa
      object$data <- data.pa
      # set the bias coefficients to zero for predicting mean abundance
      object$fixed.effects[is.na(object$fixed.effects)] <- 0 # in case of NA, NaN standard errors
      object$fixed.effects[rownames(object$fixed.effects) %in% c("(Bias Intercept)", bias.preds), ] <- NA
      object$fixed.effects <- stats::na.omit(object$fixed.effects)
    } else {
      object$formula <- form.po
    }
  }
  if (missing(newdata)) {
    newdata <- object$data
  }
  if (!all(all.vars(object$formula[[3]]) %in% colnames(newdata))) {
    stop("Not all predictors of the model are found in new data provided")
  }
  if (!all(object$coord.names %in% colnames(newdata))) {
    stop("Not all coordinates of the model are found in new data provided")
  }
  ###

  # Obtain the fixed effect matrix and basis function matrix (at new data or fitted)
  X <- get.design.matrix(object$formula, newdata) # this gives a single column for intercept (whether POPA model has multiple or not)
  if (!is.na(object$approx.type)) {
    Z <- get.bf.matrix(object, newdata[ , object$coord.names], bf.matrix.type = object$bf.matrix.type)
  }

  # Calculate fixed effect components of the linear predictor
  betas <- as.numeric(object$fixed.effects[ , 1L])
  intercept.term.id <- grepl("Intercept)", rownames(object$fixed.effects), fixed = T)
  # If there are two intercepts (due to combined data model) then we need to combine the intercepts now
  if (sum(intercept.term.id) > 1) {
    combined.intercepts <- sum(object$fixed.effects[intercept.term.id , 1L])
    betas[which(intercept.term.id)] <- NA
    betas[min(which(intercept.term.id))] <- combined.intercepts
    betas <- as.vector(stats::na.omit(betas))
  }
  Xb <- X %*% betas
  # Calculate random components of the linear predictor based on density required (and if the model is not an IPP)
  if (!is.na(object$approx.type) & dens == "posterior") {
    mu <- as.numeric(object$random.effects[grepl(" Mean ", row.names(object$random.effects), fixed = T), 1L])
    Zmu <- Z %*% mu
    # add in the second (PO biasing) latent field if present and the process required is "intensity"
    if (!is.null(object$bias.field) & process == "intensity") {
      mu2 <- as.numeric(object$bias.field[grepl(" Mean ", row.names(object$bias.field), fixed = T), 1L])
      Zmu2 <- Z %*% mu2
      Zmu <- Zmu + Zmu2
    }
  } else {
    Zmu <- 0
  }

  # Calculate the variance correction of the exponential expectation if required
  if (type == "response") {
    if (!is.na(object$approx.type) & dens == "prior") {
      sigma2 <- object$random.effects[grepl("Prior Var", row.names(object$random.effects), fixed = T), 1L]
      vars <- rep(sigma2, object$basis.per.res)
      Sigma <- matrix(0, sum(object$basis.per.res), sum(object$basis.per.res))
      diag(Sigma) <- vars
      if (object$bf.matrix.type == "sparse") {
        Sigma <- methods::as(Sigma, "sparseMatrix")
      }
      tmp.Z <- as.matrix(Z)
      ZSigZ <- NULL
      for (n in 1:nrow(tmp.Z)) {
        ZSigZ[n] <- as.numeric((Z[n, ] %*% Sigma) %*% Z[n, ])
      }
    } else if (!is.na(object$approx.type) & dens == "posterior") {
      if (object$approx.type == "variational") {
        vars <- object$random.effects[grepl("Posterior Var", row.names(object$random.effects), fixed = T), 1L]
        Sigma <- matrix(0, sum(object$basis.per.res), sum(object$basis.per.res))
        diag(Sigma) <- vars
        if (object$bf.matrix.type == "sparse") {
          Sigma <- methods::as(Sigma, "sparseMatrix")
        }
        tmp.Z <- as.matrix(Z)
        ZSigZ <- NULL
        for (n in 1:nrow(tmp.Z)) {
          ZSigZ[n] <- as.numeric((Z[n, ] %*% Sigma) %*% Z[n, ])
        }
      } else {
        fullCovMat <- vcov.scampr(object)
        Sigma <- fullCovMat[rownames(fullCovMat) == "random", colnames(fullCovMat) == "random"]
        if (object$bf.matrix.type == "sparse") {
          Sigma <- methods::as(Sigma, "sparseMatrix")
        }
        tmp.Z <- as.matrix(Z)
        ZSigZ <- NULL
        for (n in 1:nrow(tmp.Z)) {
          ZSigZ[n] <- as.numeric((Z[n, ] %*% Sigma) %*% Z[n, ])
        }
      }
    } else {
      ZSigZ <- 0
    }
  }

  pred <- switch(type,
                 link = as.numeric(Xb + Zmu),
                 response = as.numeric(exp(Xb + Zmu + 0.5 * ZSigZ))
                 )
  return(pred)
}
