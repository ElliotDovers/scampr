#' Predict function for objects of class 'scampr'
#'
#' @description Functions the same as predict.glm with additional functionality. Can select whether predictions come from the realised latent field (dens = "posterior") or unrealised (dens = "prior") - depending on whether a user wants to make data specific predictions or be more broad. Particular to combined data models the user can specify whether to return the rate at which presence records occur (use.formula = "presence-only") or the underlying abundance rate (use.formula = "presence-absence").
#'
#' @param object a scampr model object
#' @param newdata a data frame of point locations to predict over as well as predictors involved in the model
#' @param type a character string , one of 'link' or 'response', indicating the type of prediction to be returned. Either log-intensity or intensity respectively.
#' @param dens a character string, one of 'posterior' or 'prior', indicating the probability density of the random effects to take the expectation from.
#' @param include.bias.accounting a logical indicating if biasing effects (random or fixed) should be included in the predictions. Default is FALSE.
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
#'
#' # obtain a sample of 10,000 quadrature points for the point process model
#' set.seed(1)
#' quad.pts <- flora$quad[sample(1:nrow(flora$quad), 10000, replace = F), ]
#' set.seed(NULL)
#'
#' # Attach the quadrature points to the presence-only data
#' dat_po <- rbind.data.frame(dat_po, quad.pts)
#'
#' # Point Process Model
#' m <- scampr(pres ~ MNT + D.Main, dat_po, include.sre = F)
#'
#' # Set up some prediction points
#' newdat <- flora$quad[sample(1:nrow(flora$quad), 100, replace = F), ]
#'
#' # Make predictions
#' preds <- predict(m, newdat)
predict.scampr <- function(object, ..., newdata, type = c("link", "response"), dens = c("posterior", "prior"), include.bias.accounting = FALSE) {

  ## checks ####################################################################
  type <- match.arg(type)
  dens <- match.arg(dens)

  # if newdata is missing, use the model data
  if (missing(newdata)) {
    newdata <- object$data
  }
  if (!all(all.vars(object$formula[[3]]) %in% colnames(newdata))) {
    stop("Not all predictors of the model are found in new data provided")
  }
  if (!all(object$coord.names %in% colnames(newdata))) {
    stop("Not all coordinates of the model are found in new data provided")
  }
  # if (include.bias.accounting) {
  #   if (is.null(object$fixed.bias.effects) & is.null(object$random.bias.effects)) {
  #     warning("the model provided does not account for presence-only biasing, yet 'include.bias.accounting' = TRUE. This will be ignored")
  #   }
  # }
  ##############################################################################

  ## Fixed Effects

  # Obtain the fixed effect design matrix
  X <- get.design.matrix(object$formula, newdata)
  if (!all(colnames(X) == rownames(object$fixed.effects))) {
    stop("names of 'newdata' fixed effect design matrix do not match those of the model fixed effect coefficients")
  }
  # obtain the fixed effect coefficients
  beta <- object$fixed.effects[ , 1L]
  # calculate fixed effect components of the linear predictor
  Xbeta <- as.numeric(X %*% beta)

  ## Bias Effects

  if (include.bias.accounting) {
    # fixed biasing effects:
    if (!is.null(object$fixed.bias.effects)) {
      # Obtain the fixed biasing effect design matrix
      B <- get.design.matrix(attr(object$formula, "bias"), newdata)
      # adjust the intercept name if present
      if (any(grepl("(Intercept)", colnames(B), fixed = T))) {
        colnames(B)[grepl("(Intercept)", colnames(B), fixed = T)] <- "(Bias Intercept)"
      }
      if (!all(colnames(B) == rownames(object$fixed.bias.effects))) {
        stop("names of 'newdata' biasing fixed effect design matrix do not match those of the model biasing fixed effect coefficients")
      }
      # obtain the fixed biasing effect coefficients
      tau <- object$fixed.bias.effects[ , 1L]
      # calculate fixed biasing effect components of the linear predictor
      Btau <- as.numeric(B %*% tau)
    } else {
      Btau <- NULL
    }
    # random biasing effects:
    if (!is.null(object$po.biasing.basis.functions) & dens == "posterior") {
      # Obtain the random biasing effect basis function matrix
      Z2 <- get.bf.matrix(object$po.biasing.basis.functions, newdata[ , object$coord.names], bf.matrix.type = object$bf.matrix.type)
      # check the dimensions
      if (ncol(Z2) != nrow(object$random.bias.effects)) {
        stop("dimension of 'newdata' random biasing effect basis function matrix does not match the corresponding number of random biasing parameters in the model")
      }
      # obtain the random biasing effect coefficients
      mu2 <- object$random.bias.effects[ , 1L]
      # calculate random biasing effect components of the linear predictor
      Z2mu2 <- as.numeric(Z2 %*% mu2)
    } else {
      Z2mu2 <- NULL
    }
  } else {
    Btau <- NULL
    Z2mu2 <- NULL
  }

  ## Random Effects

  if (object$approx.type != "not_sre" & dens == "posterior") {
    # Obtain the random effect basis function matrix
    Z <- get.bf.matrix(object, newdata[ , object$coord.names], bf.matrix.type = object$bf.matrix.type)
    if (ncol(Z) != nrow(object$random.effects)) {
      stop("dimension of 'newdata' random effect basis function matrix does not match the corresponding number of random effect parameters in the model")
    }
    # obtain the random effect means
    mu <- object$random.effects[ , 1L]
    # calculate random biasing effect components of the linear predictor
    Zmu <- as.numeric(Z %*% mu)
  } else {
    # if there are no random effects OR we are interested in the prior density (where the means are zero)
    Zmu <- NULL
  }

  # ##############################################################################
  #  NEED TO FIX TO ADD VARIANCE CORRECTION FOR BOTH THE RANDOM AND RANDOM_BIASING COMPONENTS AS NEEDED
  # # Calculate the variance correction of the exponential expectation if required
  # if (type == "response") {
  #
  #   if (object$approx.type != "not_sre" & dens == "prior") {
  #     # For random effects (not PO biasing)
  #     sigma2 <- object$variances[!grepl("PriorVar_bias", row.names(object$variances), fixed = T), 1L]
  #     vars <- rep(sigma2, object$basis.per.res)
  #     Zsquared <- Z^2
  #     ZSigZ <- Zsquared %*% vars
  #
  #   } else if (object$approx.type != "not_sre" & dens == "posterior") {
  #     if (object$approx.type == "variational") {
  #       vars <- object$random.effects[grepl("Posterior Var", row.names(object$random.effects), fixed = T), 1L]
  #       Sigma <- matrix(0, sum(object$basis.per.res), sum(object$basis.per.res))
  #       diag(Sigma) <- vars
  #       if (object$bf.matrix.type == "sparse") {
  #         Sigma <- methods::as(Sigma, "sparseMatrix")
  #       }
  #       tmp.Z <- as.matrix(Z)
  #       ZSigZ <- NULL
  #       for (n in 1:nrow(tmp.Z)) {
  #         ZSigZ[n] <- as.numeric((Z[n, ] %*% Sigma) %*% Z[n, ])
  #       }
  #     } else {
  #       fullCovMat <- vcov.scampr(object)
  #       Sigma <- fullCovMat[rownames(fullCovMat) == "random", colnames(fullCovMat) == "random"]
  #       if (object$bf.matrix.type == "sparse") {
  #         Sigma <- methods::as(Sigma, "sparseMatrix")
  #       }
  #       tmp.Z <- as.matrix(Z)
  #       ZSigZ <- NULL
  #       for (n in 1:nrow(tmp.Z)) {
  #         ZSigZ[n] <- as.numeric((Z[n, ] %*% Sigma) %*% Z[n, ])
  #       }
  #     }
  #   } else {
  #     ZSigZ <- 0
  #   }
  # }
  # FOR NOW SETTING COMPONENTS TO ZERO
  ZSigZ <- NULL
  Z2Sig2Z2 <- NULL

  # make the predictions based on link type
  pred <- switch(type,
                 link = rowSums(cbind(Xbeta, Btau, Zmu, Z2mu2)),
                 response = exp(rowSums(cbind(Xbeta, Btau, Zmu, Z2mu2, (0.5 * ZSigZ), (0.5 * Z2Sig2Z2))))
  )
  attr(pred, "Xbeta") <- Xbeta
  attr(pred, "Btau") <- Btau
  attr(pred, "Zmu") <- Zmu
  attr(pred, "Z2mu2") <- Z2mu2
  attr(pred, "ZSigZ") <- ZSigZ
  attr(pred, "Z2Sig2Z2") <- Z2Sig2Z2
  return(pred)
}
