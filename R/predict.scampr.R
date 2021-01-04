#' Predict function for objects of class 'scampr'
#'
#' @param object a scampr model object
#' @param newdata a data frame of point locations to predict over as well as predictors involved in the model
#' @param type a character string indicating the type of linear predictor to be returned. One of 'link' or 'response'
#' @param dens a character string indicating the pdf of the random effects to draw from
#' @param process a character string indictating the process to be estimated. one of 'intensity' or 'abundance'. Only available for combined data models.
#'
#' @return
#' @export
#'
#' @examples
predict.scampr <- function(object, newdata, type = c("link", "response"), dens = c("posterior", "prior"), process = c("intensity", "abundance")) {

  ## checks ##
  type <- match.arg(type)
  dens <- match.arg(dens)
  process <- match.arg(process)
  # Adjust the calculation based on required prediction (abundance or intensity)
  if (object$data.model.type == "popa") {
    data.po <- object$data
    data.pa <- attr(object$data, "pa")
    forms <- strsplit(object$formula, " |&| ", fixed = T)
    form.po <- as.formula(forms[[1]][1])
    form.pa <- as.formula(forms[[1]][2])
    if (process == "abundance") {
      object$formula <- form.pa
      object$data <- data.pa
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
  X <- scampr:::get.desgin.matrix(object$formula, newdata)
  if (!is.na(object$approx.type)) {
    Z <- scampr:::get.bf.matrix(object, newdata[ , object$coord.names])
  }

  # Calculate components of the linear predictor based on density required
  betas <- as.numeric(object$fixed.effects[ , 1L])
  intercept.term.id <- grepl("Intercept)", rownames(object$fixed.effects), fixed = T)
  if (sum(intercept.term.id) > 1) {
    combined.intercepts <- sum(object$fixed.effects[intercept.term.id , 1L])
    betas[which(intercept.term.id)] <- NA
    betas[min(which(intercept.term.id))] <- combined.intercepts
    betas <- as.vector(na.omit(betas))
  }
  Xb <- X %*% betas
  if (!is.na(object$approx.type) & dens == "posterior") {
    mu <- as.numeric(object$random.effects[grepl(" Mean ", row.names(object$random.effects), fixed = T), 1L])
    Zmu <- Z %*% mu
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
        Sigma <- as(Sigma, "sparseMatrix")
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
          Sigma <- as(Sigma, "sparseMatrix")
        }
        tmp.Z <- as.matrix(Z)
        ZSigZ <- NULL
        for (n in 1:nrow(tmp.Z)) {
          ZSigZ[n] <- as.numeric((Z[n, ] %*% Sigma) %*% Z[n, ])
        }
      } else {
        fullCovMat <- vcov(object)
        Sigma <- fullCovMat[rownames(fullCovMat) == "random", colnames(fullCovMat) == "random"]
        if (object$bf.matrix.type == "sparse") {
          Sigma <- as(Sigma, "sparseMatrix")
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
