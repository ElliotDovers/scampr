#' Predict function for objects of class 'scampr'
#'
#' @param object a scampr model object
#' @param newdata a data frame of point locations to predict over as well as predictors involved in the model
#' @param type character string indicating the type of linear predictor to be returned. One of 'link' or 'response'
#' @param dens character string indicating the pdf of the random effects to draw from
#'
#' @return
#' @export
#'
#' @examples
predict.scampr <- function(object, newdata, type = c("link", "response"), dens = c("posterior", "prior"), ...) {
  ## checks ##
  type <- match.arg(type)
  dens <- match.arg(dens)
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
  Xb <- X %*% betas
  if (!is.na(object$approx.type) & dens == "posterior") {
    mu <- as.numeric(object$random.effects[row.names(object$random.effects) == "random", 1L])
    Zmu <- Z %*% mu
  } else {
    Zmu <- 0
  }

  # Calculate the variance correction of the exponential expectation if required
  if (type == "response") {
    if (!is.na(object$approx.type) & dens == "prior") {
      sigma2 <- object$random.effects[row.names(object$random.effects) == "PriorVar", 1L]
      vars <- rep(sigma2, object$basis.per.res)
      Sigma <- matrix(0, sum(object$basis.per.res), sum(object$basis.per.res))
      diag(Sigma) <- vars
      Sigma <- as(Sigma, "sparseMatrix")
      tmp.Z <- as.matrix(Z)
      ZSigZ <- NULL
      for (n in 1:nrow(tmp.Z)) {
        ZSigZ[n] <- as.numeric((Z[n, ] %*% Sigma) %*% Z[n, ])
      }
    } else if (!is.na(object$approx.type) & dens == "posterior") {
      if (object$approx.type == "variational") {
        vars <- object$random.effects[row.names(object$random.effects) == "PosteriorVar", 1L]
        Sigma <- matrix(0, sum(object$basis.per.res), sum(object$basis.per.res))
        diag(Sigma) <- vars
        Sigma <- as(Sigma, "sparseMatrix")
        tmp.Z <- as.matrix(Z)
        ZSigZ <- NULL
        for (n in 1:nrow(tmp.Z)) {
          ZSigZ[n] <- as.numeric((Z[n, ] %*% Sigma) %*% Z[n, ])
        }
      } else {
        # Need to re-fit the model to obtain the covariance matrix in the laplace case
        des.mat <- scampr:::get.desgin.matrix(object$formula, object$data)
        pt.quad.id <- object$pt.quad.id
        fixed.names <- colnames(des.mat)
        bf.matrix <- scampr:::get.bf.matrix(object, object$data[ , object$coord.names])
        dat.list <- list(
          X_pres = des.mat[pt.quad.id == 1, ],
          Z_pres = as(bf.matrix[pt.quad.id == 1, ], "sparseMatrix"),
          X_quad = des.mat[pt.quad.id == 0, ],
          Z_quad = as(bf.matrix[pt.quad.id == 0, ], "sparseMatrix"),
          quad_size = object$data[ , object$quad.weights.name][pt.quad.id == 0],
          bf_per_res = object$basis.per.res,
          mod_type = 2
        )
        start.pars <- lapply(split(object$par, names(object$par)), unname)
        start.pars$random <- object$random.effects[row.names(object$random.effects) == "random", 1L]
        map.in <- start.pars
        for (p in 1:length(map.in)) {
          map.in[[p]] <- factor(rep(NA, length(start.pars[[p]])))
        }
        obj <- TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", silent = T)
        tmp <- TMB::sdreport(obj, getJointPrecision = T)
        fullPrecision <- tmp$jointPrecision
        fullCovMat <- solve(fullPrecision)
        Sigma <- fullCovMat[rownames(fullCovMat) == "random", colnames(fullCovMat) == "random"]
        Sigma <- as(Sigma, "sparseMatrix")
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
                 link = Xb + Zmu,
                 response = exp(Xb + Zmu + 0.5 * ZSigZ)
                 )
  return(pred)
}
