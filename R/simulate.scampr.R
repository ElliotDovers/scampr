#' Simulate functiion for objects of class scampr
#'
#' @description Simulates point patterns from a specified scampr model.
#'
#' @param object a scampr model
#' @param nsim number of point patterns to simulate. Defaults to 1.
#' @param seed an integer for setting the random seed.
#' @param ... NA
#' @param domain.data optionally, a data frame of point locations that adequately cover the domain of interest, as well as, predictors involved in the model. If missing model quadrature is used. NOTE: if model quadrature is used and is an inadequate cover of the domain, simulation may be effected as interpolated variables will contain  many NA values.
#' @param rcoef.density a character string indicating the probability density of the random effects to draw from.
#' @param which.intensity a character string, one of 'expected' or 'sample', indicating whether to simulate from the expected intensity (with respect to 'type' and 'rcoef.density') or from a randomly sampled intensity (with respect to 'rcoef.density') respectively. This will work in conjuction with 'type' (changes expectation) and 'rcoef.density' (N(0,prior_var) or N(posterior_mean, posterior_var)). Toggling this depends on the use of simulation, e.g. 'expected' may be useful in validiating a model while 'sample' may be usefu for making predictions.
#' @param return.type a character string, one of 'data.frame' or 'ppp', indicating the object type of the simulation to be returned. Default is data frame for use in scampr models. \code{ppp} object is useful for interfacing with \code{spatstat::envelope} e.g.
#' @param nsurv an optional integer describing the number of survey sites to be included in the simulated presence/absenece data. Only applies to combined data models (popa), if null the survey sites of the original model are used.
#' @param log.expected a logical indicating whether to take the expectation of the log-intensity (ignores the variance correction). Only relevant to LGCP models for which \code{which.intensity} is 'expected'.
#'
#' @return Depends on return.type. Default is a point pattern object of class 'ppp' from spatstat. Otherwise can be set to return a data.frame describing the point pattern (with corresponding quadrature as an attribute). If nsim > 1 then returned as a list of either 'return.type'.
#' @exportS3Method stats::simulate scampr
#'
#' @importFrom spatstat rpoispp interp.im
#' @importFrom methods as
#' @importFrom MASS mvrnorm
#' @importFrom fields rdist
#'
#' @examples
#' dat <- scampr::gorillas
#' dat$elev <- scale(dat$elevation)
#' mod <- po(pres ~ elev, dat, model.type = "ipp")
#' \dontrun{pp <- simulate(mod)}
simulate.scampr <- function(object, nsim = 1, seed = NULL, ..., domain.data, rcoef.density = c("posterior", "prior"), which.intensity = c("expected", "sample"), return.type = c("data.frame", "ppp"), nsurv, log.expected = T) {

  ## checks ##
  if (log.expected) {
    type <- "link"
  } else {
    type <- "response"
  }
  rcoef.density <- match.arg(rcoef.density)
  which.intensity <- match.arg(which.intensity)
  return.type <- match.arg(return.type)
  nsim <- as.integer(nsim)

  if (object$data.model.type == "pa") {
    stop("No functionality to simulate from a presence/absence data model... yet")
  }
  if (nsim < 1) {
    stop("nsim must be a postive integer")
  }

  # If no region data supplied use model quadrature
  if (missing(domain.data)) {
    domain.data <- object$data[object$pt.quad.id == 0, ]
  }
  if (!all(all.vars(object$formula[[3]]) %in% colnames(domain.data))) {
    stop("Not all predictors of the model are found in new data provided")
  }
  if (!all(object$coord.names %in% colnames(domain.data))) {
    stop("Not all coordinates of the model are found in new data provided")
  }
  ###

  # Set a seed as required
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # For IPP models
  if (is.na(object$approx.type)) {
    intens <- predict.scampr(object = object, newdata = domain.data, type = "response")
    # Convert to spatstat image
    intens_im <- vec2im(intens, domain.data[ , object$coord.names[1]], domain.data[ , object$coord.names[2]])
    pp <- spatstat::rpoispp(lambda = intens_im, nsim = nsim)
  } else { # For LGCP models
    if (which.intensity == "expected") { # predict() handles the various posterior, prior, resposne and link cases
      # can use predict here. 'type' == "link" allows us to ignore the expectation correction
      intens <- predict.scampr(object = object, newdata = domain.data, type = type, dens = rcoef.density)
      # but we need to exponentiate intensity if using "link"
      if (type == "link") {
        intens <- exp(intens)
      }
      # Convert to spatstat image
      intens_im <- vec2im(intens, domain.data[ , object$coord.names[1]], domain.data[ , object$coord.names[2]])
      pp <- spatstat::rpoispp(lambda = intens_im, nsim = nsim)
    } else { # FOR SAMPLING # Need to construct the log-intensity from scratch using random coefficients

      # can cheat to get the fixed effects component by exploiting predict function
      Xb <- predict.scampr(object = object, newdata = domain.data, type = "link", dens = "prior")

      # get the coefficient's mean vector
      if (!is.na(object$approx.type) & rcoef.density == "posterior") {
        mu <- as.numeric(object$random.effects[grepl(" Mean ", row.names(object$random.effects), fixed = T), 1L])
      } else {
        mu <- rep(0, sum(object$basis.per.res))
      }
      # get the coefficient's variance-covariance matrix
      if (rcoef.density == "prior") {
        sigma2 <- object$random.effects[grepl("Prior Var", row.names(object$random.effects), fixed = T), 1L]
        vars <- rep(sigma2, object$basis.per.res)
        Sigma <- matrix(0, sum(object$basis.per.res), sum(object$basis.per.res))
        diag(Sigma) <- vars
        Sigma <- methods::as(Sigma, "sparseMatrix")
      } else if (rcoef.density == "posterior") {
        if (object$approx.type == "variational") {
          vars <- object$random.effects[grepl("Posterior Var", row.names(object$random.effects), fixed = T), 1L]
          Sigma <- matrix(0, sum(object$basis.per.res), sum(object$basis.per.res))
          diag(Sigma) <- vars
          Sigma <- methods::as(Sigma, "sparseMatrix")
        } else { # for Laplace model posterior
          fullCovMat <- vcov.scampr(object)
          Sigma <- fullCovMat[rownames(fullCovMat) == "random", colnames(fullCovMat) == "random"]
        }
      }
      # Calculate the basis function matrix
      Z <- get.bf.matrix(object, domain.data[ , object$coord.names])
      # Cannot use rpoispp() nsim argument as we need to re-sample the latent field each time
      if (nsim == 1) {
        # generate the coefficients
        u <- MASS::mvrnorm(mu = mu, Sigma = Sigma)
        Zu <- as.numeric(Z %*% u)
        intens <- exp(Xb + Zu)
        # Convert to spatstat image
        intens_im <- vec2im(intens, domain.data[ , object$coord.names[1]], domain.data[ , object$coord.names[2]])
        pp <- spatstat::rpoispp(lambda = intens_im, nsim = 1L)
      } else { # For more than one simulation
        pp <- list()
        u <- list()
        Zu <- list()
        intens <- list()
        for (sim in 1:nsim) {
          # generate the coefficients
          u[[sim]] <- MASS::mvrnorm(mu = mu, Sigma = Sigma)
          Zu[[sim]] <- as.numeric(Z %*% u[[sim]])
          intens[[sim]] <- exp(Xb + Zu[[sim]])
          # Convert to spatstat image
          intens_im <- vec2im(intens[[sim]], domain.data[ , object$coord.names[1]], domain.data[ , object$coord.names[2]])
          pp[[sim]] <- spatstat::rpoispp(lambda = intens_im, nsim = 1)
        }
        names(pp) <- paste("Simulation", 1:nsim)
        class(pp) <- c("ppplist", "solist", "anylist", "listof", "list")
      }
    }
  }

  # Additionally simulate presence/absence data if using a combined data model (popa)
  if (object$data.model.type == "popa") {
    # Collect the relevant info

    survey.form <- attr(object$formula, "pa")
    coord.names <- object$coord.names
    survey.resp <- all.vars(survey.form[[2]])
    survey.preds <- all.vars(survey.form[[3]])
    if (missing(nsurv)) {
      # get model's survey data
      survey.data <- attr(object$data, "pa")[ , c(coord.names, survey.resp, survey.preds)]
      # can cheat to get the fixed effects component over the survey sites by exploiting predict function
      Xb.surv <- predict.scampr(object = object, newdata = survey.data, type = "link", dens = "prior", process = "abund")
      if (nsim == 1) {
        # Set the latent field depending on if it was explicitly calculated previously
        if (which.intensity == "sample" & !is.na(object$approx.type)) {
          Zu_im <- vec2im(Zu, domain.data[ , object$coord.names[1]], domain.data[ , object$coord.names[2]])
        } else {
          Xb.pp <- predict.scampr(object = object, newdata = domain.data, type = "link", dens = "prior")
          Zu_im <- vec2im(log(intens) -  Xb.pp, domain.data[ , object$coord.names[1]], domain.data[ , object$coord.names[2]])
        }
        survey.data$latent.field <- spatstat::interp.im(Zu_im, survey.data[ , object$coord.names[1]], survey.data[ , object$coord.names[2]])
        survey.data$abund <- exp(Xb.surv + survey.data$latent.field)
        survey.data$pprob <- 1 -exp(-survey.data$abund)
        survey.data$`sim Y` <- stats::rbinom(length(survey.data$pprob), 1, survey.data$pprob)
      } else { # Multiple simulations
        for (i in 1:nsim) {
          # Set the latent field depending on if it was explicitly calculated previously
          if (which.intensity == "sample" & !is.na(object$approx.type)) {
            Zu_im <- vec2im(Zu[[i]], domain.data[ , object$coord.names[1]], domain.data[ , object$coord.names[2]])
          } else {
            Xb.pp <- predict.scampr(object = object, newdata = domain.data, type = "link", dens = "prior")
            Zu_im <- vec2im(log(intens) -  Xb.pp, domain.data[ , object$coord.names[1]], domain.data[ , object$coord.names[2]])
          }
          latent.field <- spatstat::interp.im(Zu_im, survey.data[ , object$coord.names[1]], survey.data[ , object$coord.names[2]])
          abund <- exp(Xb.surv + latent.field)
          pprob <- 1 -exp(-abund)
          pa <- stats::rbinom(length(pprob), 1, pprob)
          survey.data <- cbind(survey.data, latent.field, abund, pprob, pa)
        }
        colnames(survey.data) <- c(coord.names, survey.resp, survey.preds, paste(rep(c("latent.field", "abund", "pprob", "sim Y"), nsim), rep(1:nsim, each = 4)))
      }
    } else { # When we want new surveys
      survey.data <- NA
      stop("sampling new survey sites from 'nsurv' is yet to be built")
      # # randomly sample the survey sites
      # survey.data <- data.frame(rpoispp(nsurv / nrow(domain.data), win = pp$window))
      # preds.frame <- domain.data[ , survey.preds]
      # for (p in 1:ncol(preds.frame)) {
      #   # if a factor then assign nearest value
      #   if (is.factor(preds.frame[ , p])) {
      #     pres.dists <- fields::rdist(domain.data[ , coord.names], survey.data[ , 1:2]) # distances to nearest quad points
      #     nearest.quad <- apply(pres.dists, 2, which.min) # finds the closest
      #     survey.data <- cbind(survey.data, preds.frame[nearest.quad, p])
      #   } else { # otherwise interpolate
      #     survey.data <- cbind(survey.data, spatstat::interp.im(vec2im(preds.frame[ , p], domain.data[ , coord.names[1]], domain.data[ , coord.names[2]]), x = survey.data[ , 1], y = survey.data[ , 2]))
      #   }
      # }
      # colnames(survey.data) <- c(coord.names, survey.preds) # names consistent with model
    }
  }

  # reset the random seed
  if (!is.null(seed)) {
    set.seed(NULL)
  }

  # Change the object to be returned accordingly
  if (return.type == "data.frame") {
    # get the names of columns that don't need to be interpolated to presence points
    coord.names <- object$coord.names
    resp.name <- all.vars(object$formula[[2]])
    quad.wts.name <- object$quad.weights.name
    # get the predictors included in the model
    pred.names <- all.vars(object$formula[[3]])
    # set the quadrature data frame
    quad <- domain.data[ , c(coord.names, resp.name, quad.wts.name, pred.names)]
    # get the predictors frame for looping through
    preds.frame <- quad[ , !colnames(quad) %in% c(coord.names, resp.name, quad.wts.name)]
    # Act according to number of simulations
    if (nsim == 1) {
      # set the quadrature intensity
      quad$process.intensity <- intens
      # Set the latent field depending on if it was explicitly calculated previously
      if (which.intensity == "sample" & !is.na(object$approx.type)) {
        quad$latent.field <- Zu
      } else {
        quad$latent.field <- log(intens) -  predict.scampr(object = object, newdata = domain.data, type = "link", dens = "prior")
      }
      pres <- cbind(data.frame(pp), 0, 1) # make df
      for (p in 1:ncol(preds.frame)) {
        # if a factor then assign nearest value
        if (is.factor(preds.frame[ , p])) {
          pres.dists <- fields::rdist(quad[ , coord.names], pres[ , 1:2]) # distances to nearest quad points
          nearest.quad <- apply(pres.dists, 2, which.min) # finds the closest
          pres <- cbind(pres, preds.frame[nearest.quad, p])
        } else { # otherwise interpolate
          pres <- cbind(pres, spatstat::interp.im(vec2im(preds.frame[ , p], quad[ , coord.names[1]], quad[ , coord.names[2]]), x = pres[ , 1], y = pres[ , 2]))
        }
      }
      colnames(pres) <- c(coord.names, quad.wts.name, resp.name, pred.names) # names consistent with model
      pres$process.intensity <- spatstat::interp.im(vec2im(quad$process.intensity, quad[ , coord.names[1]], quad[ , coord.names[2]]), x = pres[ , 1], y = pres[ , 2])
      pres$latent.field <- spatstat::interp.im(vec2im(quad$latent.field, quad[ , coord.names[1]], quad[ , coord.names[2]]), x = pres[ , 1], y = pres[ , 2])
      attr(pres, "quad") <- quad
      pp <- pres
    } else { # adjust for multiple simulations
      tmp <- list()
      for (i in 1:nsim) {
        tmp.quad <- quad
        # Set the intensity and latent field depending on if it was explicitly calculated previously
        if (which.intensity == "sample" & !is.na(object$approx.type)) {
          tmp.quad$process.intensity <- intens[[i]]
          tmp.quad$latent.field <- Zu[[i]]
        } else {
          tmp.quad$process.intensity <- intens
          tmp.quad$latent.field <- log(intens) -  predict.scampr(object = object, newdata = domain.data, type = "link", dens = "prior")
        }
        # set up the point pattern frame
        pres <- cbind(data.frame(pp[[i]]), 0, 1) # make df
        for (p in 1:ncol(preds.frame)) {
          # if a factor then assign nearest value
          if (is.factor(preds.frame[ , p])) {
            pres.dists <- fields::rdist(quad[ , coord.names], pres[ , 1:2]) # distances to nearest quad points
            nearest.quad <- apply(pres.dists, 2, which.min) # finds the closest
            pres <- cbind(pres, preds.frame[nearest.quad, p])
          } else { # otherwise interpolate
            pres <- cbind(pres, spatstat::interp.im(vec2im(preds.frame[ , p], quad[ , coord.names[1]], quad[ , coord.names[2]]), x = pres[ , 1], y = pres[ , 2]))
          }
        }
        colnames(pres) <- c(coord.names, quad.wts.name, resp.name, pred.names) # names consistent with model
        pres$process.intensity <- spatstat::interp.im(vec2im(tmp.quad$process.intensity, tmp.quad[ , coord.names[1]], tmp.quad[ , coord.names[2]]), x = pres[ , 1], y = pres[ , 2])
        pres$latent.field <- spatstat::interp.im(vec2im(tmp.quad$latent.field, tmp.quad[ , coord.names[1]], tmp.quad[ , coord.names[2]]), x = pres[ , 1], y = pres[ , 2])
        attr(pres, "quad") <- tmp.quad
        tmp[[i]] <- pres
        rm(pres)
      }
      pp <- tmp
      names(pp) <- paste("Simulation", 1:nsim)
    }
  }
  # Include the survey data if necessary
  if (object$data.model.type == "popa") {
    attr(pp, "survey") <- survey.data
  }
  return(pp)
}
