#' Simulate a Point Pattern from a scampr model
#'
#' @param model a scampr model
#' @param domain.data optionally, data frame of point locations that adequately cover the domain of interest, as well as, predictors involved in the model. If missing model quadrature is used.
#' @param type a character string indicating the type of prediction to be returned. One of 'link' or 'response' for log-intensity or intensity respectively.
#' @param dens a character string indicating the probability density of the random effects to draw from.
#' @param expected.or.sampled.intensity a character string indicating whether to simulate from the expected intensity (with respect to 'type' dens') or a randomly sampled intensity (with respect to 'dens'). One of 'expected' or 'sampled' respectively. This will work in conjuction with 'type' (changes expectation) and 'dens' (N(0,prior_var) or N(posterior_mean, posterior_var)). Toggling this depends on the use of simulation, e.g. 'sample' may be relevant to temporal models or for simulating future scenarios.
#' @param rseed an integer for setting the random seed.
#'
#' @return point pattern object of class 'ppp' from spatstat
#' @noRd
#'
#' @importFrom spatstat rpoispp
#' @importFrom methods as
#' @importFrom MASS mvrnorm
#' @examples
#' dat <- scampr::gorillas
#' dat$elev <- scale(dat$elevation)
#' mod <- po(pres ~ elev, dat, model.type = "ipp")
#' pp <- simulate_model_pp(mod)
simulate_model_pp <- function(model, domain.data, type = c("link", "response"), dens = c("posterior", "prior"), rseed = NA, expected.or.sampled.intensity = c("expected", "sampled")) {

  ## checks ##
  type <- match.arg(type)
  dens <- match.arg(dens)
  expected.or.sampled.intensity <- match.arg(expected.or.sampled.intensity)

  if (model$data.model.type == "pa") {
    stop("Cannot simulate a point pattern from a presence/absence data model")
  }

  # If no region data supplied use model quadrature
  if (missing(domain.data)) {
    domain.data <- model$data[model$pt.quad.id == 0, ]
  }
  if (!all(all.vars(model$formula[[3]]) %in% colnames(domain.data))) {
    stop("Not all predictors of the model are found in new data provided")
  }
  if (!all(model$coord.names %in% colnames(domain.data))) {
    stop("Not all coordinates of the model are found in new data provided")
  }
  ###

  # Set a seed as required
  if (!is.na(rseed)) {
    set.seed(rseed)
  }

  if (is.na(model$approx.type)) {
    intens <- predict.scampr(object = model, newdata = domain.data, type = "response")
    # Convert to spatstat image
    intens_im <- vec2im(intens, domain.data[ , model$coord.names[1]], domain.data[ , model$coord.names[2]])
    pp <- spatstat::rpoispp(lambda = intens_im)
  } else {
    if (expected.or.sampled.intensity == "expected") { # predict() handles the various posterior, prior, resposne and link cases
      # can use predict here. 'type' == "link" allows us to ignore the expectation correction
      intens <- predict.scampr(object = model, newdata = domain.data, type = type, dens = dens)
      # but we need to exponentiate intensity if using "link"
      if (type == "link") {
        intens <- exp(intens)
      }
      # Convert to spatstat image
      intens_im <- vec2im(intens, domain.data[ , model$coord.names[1]], domain.data[ , model$coord.names[2]])
      pp <- spatstat::rpoispp(lambda = intens_im)
    } else { # FOR SAMPLING # Need to construct the log-intensity from scratch using random coefficients

      # can cheat to get the fixed effects component by exploiting predict function
      Xb <- predict.scampr(object = model, newdata = domain.data, type = "link", dens = "prior")

      # get the coefficient's mean vector
      if (!is.na(model$approx.type) & dens == "posterior") {
        mu <- as.numeric(model$random.effects[grepl(" Mean ", row.names(model$random.effects), fixed = T), 1L])
      } else {
        mu <- rep(0, sum(model$basis.per.res))
      }
      # get the coefficient's variance-covariance matrix
      if (dens == "prior") {
        sigma2 <- model$random.effects[grepl("Prior Var", row.names(model$random.effects), fixed = T), 1L]
        vars <- rep(sigma2, model$basis.per.res)
        Sigma <- matrix(0, sum(model$basis.per.res), sum(model$basis.per.res))
        diag(Sigma) <- vars
        Sigma <- methods::as(Sigma, "sparseMatrix")
      } else if (dens == "posterior") {
        if (model$approx.type == "variational") {
          vars <- model$random.effects[grepl("Posterior Var", row.names(model$random.effects), fixed = T), 1L]
          Sigma <- matrix(0, sum(model$basis.per.res), sum(model$basis.per.res))
          diag(Sigma) <- vars
          Sigma <- methods::as(Sigma, "sparseMatrix")
        } else { # for Laplace model posterior
          fullCovMat <- vcov.scampr(model)
          Sigma <- fullCovMat[rownames(fullCovMat) == "random", colnames(fullCovMat) == "random"]
        }
      }

      # generate the coefficients
      u <- MASS::mvrnorm(mu = mu, Sigma = Sigma)
      # get the bf matrix
      Z <- get.bf.matrix(model, domain.data[ , model$coord.names])
      Zu <- as.numeric(Z %*% u)
      intens <- exp(Xb + Zu)
      # Convert to spatstat image
      intens_im <- vec2im(intens, domain.data[ , model$coord.names[1]], domain.data[ , model$coord.names[2]])
      pp <- spatstat::rpoispp(lambda = intens_im)
    }
  }
  if (!is.na(rseed)) {
    set.seed(NULL)
  }

  return(pp)
}
