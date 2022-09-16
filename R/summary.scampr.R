#' Summary for objects of class 'scampr'
#'
#' @param object A scampr model object
#' @param ... NA
#'
#' @return printed argument displaying summary information of the model fit.
#' @exportS3Method base::summary scampr
#'
#' @importFrom stats pnorm aggregate quantile printCoefmat
#'
#' @examples
#' # Get the gorilla nesting data
#' data(gorillas, package = "scampr")
#' dat <- gorillas
#'
#' # Standardise the elevation covariate
#' dat$elev.std <- scale(dat$elevation)
#'
#' # Fit a scampr model to the point pattern
#' m <- scampr(pres ~ elev.std, data = dat, model.type = "ipp")
#'
#' summary(m)
summary.scampr <- function(object, ...) {

  # Collect elements for reporting #

  # Model used
  if (object$data.model.type == "popa") {
    po.formula <- object$formula
    pa.formula <- attr(object$formula, "pa")
    pa.resp <- all.vars(pa.formula[[2]])
    pa.pred <- all.vars(pa.formula[[3]])
    po.resp <- all.vars(po.formula[[2]])
    po.pred <- all.vars(po.formula[[3]])
    tmp.formula <- paste0(po.resp, " ~ ", paste(po.pred, collapse = " + "), " |&| ", pa.resp, " ~ ", paste(pa.pred, collapse = " + "))
  } else {
    tmp.formula <- as.character(object$formula)
  }

  # get an identifier for the model type
  mod.id <- NULL
  if (is.na(object$approx.type)) {
    mod.id <- "ipp"
  } else {
    mod.id <- object$approx.type
  }

  # Set the model description
  if (object$data.model.type == "popa") {
    Model.Desc <- switch(mod.id,
                         ipp = "Combined data inhomogeneous process model",
                         variational = "Combined data model\n w. spatially correlated errors - Variational approx.",
                         laplace = "Combined data model\n w. spatially correlated errors - Laplace approx.")
  } else if (object$data.model.type == "po") {
    Model.Desc <- switch(mod.id,
                         ipp = "Inhomogeneous Poisson process",
                         variational = "Log-Gaussian Cox process - Variational Approx.",
                         laplace = "Log-Gaussian Cox process - Laplace Approx.")
  } else if (object$data.model.type == "pa") {
    Model.Desc <- switch(mod.id,
                         ipp = "Binary regression model\nw. complimentary log-log link function",
                         variational = "Spatially correlated binary regression model\nw. complimentary log-log link function - Variational approx.",
                         laplace = "Spatially correlated binary regression model\nw. complimentary log-log link function - Laplace approx.")
  } else (
    stop("unknown 'data.model.type' in scampr model")
  )


  # Fixed effects of the model
  tmp.fixed_effects <- as.data.frame(object$fixed.effects)
  tmp.fixed_effects$`z value` <- tmp.fixed_effects$Estimate / tmp.fixed_effects$`Std. Error`
  tmp.fixed_effects$`Pr(>|z|)` <- stats::pnorm(abs(tmp.fixed_effects$`z value`), lower.tail = FALSE)
  # Random effects of the model
  if (mod.id == "ipp") {
    post.means.summary <- NA
    prior.variance <- NA
  } else {
    tmp.random_effects <- object$random.effects
    post.means <- tmp.random_effects[grepl(" Mean ", row.names(tmp.random_effects), fixed = T), 1]
    prior.variance <- as.data.frame(t(formatC(tmp.random_effects[grepl("Prior Var", row.names(tmp.random_effects), fixed = T), 1], digits = 2)))
    colnames(prior.variance) <- paste0("res. ", 1:length(prior.variance))
    rownames(prior.variance) <- ""
    tmp.min <- as.data.frame(stats::aggregate(post.means, by = list(object$basis.fn.info$res), FUN = function(x){stats::quantile(x, probs = 0)}))
    tmp.25 <- as.data.frame(stats::aggregate(post.means, by = list(object$basis.fn.info$res), FUN = function(x){stats::quantile(x, probs = 0.25)}))
    tmp.50 <- as.data.frame(stats::aggregate(post.means, by = list(object$basis.fn.info$res), FUN = function(x){stats::quantile(x, probs = 0.5)}))
    tmp.mean <- as.data.frame(stats::aggregate(post.means, by = list(object$basis.fn.info$res), FUN = function(x){mean(x)}))
    tmp.75 <- as.data.frame(stats::aggregate(post.means, by = list(object$basis.fn.info$res), FUN = function(x){stats::quantile(x, probs = 0.75)}))
    tmp.max <- as.data.frame(stats::aggregate(post.means, by = list(object$basis.fn.info$res), FUN = function(x){stats::quantile(x, probs = 1)}))
    post.means.summary <- as.data.frame(cbind(tmp.min[ , -1], tmp.25[ , -1], tmp.50[ , -1], tmp.mean[ , -1], tmp.75[ , -1], tmp.max[ , -1]))
    colnames(post.means.summary) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
    # if the model contains a secondary latent field capturing bias in the presence-only data
    if (!is.null(object$bias.field)) {
      tmp.bias_field <- object$bias.field
      post.means <- tmp.bias_field[grepl(" Mean ", row.names(tmp.bias_field), fixed = T), 1]
      prior.variance_bias.field <- as.data.frame(t(formatC(tmp.bias_field[grepl("Prior Var", row.names(tmp.bias_field), fixed = T), 1], digits = 2)))
      colnames(prior.variance_bias.field) <- paste0("res. ", 1:length(prior.variance_bias.field))
      rownames(prior.variance_bias.field) <- ""
      tmp.min <- as.data.frame(stats::aggregate(post.means, by = list(object$basis.fn.info$res), FUN = function(x){stats::quantile(x, probs = 0)}))
      tmp.25 <- as.data.frame(stats::aggregate(post.means, by = list(object$basis.fn.info$res), FUN = function(x){stats::quantile(x, probs = 0.25)}))
      tmp.50 <- as.data.frame(stats::aggregate(post.means, by = list(object$basis.fn.info$res), FUN = function(x){stats::quantile(x, probs = 0.5)}))
      tmp.mean <- as.data.frame(stats::aggregate(post.means, by = list(object$basis.fn.info$res), FUN = function(x){mean(x)}))
      tmp.75 <- as.data.frame(stats::aggregate(post.means, by = list(object$basis.fn.info$res), FUN = function(x){stats::quantile(x, probs = 0.75)}))
      tmp.max <- as.data.frame(stats::aggregate(post.means, by = list(object$basis.fn.info$res), FUN = function(x){stats::quantile(x, probs = 1)}))
      post.means.summary_bias.field <- as.data.frame(cbind(tmp.min[ , -1], tmp.25[ , -1], tmp.50[ , -1], tmp.mean[ , -1], tmp.75[ , -1], tmp.max[ , -1]))
      colnames(post.means.summary_bias.field) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
    }
  }

  # Print the report #

  cat(
    "Model Type: ", Model.Desc, "\n\nFormula: ", if (length(tmp.formula) > 1) { tmp.formula[c(2, 1, 3)] } else { tmp.formula }, "\n\nAIC: ",
    AIC.scampr(object), "\ ", " approx.", if (mod.id != "ipp") {"marginal"} else {""}, "logLik: ", logLik.scampr(object),
    if (mod.id != "ipp") {paste0("\n\nBasis functions per res. ", paste(object$basis.per.res, collapse = ", "))} else {""}, "\n\nFixed Effects:\n\n"
  )
  stats::printCoefmat(tmp.fixed_effects)
  if (mod.id != "ipp") {
    cat(
      "---\n\nSpatial Random Effects:\n\nPosterior Means per Spatial Resolution(s):\n"
    )
    print(format(post.means.summary, digit = 2))
    cat(
      "\nPrior Variance(s):\n"
    )
    print(prior.variance)
    if (!is.null(object$bias.field)) {
      cat(
        "---\n\nSpatial Random Effects on the second latent field:\n\nPosterior Means per Spatial Resolution(s):\n"
      )
      print(format(post.means.summary_bias.field, digit = 2))
      cat(
        "\nPrior Variance(s):\n"
      )
      print(prior.variance_bias.field)
    }
  }
}
