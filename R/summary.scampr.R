#' Summary for objects of class 'scampr'
#'
#' @param object A scampr model object
#' @param ... NA
#'
#' @return a data.frame of the conditional model summary.
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
#' m <- scampr(pres ~ elev.std, data = dat, include.sre = F)
#'
#' summary(m)
summary.scampr <- function(object, ...) {

  # Collect elements for reporting #

  # Model used
  if (!is.null(attr(object$formula, "bias"))) {
    fixed.formula <- object$formula
    bias.formula <- attr(object$formula, "bias")
    tmp.formula <- paste0(Reduce(paste, deparse(fixed.formula)), ", accounting for presence-only biasing with: ", Reduce(paste, deparse(bias.formula)))
  } else {
    tmp.formula <- as.character(object$formula)
  }

  # get an identifier for the model type
  mod.id <- NULL
  if (object$approx.type == "not_sre" & !is.null(object$random.bias.effects)) {
    mod.id <- "laplace"
  } else {
    mod.id <- object$approx.type
  }

  # Set the model description
  if (object$model.type == "IDM") {
    Model.Desc <- switch(mod.id,
                         not_sre = "Integrated Data Model without spatial random effects",
                         variational = "Integrated Data Model with spatial random effects - Variational approx.",
                         laplace = "Integrated Data Model with spatial random effects - Laplace approx.")
  } else if (object$model.type == "PO") {
    Model.Desc <- switch(mod.id,
                         not_sre = "Inhomogeneous Poisson process",
                         variational = "Log-Gaussian Cox process - Variational Approx.",
                         laplace = "Log-Gaussian Cox process - Laplace Approx.")
  } else if (object$model.type == "PA") {
    Model.Desc <- switch(mod.id,
                         not_sre = "Binary regression model\nw. complimentary log-log link function (without spatial random effects)",
                         variational = "Binary regression model with spatial random effects\nw. complimentary log-log link function - Variational approx.",
                         laplace = "Binary regression model with spatial random effects\nw. complimentary log-log link function - Laplace approx.")
  } else (
    stop("unknown 'model.type' in scampr model")
  )


  # Fixed effects of the model
  tmp.fixed_effects <- as.data.frame(object$fixed.effects)
  tmp.fixed_effects$`z value` <- tmp.fixed_effects$Estimate / tmp.fixed_effects$`Std. Error`
  tmp.fixed_effects$`Pr(>|z|)` <- stats::pnorm(abs(tmp.fixed_effects$`z value`), lower.tail = FALSE)
  # save as object to return
  ret.obj <- tmp.fixed_effects
  # PO biasing fixed effects if present
  if (!is.null(object$fixed.bias.effects)) {
    tmp.fixed.bias.effects <- as.data.frame(object$fixed.bias.effects)
    tmp.fixed.bias.effects$`z value` <- tmp.fixed.bias.effects$Estimate / tmp.fixed.bias.effects$`Std. Error`
    tmp.fixed.bias.effects$`Pr(>|z|)` <- stats::pnorm(abs(tmp.fixed.bias.effects$`z value`), lower.tail = FALSE)
    # add to return object
    ret.obj <- rbind(ret.obj, tmp.fixed.bias.effects)
  }
  # Random effects of the model
  if (object$approx.type == "not_sre") {
    post.means.summary <- NA
    prior.variance <- NA
  } else {
    tmp.random_effects <- object$random.effects
    post.means <- tmp.random_effects[ , 1]
    prior.variance <- as.data.frame(t(formatC(object$prior.variances[!grepl("_bias", row.names(object$prior.variances), fixed = T), 1L], digits = 2)))
    colnames(prior.variance) <- paste0("res. ", 1:length(prior.variance))
    rownames(prior.variance) <- ""
    tmp.min <- as.data.frame(stats::aggregate(post.means, by = list(object$random.effects$res), FUN = function(x){stats::quantile(x, probs = 0)}))
    tmp.25 <- as.data.frame(stats::aggregate(post.means, by = list(object$random.effects$res), FUN = function(x){stats::quantile(x, probs = 0.25)}))
    tmp.50 <- as.data.frame(stats::aggregate(post.means, by = list(object$random.effects$res), FUN = function(x){stats::quantile(x, probs = 0.5)}))
    tmp.mean <- as.data.frame(stats::aggregate(post.means, by = list(object$random.effects$res), FUN = function(x){mean(x)}))
    tmp.75 <- as.data.frame(stats::aggregate(post.means, by = list(object$random.effects$res), FUN = function(x){stats::quantile(x, probs = 0.75)}))
    tmp.max <- as.data.frame(stats::aggregate(post.means, by = list(object$random.effects$res), FUN = function(x){stats::quantile(x, probs = 1)}))
    post.means.summary <- as.data.frame(cbind(tmp.min[ , -1], tmp.25[ , -1], tmp.50[ , -1], tmp.mean[ , -1], tmp.75[ , -1], tmp.max[ , -1]))
    colnames(post.means.summary) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  }
  # if the model is an IDM and contains a secondary latent field capturing bias in the presence-only data
  if (object$random.bias.type %in% c("field1", "field2")) {
    tmp.bias_field <- object$random.bias.effects
    post.means <- tmp.bias_field[, 1L]
    prior.variance_bias.field <- as.data.frame(t(formatC(object$prior.variances[grepl("_bias", row.names(object$prior.variances), fixed = T), 1L], digits = 2)))
    colnames(prior.variance_bias.field) <- paste0("res. ", 1:length(prior.variance_bias.field))
    rownames(prior.variance_bias.field) <- ""
    tmp.min <- as.data.frame(stats::aggregate(post.means, by = list(object$random.bias.effects$res), FUN = function(x){stats::quantile(x, probs = 0)}))
    tmp.25 <- as.data.frame(stats::aggregate(post.means, by = list(object$random.bias.effects$res), FUN = function(x){stats::quantile(x, probs = 0.25)}))
    tmp.50 <- as.data.frame(stats::aggregate(post.means, by = list(object$random.bias.effects$res), FUN = function(x){stats::quantile(x, probs = 0.5)}))
    tmp.mean <- as.data.frame(stats::aggregate(post.means, by = list(object$random.bias.effects$res), FUN = function(x){mean(x)}))
    tmp.75 <- as.data.frame(stats::aggregate(post.means, by = list(object$random.bias.effects$res), FUN = function(x){stats::quantile(x, probs = 0.75)}))
    tmp.max <- as.data.frame(stats::aggregate(post.means, by = list(object$random.bias.effects$res), FUN = function(x){stats::quantile(x, probs = 1)}))
    post.means.summary_bias.field <- as.data.frame(cbind(tmp.min[ , -1], tmp.25[ , -1], tmp.50[ , -1], tmp.mean[ , -1], tmp.75[ , -1], tmp.max[ , -1]))
    colnames(post.means.summary_bias.field) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  }

  # Print the report #

  cat(
    "Model Type: ", Model.Desc, "\n\nFormula: ", if (length(tmp.formula) > 1) { tmp.formula[c(2, 1, 3)] } else { tmp.formula }, "\n\nAIC: ",
    AIC.scampr(object), "\ ", " approx.", if (mod.id != "not_sre") {"marginal"} else {""}, "logLik: ", logLik.scampr(object), "\n\nFixed Effects:\n\n"
  )
  stats::printCoefmat(tmp.fixed_effects)
  if (!is.null(object$fixed.bias.effects)) {
    cat("\nFixed Biasing Effects:\n\n")
    stats::printCoefmat(tmp.fixed.bias.effects)
  }
  if (object$approx.type != "not_sre") {
    cat(
      "---\n\nSpatial Random Effects:\n", paste0("\nBasis functions per res. ", paste(object$basis.per.res, collapse = ", ")),
      "\n\nPosterior Means per Spatial Resolution(s):\n"
    )
    print(format(post.means.summary, digit = 2))
    cat(
      "\nPrior Variance(s):\n"
    )
    print(prior.variance)
  }
  if (object$random.bias.type %in% c("field1", "field2")) {
    cat(
      "---\n\nSpatial Random Effects used for Presence-only biasing:\n",
      paste0("\nBasis functions per res. ", paste(attr(object$basis.per.res, "bias"), collapse = ", ")),
      "\n\nPosterior Means per Spatial Resolution(s):\n"
    )
    print(format(post.means.summary_bias.field, digit = 2))
    cat(
      "\nPrior Variance(s):\n"
    )
    print(prior.variance_bias.field)
  }
  invisible(ret.obj)
}
