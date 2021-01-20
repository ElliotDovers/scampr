#' Summary for objects of class 'scampr'
#'
#' @param object
#'
#' @return
#' @export
#'
summary.scampr <- function(object) {

  # Collect elements for reporting #
  # Model used
  tmp.formula <- as.character(object$formula)

  # get an identifier for the model type
  mod.id <- NULL
  if (is.na(object$approx.type)) {
    mod.id <- "ipp"
  } else {
    mod.id <- object$approx.type
  }
  Mod_Type <- switch(mod.id,
                     ipp = "IPP",
                     variational = "LGCP with Variational Approx.",
                     laplace = "LGCP with Laplace Approx.")
  # Set the model description
  if (object$data.model.type == "popa") {
    Model.Desc <- switch(mod.id,
                         ipp = "Combined data inhomogeneous process model",
                         variational = "Combined data model\n w. spatially correlated errors - Variational approx.",
                         laplace = "Combined data model\n w. spatially correlated errors - Laplace approx.")
  } else if (object$data.model.type == "po") {
    Model.Desc <- switch(mod.id,
                         ipp = "Inhomogeneous Poisson process",
                         variational = "LGCP with Variational Approx.",
                         laplace = "LGCP with Laplace Approx.")
  } else if (object$data.model.type == "pa") {
    Model.Desc <- switch(mod.id,
                         ipp = "Logistic Regression model\nw. complimentary log-log link function",
                         variational = "Spatially correlated logistic Regression model\nw. complimentary log-log link function - Variational approx.",
                         laplace = "Spatially correlated logistic Regression model\nw. complimentary log-log link function - Laplace approx.")
  } else (
    stop("unknown 'data.model.type' in scampr model")
  )


  # Fixed effects of the model
  tmp.fixed_effects <- as.data.frame(object$fixed.effects)
  tmp.fixed_effects$`z value` <- tmp.fixed_effects$Estimate / tmp.fixed_effects$`Std. Error`
  tmp.fixed_effects$`Pr(>|z|)` <- pnorm(abs(tmp.fixed_effects$`z value`), lower.tail = FALSE)
  # Random effects of the model
  if (length(object$random.effects) == 1) {
    if (is.na(object$random.effects)) {
      post.means.summary <- NA
      prior.variance <- NA
    }
  } else {
    tmp.random_effects <- object$random.effects
    post.means <- tmp.random_effects[grepl(" Mean ", row.names(tmp.random_effects), fixed = T), 1]
    prior.variance <- as.data.frame(t(formatC(tmp.random_effects[grepl("Prior Var", row.names(tmp.random_effects), fixed = T), 1], digits = 2)))
    colnames(prior.variance) <- paste0("res. ", 1:length(prior.variance))
    rownames(prior.variance) <- ""
    tmp.min <- as.data.frame(aggregate(post.means, by = list(object$basis.fn.info$res), FUN = function(x){quantile(x, probs = 0)}))
    tmp.25 <- as.data.frame(aggregate(post.means, by = list(object$basis.fn.info$res), FUN = function(x){quantile(x, probs = 0.25)}))
    tmp.50 <- as.data.frame(aggregate(post.means, by = list(object$basis.fn.info$res), FUN = function(x){quantile(x, probs = 0.5)}))
    tmp.mean <- as.data.frame(aggregate(post.means, by = list(object$basis.fn.info$res), FUN = function(x){mean(x)}))
    tmp.75 <- as.data.frame(aggregate(post.means, by = list(object$basis.fn.info$res), FUN = function(x){quantile(x, probs = 0.75)}))
    tmp.max <- as.data.frame(aggregate(post.means, by = list(object$basis.fn.info$res), FUN = function(x){quantile(x, probs = 1)}))
    post.means.summary <- as.data.frame(cbind(tmp.min[ , -1], tmp.25[ , -1], tmp.50[ , -1], tmp.mean[ , -1], tmp.75[ , -1], tmp.max[ , -1]))
    colnames(post.means.summary) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  }

  # Print the report #

  cat(
    "Model Type: ", Model.Desc, "\n\nFormula: ", if (length(tmp.formula) > 1) { tmp.formula[c(2, 1, 3)] } else { tmp.formula }, "\n\nAIC: ",
    AIC(object), "\ ", " approx.", if (mod.id != "ipp") {"marginal"} else {""}, "logLik: ", logLik(object),
    "\n\nBasis functions per res. ", object$basis.per.res, "\n\nFixed Effects:\n\n"
  )
  printCoefmat(tmp.fixed_effects)
  cat(
    "---\n\nSpatial Random Effects:\n\nPosterior Means per Spatial Resolution(s):\n"
  )
  print(format(post.means.summary, digit = 2))
  cat(
    "\nPrior Variance(s):\n"
  )
  print(prior.variance)
}
