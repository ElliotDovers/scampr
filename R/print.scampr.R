#' Print for objects of class 'scampr'
#'
#' @param object a scampr model
#'
#' @return
#' @export
#'
print.scampr <- function(object) {

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
  # Set the model description
  if (object$data.model.type == "popa") {
    Model.Desc <- switch(mod.id,
                        ipp = "Combined data inhomogeneous process model",
                        variational = "Combined data model w. spatially correlated errors - Variational approx.",
                        laplace = "Combined data model w. spatially correlated errors - Laplace approx.")
  } else if (object$data.model.type == "po") {
    Model.Desc <- switch(mod.id,
                         ipp = "Inhomogeneous Poisson process",
                         variational = "LGCP with Variational Approx.",
                         laplace = "LGCP with Laplace Approx.")
  } else if (object$data.model.type == "pa") {
    Model.Desc <- switch(mod.id,
                         ipp = "Logistic Regression model w. complimentary log-log link function",
                         variational = "Spatially correlated logistic Regression model w. complimentary log-log link function - Variational approx.",
                         laplace = "Spatially correlated logistic Regression model w. complimentary log-log link function - Laplace approx.")
  } else (
    stop("unknown 'data.model.type' in scampr model")
  )

  # Fixed effects of the model
  tmp.fixed_effects <- as.data.frame(object$fixed.effects)
  tmp.fixed_effects$`z value` <- tmp.fixed_effects$Estimate / tmp.fixed_effects$`Std. Error`
  tmp.fixed_effects$`Pr(>|z|)` <- pnorm(abs(tmp.fixed_effects$`z value`), lower.tail = FALSE)

  # Print the report #

  cat(
    "Model Type: ", Model.Desc, "\n\nFormula: ", if (length(tmp.formula) > 1) { tmp.formula[c(2, 1, 3)] } else { tmp.formula }, "\n\napprox.", if (mod.id != "ipp") {"marginal"} else {""}, "logLik: ", logLik(object),
    "\n\nFixed Effects Ceofficients:\n\n"
  )
  print(coef(object)[1:nrow(object$fixed.effects)])
}
