#' Print for objects of class 'scampr'
#'
#' @param object a scampr model
#'
#' @return
#' @export
#'
#' @examples
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
  Mod_Type <- switch(mod.id,
                     ipp = "IPP",
                     variational = "LGCP with Variational Approx.",
                     laplace = "LGCP with Laplace Approx.")

  # Fixed effects of the model
  tmp.fixed_effects <- as.data.frame(object$fixed.effects)
  tmp.fixed_effects$`z value` <- tmp.fixed_effects$Estimate / tmp.fixed_effects$`Std. Error`
  tmp.fixed_effects$`Pr(>|z|)` <- pnorm(abs(tmp.fixed_effects$`z value`), lower.tail = FALSE)

  # Print the report #

  cat(
    "Model Type: ", Mod_Type, "\n\nFormula: ", tmp.formula[c(2, 1, 3)], "\n\napprox.", if (mod.id != "ipp") {"marginal"} else {""}, "logLik: ", logLik(object),
    "\n\nFixed Effects Ceofficients:\n\n"
  )
  print(coef(object)[1:nrow(object$fixed.effects)])
}
