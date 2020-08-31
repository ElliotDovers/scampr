#' Print for objects of class 'scampr'
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
print.scampr <- function(object) {

  # Collect elements for reporting #

  # Model used
  tmp.formula <- as.character(object$formula)
  # Approximate marginal log-likelihood
  tmp.loglik <- -object$value

  # Fixed effects of the model
  tmp.fixed_effects <- as.data.frame(object$fixed.effects)
  tmp.fixed_effects$`z value` <- tmp.fixed_effects$Estimate / tmp.fixed_effects$`Std. Error`
  tmp.fixed_effects$`Pr(>|z|)` <- pnorm(abs(tmp.fixed_effects$`z value`), lower.tail = FALSE)

  # Print the report #

  cat(
    "Formula: ", tmp.formula[c(2, 1, 3)], "\n\napprox. marginal logLik: ", tmp.loglik,
    "\n\nFixed Effects:\n\n"
  )
  printCoefmat(tmp.fixed_effects)
}
