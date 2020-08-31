#' Summary for objects of class 'scampr'
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
summary.scampr <- function(object) {

  # Collect elements for reporting #

  # Model used
  tmp.formula <- as.character(object$formula)
  # Approximate marginal log-likelihood
  tmp.loglik <- -object$value
  # Akaike's Information Criteria
  tmp.aic <- -2*(tmp.loglik - length(unlist(object$starting.pars)))
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
    post.means <- tmp.random_effects[row.names(tmp.random_effects) == "random", 1]
    prior.variance <- tmp.random_effects[row.names(tmp.random_effects) == "PriorVar", 1]
    prior.variance <- as.data.frame(t(formatC(prior.variance, digits = 2)))
    colnames(prior.variance) <- paste0("res. ", 1:length(prior.variance))
    rownames(prior.variance) <- ""
    post.means.summary <- as.data.frame(t(formatC(c(quantile(post.means, probs = seq(0, 0.5, 0.25)), mean(post.means), quantile(post.means, probs = c(0.75, 1))), digits = 2)))
    colnames(post.means.summary) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
    rownames(post.means.summary) <- ""
  }

  # Print the report #

  cat(
    "Formula: ", tmp.formula[c(2, 1, 3)], "\n\nAIC: ",
    tmp.aic, "\ ", " approx. marginal logLik: ", tmp.loglik,
    "\n\nBasis functions per res. ", object$basis.per.res, "\n\nFixed Effects:\n\n"
  )
  printCoefmat(tmp.fixed_effects)
  cat(
    "---\n\nSpatial Random Effects:\n\nPosterior Means:\n"
  )
  print(post.means.summary)
  cat(
    "\nPrior Variance(s):\n"
  )
  print(prior.variance)
}
