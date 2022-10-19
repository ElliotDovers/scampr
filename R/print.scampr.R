#' Print for objects of class 'scampr'
#'
#' @param x a scampr model object
#' @param ... NA
#'
#' @return printed argument - subset of summary.scampr
#' @export
#'
#' @importFrom stats pnorm coef
#'
#' @examples
#' # Get the gorilla nesting data
#' dat <- gorillas
#'
#' # Standardise the elevation covariate
#' dat$elev.std <- scale(dat$elevation)
#'
#' # Fit a scampr model to the point pattern
#' m <- scampr(pres ~ elev.std, data = dat, model.type = "ipp")
#'
#' m
print.scampr <- function(x, ...) {

  # Collect elements for reporting #

  # Model used
  if (!is.null(attr(x$formula, "bias"))) {
    fixed.formula <- object$formula
    bias.formula <- attr(object$formula, "bias")
    tmp.formula <- paste0(Reduce(paste, deparse(fixed.formula)), ", accounting for presence-only biasing with: ", Reduce(paste, deparse(bias.formula)))
  } else {
    tmp.formula <- as.character(x$formula)
  }

  # get an identifier for the model type
  mod.id <- NULL
  if (is.na(x$approx.type)) {
    mod.id <- "ipp"
  } else {
    mod.id <- x$approx.type
  }
  # Set the model description
  if (x$model.type == "IDM") {
    Model.Desc <- switch(mod.id,
                        ipp = "Integrated Data Model without spatial random effects",
                        variational = "Integrated Data Model with spatially correlated errors - Variational Approx.",
                        laplace = "Integrated Data Model with spatially correlated errors - Laplace approx.")
  } else if (x$model.type == "PO") {
    Model.Desc <- switch(mod.id,
                         ipp = "Inhomogeneous Poisson process",
                         variational = "Log-Gaussian Cox process - Variational Approx.",
                         laplace = "Log-Gaussian Cox process - Laplace Approx.")
  } else if (x$model.type == "PA") {
    Model.Desc <- switch(mod.id,
                         ipp = "Binary regression model w. complimentary log-log link function (without spatial random effects)",
                         variational = "Spatially correlated, binary regression model w. complimentary log-log link function - Variational approx.",
                         laplace = "Spatially correlated, binary regression model w. complimentary log-log link function - Laplace approx.")
  } else (
    stop("unknown 'model.type' in scampr model")
  )

  # Fixed effects of the model
  tmp.fixed_effects <- as.data.frame(x$fixed.effects)
  tmp.fixed_effects$`z value` <- tmp.fixed_effects$Estimate / tmp.fixed_effects$`Std. Error`
  tmp.fixed_effects$`Pr(>|z|)` <- stats::pnorm(abs(tmp.fixed_effects$`z value`), lower.tail = FALSE)

  # Print the report #

  cat(
    "Model Type: ", Model.Desc, "\n\nFormula: ", if (length(tmp.formula) > 1) { tmp.formula[c(2, 1, 3)] } else { tmp.formula }, "\n\napprox.", if (mod.id != "ipp") {"marginal"} else {""}, "logLik: ", logLik.scampr(x),
    "\n\nFixed Effect Coefficients:\n\n"
  )
  print(stats::coef(x)[1:nrow(x$fixed.effects)])
}
