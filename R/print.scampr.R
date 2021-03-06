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
  if (x$data.model.type == "popa") {
    po.formula <- x$formula
    pa.formula <- attr(x$formula, "pa")
    pa.resp <- all.vars(pa.formula[[2]])
    pa.pred <- all.vars(pa.formula[[3]])
    po.resp <- all.vars(po.formula[[2]])
    po.pred <- all.vars(po.formula[[3]])
    tmp.formula <- paste0(po.resp, " ~ ", paste(po.pred, collapse = " + "), " |&| ", pa.resp, " ~ ", paste(pa.pred, collapse = " + "))
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
  if (x$data.model.type == "popa") {
    Model.Desc <- switch(mod.id,
                        ipp = "Combined data inhomogeneous process model",
                        variational = "Combined data model w. spatially correlated errors - Variational approx.",
                        laplace = "Combined data model w. spatially correlated errors - Laplace approx.")
  } else if (x$data.model.type == "po") {
    Model.Desc <- switch(mod.id,
                         ipp = "Inhomogeneous Poisson process",
                         variational = "Log-Gaussian Cox process - Variational Approx.",
                         laplace = "Log-Gaussian Cox process - Laplace Approx.")
  } else if (x$data.model.type == "pa") {
    Model.Desc <- switch(mod.id,
                         ipp = "Binary regression model w. complimentary log-log link function",
                         variational = "Spatially correlated, binary regression model w. complimentary log-log link function - Variational approx.",
                         laplace = "Spatially correlated, binary regression model w. complimentary log-log link function - Laplace approx.")
  } else (
    stop("unknown 'data.model.type' in scampr model")
  )

  # Fixed effects of the model
  tmp.fixed_effects <- as.data.frame(x$fixed.effects)
  tmp.fixed_effects$`z value` <- tmp.fixed_effects$Estimate / tmp.fixed_effects$`Std. Error`
  tmp.fixed_effects$`Pr(>|z|)` <- stats::pnorm(abs(tmp.fixed_effects$`z value`), lower.tail = FALSE)

  # Print the report #

  cat(
    "Model Type: ", Model.Desc, "\n\nFormula: ", if (length(tmp.formula) > 1) { tmp.formula[c(2, 1, 3)] } else { tmp.formula }, "\n\napprox.", if (mod.id != "ipp") {"marginal"} else {""}, "logLik: ", logLik.scampr(x),
    "\n\nFixed Effects Ceofficients:\n\n"
  )
  print(stats::coef(x)[1:nrow(x$fixed.effects)])
}
