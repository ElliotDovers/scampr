#' Internal scampr function that creates a fixed effect offset vector from a formula and data
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of a model. See GLM function for further formula details.
#' @param data a data frame containing predictors and response of the formula parameter.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action a function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset. The ‘factory-fresh’ default is na.omit. Another possible value is NULL, no action. Value na.exclude can be useful.
#'
#' @return a data.frame
#' @export
#'
#' @importFrom stats model.frame model.matrix update.formula
#'
#' @examples
#' des.mat <- scampr:::get.offset(Petal.Length ~ Petal.Width + Species + offset(Petal.Length), iris)
#' head(des.mat)
get.offset <- function(formula, data, subset, na.action) {
  mf <- match.call(expand.dots = FALSE)
  m <- match(x = c("formula", "data", "subset", "na.action"),
             table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  # Update formula so that this will work if formula response is missing (for predictions at new data)
  mf$formula <- stats::update.formula(formula, NULL ~ .)
  mf <- eval(expr = mf, envir = parent.frame())
  offset <- as.vector(model.offset(mf))
  if (is.null(offset)) {
    offset <- rep(NA, nrow(data))
    attr(offset, "check") <- FALSE
  } else {
    attr(offset, "check") <- TRUE
  }
  return(offset)
}
