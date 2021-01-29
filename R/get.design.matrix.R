#' Internal scampr function that creates a fixed effect design matrix from a formula and data
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of a model. See GLM function for further formula details.
#' @param data a data frame containing predictors and response of the formula parameter.
#'
#' @return
#' @noRd
#'
#' @importFrom stats model.frame model.matrix
#'
#' @examples
#' des.mat <- scampr:::get.desgin.matrix(Petal.Length ~ Petal.Width + Species, iris)
#' head(des.mat)
get.desgin.matrix <- function(formula, data) {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(x = c("formula", "data"),
             table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(expr = mf, envir = parent.frame())

  mt <- attr(x = mf, which = "terms")
  des.mat <- stats::model.matrix(object = mt, data = mf)
  return(des.mat)
}
