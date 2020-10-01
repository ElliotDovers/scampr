#' Coefficients for objects of class 'scampr'
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
coefficients.scampr <- function(object) {
  tmp <- object$random.effects
  if (length(tmp) ==1L) {
    if (is.na(tmp)) {
      return(list(fixed = object$fixed.effects))
    } else {
      return(list(fixed = object$fixed.effects, random = tmp))
    }
  } else {
    tmp <- tmp[rownames(tmp) != "log_variance_component", ]
    return(list(fixed = object$fixed.effects, random = tmp))
  }
}
