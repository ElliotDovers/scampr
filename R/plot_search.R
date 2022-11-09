#' Plot basis function search results - for objects of class 'basis.search' ### WIP ###
#'
#' @description Plots the results of a simple basis function search created by \code{scampr::basis.search()}.
#'
#' @param x an object of class 'basis.search'
#' @param ... additional plotting arguments
#'
#' @return
#' @export
#'
#' @importFrom graphics plot.default points abline
#'
#' @examples
#' # Get the flora data for one of the species
#' dat_po <- flora$po$sp1
#' dat_pa <- flora$pa
#'
#' # Create the basis functions on the flora quadrature
#' bfs <-  simple_basis(9, flora$quad)
#'
#' # Plot the bi-square basis functions
#' plot(bfs)
#'
#' # Fit a combined data model
#' m.comb <- scampr(pres ~ MNT + D.Main, dat_po, sp1 ~ MNT,
#' dat_pa, model.type = "ipp")
#'
#' # Fit presence/absence model
#' m.pa <- scampr(pa.formula = sp1 ~ MNT, pa.data = dat_pa, model.type = "ipp")
#'
#' plot(m.ipp)
#' plot(m.comb)
#' plot(m.pa)
plot_search <- function(x, ...) {
  if (!is(x, "basis.search")) {
    stop("Please provide an object of class 'basis.search'")
  }
  metrics <- colnames(x)[colnames(x) %in% c("fitted.ll", "predicted.ll_po", "predicted.ll_pa", "AUC")]
  nice.names <- c("Fitted approx. log-Likelihood", "Predicted Conditional logLik. (PPM)", "Predicted Conditional logLik. (BINOMIAL)", "AUC")
  nplots <- length(metrics)
  image_flag <- all(c("k", "k_bias") %in% colnames(x))
  # sort out the plot dimensions
  if (nplots == 4) {
    layout(mat = matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE),
           heights = c(0.5,0.5), widths = c(0.5,0.5))
  } else if (nplots == 3) {
    layout(mat = matrix(1:3, nrow = 1, ncol = 3, byrow = TRUE),
                  heights = 1, widths = rep(1/3,3))
  } else if (nplots == 2) {
    layout(mat = matrix(1:2, nrow = 1, ncol = 2, byrow = TRUE),
           heights = 1, widths = c(0.5,0.5))
  } else {
    layout(mat = matrix(1, nrow = 1, ncol = 1, byrow = TRUE),
           heights = 1, widths = 1)
  }
  for (p in 1:length(metrics)) {
    if (image_flag) {
      fields::image.plot(x = unique(x$k), y = unique(x$k_bias), z = vec2mat(x[,metrics[p]], x$k, x$k_bias),
                         xlab = "# basis functions in shared field",
                         ylab = "# basis functions in PO bias field",
                         main = nice.names[p])
    } else {
      graphics::plot.default(x$k, x[,metrics[p]], type = "l",
                             xlab = "# basis functions in shared field", ylab = "",
                             main = nice.names[p])
      graphics::points(x$k, x[,metrics[p]])
      graphics::abline(v = x$k[which.max(x[,metrics[p]])], col = "red", lty = "dashed")
    }
  }
  layout(mat = matrix(1, nrow = 1, ncol = 1, byrow = TRUE),
         heights = 1, widths = 1)
}
