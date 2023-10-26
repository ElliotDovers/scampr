#' Plot basis function search results - for objects of class 'basis.search'
#'
#' @description Plots the results of a simple basis function search created by \code{scampr::basis.search()}.
#'
#' @param x a scampr model fitted via \code{basis.search()} or object of class 'basis.search'
#' @param ... additional plotting arguments
#' @param plot.metric one or more numbers from 1 to 3 indicating the metric(s) to be plotted. 1: Fitted approximate log-likelihood. 2: Akaike's Information Criterion. 3: Bayesian Information Criterion
#'
#' @return NULL
#' @export
#'
#' @importFrom graphics plot.default points abline title
#'
#' @examples
#' # Get the flora data for one of the species
#' dat_po <- flora$po$sp1
#' dat_pa <- flora$pa
#'
#' # obtain a sample of 10,000 quadrature points for the point process model
#' set.seed(1)
#' quad.pts <- flora$quad[sample(1:nrow(flora$quad), 10000, replace = F), ]
#' set.seed(NULL)
#'
#' # Attach the quadrature points to the presence-only data
#' dat_po <- rbind.data.frame(dat_po, quad.pts)
#'
#' # Ensure the "response" variable in each data set shares the same name
#' dat_po$presence <- dat_po$pres
#' dat_pa$presence <- dat_pa$sp1
#'
#' # Fit an integrated data model without a latent effects to initialise the basis search
#' base.idm <- scampr(presence ~ MNT, dat_po, bias.formula = ~ D.Main,
#' pa.data = dat_pa, include.sre = F, model.type = "IDM", latent.po.biasing = F)
#'
#' # Perform the basis search and fit the optimised model
#' idm <- basis.search(base.idm, return.model = TRUE)
#'
#' # Plot the search results
#' plot_search(idm)
plot_search <- function(x, ..., plot.metric = 1) {

  # adjust for an optimised model
  if (is(x, "scampr")) {
    x <- attr(x, "search.res")
  }

  # adjust for dual optimised GRFs
  if (sum(x$opt) == 2) {
    par(mfrow = c(1,2), mar = c(4.1, 4.1, 4.1, 4.1))
    graphics::plot.default(x$k_full[x$data == "PA"], x[x$data == "PA", c("ll", "AIC", "BIC")][,plot.metric],
         xlab = "# Basis Functions", ylab = paste0("PA", c(" log-Likelihood", " AIC", " BIC")[plot.metric]),
         type = "l", main = "", col = "black", bty = "n",
         ylim = range(x[x$data == "PA",c("ll", "AIC", "BIC")][,plot.metric]))
    graphics::points(x$k_full[x$data == "PA"], x[x$data == "PA", c("ll", "AIC", "BIC")][,plot.metric])
    graphics::abline(v = x$k_full[x$data == "PA"][x$opt[x$data == "PA"]], col = "red", lty = "dashed")
    graphics::plot.default(x$k_full[x$data == "PO"], x[x$data == "PO", c("ll", "AIC", "BIC")][,plot.metric],
         xlab = "# Basis Functions", ylab = paste0("PO", c(" log-Likelihood", " AIC", " BIC")[plot.metric]),
         type = "l", main = "", col = "black", bty = "n",
         ylim = range(x[x$data == "PO",c("ll", "AIC", "BIC")][,plot.metric]))
    graphics::points(x$k_full[x$data == "PO"], x[x$data == "PO", c("ll", "AIC", "BIC")][,plot.metric])
    graphics::abline(v = x$k_full[x$data == "PO"][x$opt[x$data == "PO"]], col = "red", lty = "dashed")
    if (attr(x, "prune.stop")) {
      title(main = "Search stopped\ndue to PO data pruning")
    }
    par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
  } else {
    graphics::plot.default(x$k_full, x[ , c("ll", "AIC", "BIC")][,plot.metric],
         xlab = "# Basis Functions", ylab = paste0(x$data[1], c(" log-Likelihood", " AIC", " BIC")[plot.metric]),
         type = "l", main = "", col = "black", bty = "n",
         ylim = range(x[ ,c("ll", "AIC", "BIC")][,plot.metric]))
    graphics::points(x$k_full, x[, c("ll", "AIC", "BIC")][,plot.metric])
    graphics::abline(v = x$k_full[x$opt], col = "red", lty = "dashed")
    if (attr(x, "prune.stop")) {
      graphics::title(main = "Search stopped\ndue to PO data pruning")
    }
  }
}
