#' Plot the 2D fields from a scampr model ### WIP ###
#'
#' @description Plots a selected field from a fitted scampr model. Options include the fitted intensity or its components on the log scale.
#'
#' @param object a scampr model object
#' @param ... additional plotting arguments
#' @param domain.data Optional. A data frame of columns 'coord.names' that contains at least the extremities of the domain of interest. If missing the quadrature points from the model will be used - Note that if these do not form a cover of the domain the image will not work (likewise, will not work for scampr models that do not include quadrature).
#' @param field
#'
#' @return See spatstat::plot.im
#' @noRd
#'
#' @importFrom spatstat.geom plot.im
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
#' # Create the basis functions on the flora quadrature
#' bfs <-  simple_basis(9, flora$quad)
#'
#' # Combined Data Model
#' m.comb <- scampr(pres ~ MNT, dat_po, bias.formula = ~ D.Main, dat_pa, include.sre = F)
#'
#' # Plot the fitted log-intensity
#' plot_field(m.comb, domain.data = flora$quad)
plot_field <- function(object, domain.data, field = c("log Lambda", "Presence Probability", "Fixed Effects", "Fixed Bias Effects", "Latent Field", "Latent Bias Field"), ...) {

  field <- match.arg(field)

  # if domain data is missing, use the model's quadrature
  if (missing(domain.data)) {
    if (object$model.type == "PA") {
      stop("'domain.data' must be provided for a presence/absence model")
    } else {
      domain.data <- object$data[object$pt.quad.id == 0, ]
      warning("plotting will use the model's quadrature points, if this is not a complete cover of the domain please provide 'domain.data'")
    }
  }
  if (!all(all.vars(object$formula[[3]]) %in% colnames(domain.data))) {
    stop("Not all predictors of the model are found in new data provided")
  }
  if (!all(object$coord.names %in% colnames(domain.data))) {
    stop("Not all coordinates of the model are found in new data provided")
  }

  # Predict the fitted model across the domain data
  flds <- predict.scampr(object, newdata = domain.data, include.bias.accounting = TRUE)
  # select the requested field
  fld4plot <- switch(field,
                     `log Lambda` = flds,
                     `Presence Probability` = 1 - exp(-exp(attr(flds, "Xbeta") + attr(flds, "Zmu"))),
                     `Fixed Effects` = attr(flds, "Xbeta"),
                     `Fixed Bias Effects` = attr(flds, "Btau"),
                     `Latent Field` = attr(flds, "Zmu"),
                     `Latent Bias Field` = attr(flds, "Z2mu2")
  )
  # plot the desired field
  spatstat.geom::plot.im(vec2im(fld4plot, domain.data[,object$coord.names[1]], domain.data[,object$coord.names[2]]),
       main = field, box = FALSE, ...)
}
