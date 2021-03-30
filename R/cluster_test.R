#' Perform simulations from a scampr model to calculate global envelopes for the inhomogeneous K function
#'
#' @param model a scampr model object
#' @param conf.level a numeric between 0 and 1 indicating the confidence level of the global envelope test
#' @param nsims an integer describing the number of simulations to conduct
#' @param Kfun.correction a character string describing the correction type to be used. See \code{spatstat::Kinhom()} for details.
#' @param return.data a logical indicating whether to return the observed and simulated data
#'
#' @return Optionally, a list of the observed and simulated K functions - see \code{GET::create_curve_set()}. Results are plotted.
#' @export
#'
#' @importFrom spatstat owin
#' @importFrom mgcv in.out
#' @importFrom GET create_curve_set global_envelope_test
#' @importFrom graphics par polygon lines legend
#'
#' @examples
#' # Get the data
#' dat <- scampr::gorillas
#' dat$elev <- scale(dat$elevation)
#' mod <- po(pres ~ elev, dat, model.type = "ipp")
#' cluster_test(mod)
cluster_test <- function(model, conf.level = 0.95, nsims = 100, Kfun.correction = c("border", "bord.modif", "isotropic", "translate"), return.data = F) {

  # checks and data from model
  Kfun.correction = match.arg(Kfun.correction)
  if (class(model) != "scampr") {
    stop("provided model must be a scampr")
  }
  bnd <- bound(model)
  model.win <- spatstat::owin(poly = bnd)
  coord.names <- model$coord.names

  # Calculate the observed K function
  K_obs <- kfunc(model, correction = Kfun.correction, spatstat.win = model.win)

  # Calculate the simulated K functions
  K_sim <- NULL
  for (sim in 1:nsims) {
    tmp.dat <- simulate.scampr(model)
    in.bnd <- mgcv::in.out(as.matrix(bnd), as.matrix(tmp.dat[ , coord.names]))
    tmp.Kfn <- kfunc(model, point.pattern = tmp.dat[in.bnd, coord.names], correction = Kfun.correction, spatstat.win = model.win)
    K_sim <- cbind(K_sim, tmp.Kfn$Kfn)
  }

  # Use the GET functions to set up global envelopes
  datCurves <- GET::create_curve_set(list(obs = K_obs$Kfn, sim_m = K_sim))
  cr = GET::global_envelope_test(datCurves, type = "erl", alpha = 1 - conf.level)
  cr$r <- K_obs$dist

  graphics::par(mar = c(5.1, 4.1, 2.1, 2.1))
  plot.default(cr$r, cr$central, type = "n", ylim = range(c(cr$obs, cr$hi, cr$lo)),
       ylab = "inhomogeneous K", xlab = "distance")
  graphics::polygon(c(rev(cr$r), cr$r), c(rev(cr$hi), cr$lo), col = 'grey80', border = NA)
  graphics::lines(cr$r, cr$central, lty = "dashed")
  graphics::lines(cr$r, cr$obs, col = "red")
  graphics::legend("topleft", legend = c("Observed", "Sim. Central", paste0("Sim. ", conf.level*100, "% Bounds")), col = c("red", "black", "grey80"), lty = c("solid", "dashed", "solid"), cex = 1, lwd = rep(1.5, 3), bty = "n")

  if (return.data) {
    return(datCurves)
  }
}
