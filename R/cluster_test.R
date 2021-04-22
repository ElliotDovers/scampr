#' Perform simulations from a scampr model to calculate global envelopes for the inhomogeneous K function
#'
#' @param model a scampr model object
#' @param conf.level a numeric between 0 and 1 indicating the confidence level of the global envelope test
#' @param nsims an integer describing the number of simulations to conduct
#' @param Kfun.correction a character string describing the correction type to be used. See \code{spatstat::Kinhom()} for details.
#' @param return.data a logical indicating whether to return the observed and simulated data
#' @param spatstat.win Optionally, a spatstat window object.
#' @param dists a vector of values for the distances at which the inhomogeneous K function should be evaluated. Not normally given by the user; \code{spatstat} provides a sensible default.
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
#' mod <- scampr(pres ~ elev, dat, model.type = "ipp")
#' cluster_test(mod)
cluster_test <- function(model, conf.level = 0.95, nsims = 100, Kfun.correction = c("border", "bord.modif", "isotropic", "translate"), return.data = F, spatstat.win, dists = NULL) {

  # checks and data from model
  Kfun.correction = match.arg(Kfun.correction)
  if (class(model) != "scampr") {
    stop("provided model must be a scampr")
  }
  # If the window is missing try and automatically determine one
  if (missing(spatstat.win)) {
    bnd <- bound(model)
    spatstat.win <- spatstat::owin(poly = bnd)
  }
  coord.names <- model$coord.names

  # Calculate the observed K function
  K_obs <- kfunc(model, correction = Kfun.correction, spatstat.win = spatstat.win, dists = dists)

  # Get the boundary into appropriate form for mgcv::in.out
  bnd <- as.matrix(cbind(spatstat.win$bdry[[1]][[1]], spatstat.win$bdry[[1]][[2]]))
  bnd <- rbind(bnd, bnd[1, ])

  # Calculate the simulated K functions
  K_sim <- NULL
  for (sim in 1:nsims) {
    tmp.dat <- simulate.scampr(model)
    in.bnd <- mgcv::in.out(as.matrix(bnd), as.matrix(tmp.dat[ , coord.names]))
    tmp.Kfn <- kfunc(model, point.pattern = tmp.dat[in.bnd, coord.names], correction = Kfun.correction, spatstat.win = spatstat.win, dists = K_obs$r)
    K_sim <- cbind(K_sim, tmp.Kfn$Kfn)
  }

  # Check and report simulation dramas with non-finite or NA values
  sim.ok <- apply(K_sim, 2, function(x){all(is.finite(x))})
  if (sum(sim.ok) == 0) {
    stop("Each simulated point pattern has created an NaN or NA K function value - try a different border correction?")
  } else if (sum(sim.ok) < length(sim.ok)) {
    warning(paste0("Only ", sum(sim.ok), " of the ", length(sim.ok), " simulations are being used due to NaN (or NA) K function values."))
    K_sim <- K_sim[ , sim.ok]
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

  attr(cr, "sim_curves") <- datCurves$funcs
  if (return.data) {
    return(cr)
  }
}
