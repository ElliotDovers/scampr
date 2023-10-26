#' Inhomogeneous K function for scampr models (is a wrapper of \code{spatstat}'s version). Or calculate the K function from point pattern (and intensity) provided.
#'
#' @param model a scampr model object
#' @param correction a character string describing the correction type to be used. See \code{spatstat.explore::Kinhom()} for details.
#' @param point.pattern Optionally, a data frame, the first two columns of which describe the horizontal and vertical coordinates of point locations respectively.
#' @param intensity.at.pp Optionally, a vector of length \code{nrow(point.pattern)} that describes the intensity at each point.
#' @param spatstat.win Optionally, a spatstat window object.
#' @param intensity.at Optionally, a character string, one of either 'im' or 'pts'. Describes whether the intensity used in \code{spatstat.explore::Kinhom()} is an image object or just provided at the presence points.
#' @param dists a vector of values for the distances at which the inhomogeneous K function should be evaluated. Not normally given by the user; \code{spatstat} provides a sensible default.
#'
#' @return a data.frame of two columns - distance (dist) and corresponding K function values (Kfn).
#' @export
#'
#' @importFrom spatstat.explore Kinhom
#' @importFrom spatstat.geom ppp owin interp.im
#'
#' @examples
#' # Get the data
#' dat <- scampr::gorillas
#' dat$elev <- scale(dat$elevation)
#' mod <- scampr(pres ~ elev, dat, model.type = "PO", include.sre = F)
#' K_func <- kfunc(mod)
kfunc <- function(model, correction = c("border", "bord.modif", "isotropic", "translate"), point.pattern, intensity.at.pp, spatstat.win, intensity.at = c("im", "pts"), dists = NULL) {

  correction = match.arg(correction)
  intensity.at = match.arg(intensity.at)

  # Perform checks if a model is provided
  if (!missing(model)) {
    # Check the object type
    if (class(model) != "scampr") {
      stop("provided model must be a scampr")
    }
    # If the window is missing try and automatically determine one
    if (missing(spatstat.win)) {
      bnd <- bound(model)
      spatstat.win <- spatstat.geom::owin(poly = bnd)
    }
  } else if (missing(spatstat.win)) {
    spatstat.win <- NULL
  }

  # Adjust the data to be used based on the arguments provided
  if (missing(model) & missing(point.pattern)) { # BOTH MISSING
    stop("Provide either a scampr model or point.pattern")
  } else if (missing(point.pattern)) { # POINT PATTERN MISSING
    # Extract the required model objects
    pres.id <- model$pt.quad.id
    coord.names <- model$coord.names
    pres <- model$data[pres.id == 1, coord.names]
    quad <- model$data[pres.id == 0, coord.names]
    lambda <- exp(model$fitted.values)
    # Adjust the intensity if needed
    if (intensity.at == "im") {
      inten.im <- vec2im(lambda[pres.id == 0], quad[,1], quad[,2])
    }
    lambda.at.pts <- lambda[pres.id == 1]
  } else if (missing(model)) { # MODEL MISSING
    if (missing(intensity.at.pp)) { # INTENSITY MISSING
      stop("Provide the intensity at the presence points or a scampr model from which to obtain the intensity")
    } else { # INTENSITY PROVIDED
      if (nrow(point.pattern) != length(intensity.at.pp)) {
        stop("number of rows in 'point.pattern' must equal the length of the 'intensity.at.pp' vector")
      }
      coord.names <- colnames(point.pattern)[1:2]
      pres <- point.pattern[ , coord.names]
      lambda.at.pts <- intensity.at.pp
      if (intensity.at == "im") {
        intensity.at <- "pts"
        warning("intensity.at cannot equal im when calculating K function from point.pattern and intensity.at.pp")
      }
    }
  } else { # NETHER MODEL NOR POINT PATTERN MISSING (IGNORES PROVIDED INTENSITY IN THIS CASE)
    pres.id <- model$pt.quad.id
    coord.names <- model$coord.names
    quad <- model$data[pres.id == 0, coord.names]
    lambda <- exp(model$fitted.values)
    if (!all(coord.names %in% colnames(point.pattern))) {
      stop("point.pattern must be a data frame with columns named the same as model coord.names")
    }
    pres <- point.pattern[ , coord.names]
    # Adjust the intensity if needed
    inten.im <- vec2im(lambda[pres.id == 0], quad[,1], quad[,2])
    if (intensity.at == "pts") {
      if (missing(intensity.at.pp)) {
        lambda.at.pts <- spatstat.geom::interp.im(inten.im, pres[,1], pres[,2])
      } else {
        lambda.at.pts <- intensity.at.pp
      }
    }
  }

  # Create the spatstat required point pattern object
  if (is.null(spatstat.win)) {
    pres.pp <- spatstat.geom::ppp(pres[,1], pres[,2])
  } else {
    pres.pp <- spatstat.geom::ppp(pres[,1], pres[,2], window = spatstat.win)
  }

  # Check if the distances ('dists') are missing to use defaults
  # if (missing(dists)) {
  #   dists <- NULL
  # }

  res <- switch(intensity.at,
                im = spatstat.explore::Kinhom(pres.pp, lambda = inten.im, correction = correction, r = dists),
                pts = spatstat.explore::Kinhom(pres.pp, lambda = lambda.at.pts, correction = correction, r = dists)
                )

  fn <- switch(correction,
                border = res$border,
                bord.modif = res$bord.modif,
                isotropic = res$iso,
                translate = res$trans
               )

  return(cbind.data.frame(dist = res$r, Kfn = fn))
}
