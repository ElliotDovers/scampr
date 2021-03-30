#' Inhomogeneous K function for scampr models (is a wrapper of \code{spatstat}'s version). Or calculate the K function from point pattern (and intensity) provided.
#'
#' @param model a scampr model object
#' @param correction a character string describing the correction type to be used. See \code{spatstat::Kinhom()} for details.
#' @param point.pattern Optionally, a data frame, the first two columns of which describe the horizontal and vertical coordinates of point locations respectively.
#' @param intensity.at.pp Optionally, a vector of length \code{nrow(point.pattern)} that describes the intensity at each point.
#' @param spatstat.win Optionally, a spatstat window object.
#' @param intensity.at Optionally, a character string, one of either 'im' or 'pts'. Describes whether the intensity used in \code{spatstat::Kinhom()} is an image object or just provided at the presence points.
#'
#' @return a data.frame of two columns - distance (dist) and corresponding K function values (Kfn).
#' @export
#'
#' @importFrom spatstat owin ppp Kinhom interp.im
#'
#' @examples
#' # Get the data
#' dat <- scampr::gorillas
#' dat$elev <- scale(dat$elevation)
#' mod <- po(pres ~ elev, dat, model.type = "ipp")
#' K_func <- kfunc(mod)
kfunc <- function(model, correction = c("border", "bord.modif", "isotropic", "translate"), point.pattern, intensity.at.pp, spatstat.win, intensity.at = c("im", "pts")) {

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
      spatstat.win <- spatstat::owin(poly = bnd)
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
      lambda.at.pts <- spatstat::interp.im(inten.im, pres[,1], pres[,2])
    }
  }

  # Create the spatstat required point pattern object
  if (is.null(spatstat.win)) {
    pres.pp <- spatstat::ppp(pres[,1], pres[,2])
  } else {
    pres.pp <- spatstat::ppp(pres[,1], pres[,2], window = spatstat.win)
  }

  res <- switch(intensity.at,
                im = spatstat::Kinhom(pres.pp, lambda = inten.im, correction = correction),
                pts = spatstat::Kinhom(pres.pp, lambda = lambda.at.pts, correction = correction)
                )

  fn <- switch(correction,
                border = res$border,
                bord.modif = res$bord.modif,
                isotropic = res$iso,
                translate = res$trans
               )

  return(cbind.data.frame(dist = res$r, Kfn = fn))
}
