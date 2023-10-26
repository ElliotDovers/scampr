#' Plot the image of a field over a model quadrature ### WIP ###
#'
#' @description An image plotting function that uses fields::image.plot() to display the specified variable as a field over the quadrature.
#'
#' @param x a scampr model object
#' @param z Either a single characater string of the variable name (in the model data) or one of 'fitted', 'residuals'. Alternatively, a vector of numeric values to be plotted
#' @param domain.data a data frame containing (at least) the horizontal and vertical coordinates as used in the model 'x' as well as the variable name of 'z' if specified. If missing, the quadrature points of the model will be used.
#' @param residual.type an optional character string for residual type if z == 'residual'
#' @param residual.smoothing an optional numeric for the scale of residual smoothing (theta in fields::image.smooth())
#' @param ... additional plotting arguments
#'
#' @return See fields::image.plot()
#' @export
#'
#' @importFrom graphics plot.default
#' @importFrom grDevices topo.colors
#' @importFrom fields image.plot image.smooth as.image
#'
#' @examples
#' # Get the flora data for one of the species
#' dat_po <- flora$po$sp1
#' dat_pa <- flora$pa
#'
#' # Attach the quadrature points to the presence-only data
#' dat_po <- rbind.data.frame(dat_po, flora$quad) # using full quadrature for plotting
#'
#' # Ensure the "response" variable in each data set shares the same name
#' dat_po$presence <- dat_po$pres
#' dat_pa$presence <- dat_pa$sp1
#'
#' # Integrated Data Model
#' idm <- scampr(presence ~ MNT, dat_po, bias.formula = ~ D.Main,
#' pa.data = dat_pa, include.sre = F, model.type = "IDM", latent.po.biasing = F)
#'
#' image(idm, "fitted")
#' \dontrun{
#' image(idm, "MNT")
#' image(idm, "residuals")
#' }
image.scampr <- function(x, z, domain.data, residual.type, residual.smoothing = 0.5, ...) {
  xtrargs <- list(...)

  if (!is(x, "scampr")) {
    stop("x must be a fitted scampr model")
  }
  use.quad.flag <- FALSE # indicator that the quadrature is being used as domain data
  if (missing(domain.data)) {
    if (x$model.type == "PA") {
      domain.data <- x$data
    } else {
      domain.data <- x$data[x$pt.quad.id == 0, ]
      use.quad.flag <- TRUE # indicator that the quadrature is being used as domain data
    }
  }
  if (!all(x$coord.names %in% colnames(domain.data))) {
    stop("model coordinate names are not found in the domain data")
  }

  # set a residual plot identifier
  is.resid <- F
  # check if model is from PA data - no image available so just plot the data
  if (x$model.type == "PA") {
    # quad <- domain.data
    resp <- domain.data[ , all.vars(x$formula[[2L]])] > 0
    pa.col <- resp
    pa.col[resp == 0] <- "lightblue"
    pa.col[resp == 1] <- "darkblue"
    pa.pch <- resp
    pa.pch[resp == 0] <- 4
    pa.pch[resp == 1] <- 1
    # Set the default names if not supplied
    if (!"xlab" %in% names(xtrargs)) {
      xtrargs$xlab <- x$coord.names[1]
    }
    if (!"ylab" %in% names(xtrargs)) {
      xtrargs$ylab <- x$coord.names[2]
    }
    if (!"main" %in% names(xtrargs)) {
      xtrargs$main <- "Presence/Absence Sites"
    }
    if (!"col" %in% names(xtrargs)) {
      xtrargs$col <- pa.col
    }
    if (!"pch" %in% names(xtrargs)) {
      xtrargs$pch <- pa.pch
    }
    if (!"asp" %in% names(xtrargs)) {
      xtrargs$asp <- 1
    }
    xtrargs$x <- domain.data[ , x$coord.names[1]]
    xtrargs$y <- domain.data[ , x$coord.names[2]]
    do.call(graphics::plot.default, xtrargs)
  } else {
    # checks
    if (length(z) == 1) {
      if (typeof(z) == "character") {
        if (!(z %in% c(colnames(domain.data), 'fitted', 'residuals'))) {
          stop(paste0(z, " must be one of 'fitted', 'residuals' or one of the columns of the model data frame/domain data.\nAlternatively supply a vector of z values directly"))
        } else {
          if (z %in% colnames(domain.data)) {
            z.name <- z
            z <- domain.data[ , z.name]
          } else if (z == "fitted") {
            z.name <- "Fitted log-Intensity"
            if (use.quad.flag) {
              z <- x$fitted.values[x$pt.quad.id == 0]
            } else {
              z <- predict.scampr(x, newdata = domain.data)
            }
          } else {
            is.resid <- T # residual flag
            if (missing(residual.type)) { # default to raw residuals if missing
              residual.type <- "raw"
            }
            if (!residual.type %in% c("raw", "inverse", "pearson")) {
              stop("residual.type must be one of 'raw', 'inverse' or 'pearson'")
            } else {
              # create a nice name for the residual image
              res.name <- residual.type
              substr(res.name, 1, 1) <- toupper(substr(res.name, 1, 1))
              z.name <- paste0("Residuals (", res.name, ")")
              tmp.z <- residuals.scampr(x, residual.type)

              # adjust approach according to whether domain.data was supplied
              if (!use.quad.flag) {
                # set up spatstat ppp object for fast nearest neighbour alogrithm
                domain.pp <- spatstat.geom::ppp(domain.data[, x$coord.names[1]], domain.data[, x$coord.names[2]],
                                                window = spatstat.geom::owin(xrange = range(domain.data[, x$coord.names[1]]), yrange = range(domain.data[, x$coord.names[2]])))
                # predict lambda across the new domain
                lambda.domain <- exp(predict(x, newdata = domain.data))
                # get the weights for the domain data
                if (!x$quad.weights.name %in% colnames(domain.data)) {
                  warning(paste0("Quadrature weight column used in model (", x$quad.weights.name, ") is not found in the domain data.\nAssuming domain data weights are all 1..."))
                  w <- rep(1, nrow(domain.data))
                } else {
                  w <- domain.data[ , x$quad.weights.name]
                }
              } else {
                # otherwise, just use the quadrature as the domain data
                domain.pp <- spatstat.geom::ppp(x$data[x$pt.quad.id == 0, x$coord.names[1]], x$data[x$pt.quad.id == 0, x$coord.names[2]],
                                                window = spatstat.geom::owin(xrange = range(x$data[x$pt.quad.id == 0, x$coord.names[1]], na.rm = T),
                                                                             yrange = range(x$data[x$pt.quad.id == 0, x$coord.names[2]], na.rm = T)
                                                )
                )
                # get lambda over the quadrature points
                lambda.domain <- exp(x$fitted.values[x$pt.quad.id == 0])
                # get the weights for the domain data
                w <- x$data[x$pt.quad.id == 0, x$quad.weights.name]
              }
              # calculate the residuals over the domain
              z <- switch(residual.type,
                          raw = - (w * lambda.domain),
                          inverse = - w,
                          pearson = - (w * sqrt(lambda.domain)),
              )
              # set presence points up as ppp object
              pres.pp <- spatstat.geom::ppp(x$data[x$pt.quad.id == 1, x$coord.names[1]], x$data[x$pt.quad.id == 1, x$coord.names[2]], window = domain.pp$window)
              # find the nearest neighbours of domain data to presence points
              pres.to.quad <- spatstat.geom::nncross(pres.pp, domain.pp)
              # get counts of presence points within each domain cell (that are > 0)
              quad.counts <- table(pres.to.quad$which)
              # get the index of domain points for added the presence point residuals (these are all 1)
              quad.idx <- as.numeric(attr(quad.counts , "dimnames")[[1]])
              # calculate the amounts to add to quad residuals in the denominator
              quad.add.denominator <- switch(residual.type,
                                             raw = rep(1, nrow(domain.data)),
                                             inverse = lambda.domain,
                                             pearson = sqrt(lambda.domain)
              )
              # add presence point residuals to the domain cells
              z[quad.idx] <- z[quad.idx] + (as.vector(quad.counts) / quad.add.denominator[quad.idx])
            }
          }
        }
      } else {
        z.name <- quote(z) # can't seem to get this to quote the input argument
      }
    }

    # adjust for smoothing if plotting a residuals image
    if (is.resid) {
      im_field <- fields::as.image(z, x = as.matrix(domain.data[ , x$coord.names]))
      smooth.z <- fields::image.smooth(im_field, theta = residual.smoothing)
      smooth.z$z[is.na(im_field$z)] <- NA
      xtrargs$x <- smooth.z$x
      xtrargs$z <- smooth.z$z
      xtrargs$y <- smooth.z$y
      # xtrargs$residual.type <- NULL
    } else  {
      xs <- sort(unique(domain.data[ , x$coord.names[1]]))
      ys <- sort(unique(domain.data[ , x$coord.names[2]]))
      zs <- vec2mat(z, domain.data[ , x$coord.names[1]], domain.data[ , x$coord.names[2]])
      # Enforce certain plotting elements
      xtrargs$x <- xs
      xtrargs$y <- ys
      xtrargs$z <- zs
    }

    # Set the default names if not supplied
    if (!"xlab" %in% names(xtrargs)) {
      xtrargs$xlab <- x$coord.names[1]
    }
    if (!"ylab" %in% names(xtrargs)) {
      xtrargs$ylab <- x$coord.names[2]
    }
    if (!"main" %in% names(xtrargs)) {
      xtrargs$main <- z.name
    }
    # if (!"col" %in% names(xtrargs)) {
    #   xtrargs$col <- grDevices::topo.colors(100)
    # }
    if (!"asp" %in% names(xtrargs)) {
      xtrargs$asp <- 1
    }
    xtrargs$bty <- 'n'
    do.call(fields::image.plot, xtrargs)
  }
}
