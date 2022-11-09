#' Plot the image of a field over a model quadrature ### WIP ###
#'
#' @description An image plotting function that uses fields::image.plot() to display the specified variable as a field over the quadrature.
#'
#' @param x A scampr model object
#' @param z Either a single characater string of the variable name (in the model data) or one of 'fitted', 'residuals'. Alternatively, a vector of numeric values to be plotted
#' @param residual.type an optional character string for residual type if z == 'residual'
#' @param residual.smoothing an optional numeric for the scale of residual smoothing (theta in fields::image.smooth())
#' @param ... additional plotting arguments
#'
#' @return See fields::image.plot()
#' @export
#'
#' @importFrom graphics plot.default
#' @importFrom grDevices topo.colors
#' @importFrom fields image.plot rdist image.smooth as.image
#'
#' @examples
#' # Get the flora data for one of the species
#' dat_po <- flora$po$sp1
#' dat_pa <- flora$pa
#'
#' # Attach the quadrature to the PO data
#' dat_po <- rbind.data.frame(dat_po, flora$quad)
#'
#' # Fit a combined data model
#' m.comb <- scampr(pres ~ MNT + D.Main, dat_po, sp1 ~ MNT,
#' dat_pa, model.type = "ipp")
#'
#' image(m.comb, "fitted")
#' \dontrun{
#' image(m.comb, "MNT")
#' image(m.comb, "residuals")
#' }
image.scampr <- function(x, z, residual.type, residual.smoothing = 0.5, ...) {
  xtrargs <- list(...)
  # set a residual plot identifier
  is.resid <- F
  # check if model is from pa data - no image available so just plot the data
  if (x$data.model.type == "pa") {
    quad <- x$data
    resp <- quad[ , all.vars(x$formula[[2L]])] > 0
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
    xtrargs$x <- quad[ , x$coord.names[1]]
    xtrargs$y <- quad[ , x$coord.names[2]]
    # graphics::plot.default(quad[ , x$coord.names[1]], quad[ , x$coord.names[2]],
    #      col = pa.col, pch = pa.pch, asp = 1, xlab = x$coord.names[1],
    #      ylab = x$coord.names[2], main = "Presence/Absence Sites"
    # )
    do.call(graphics::plot.default, xtrargs)
    # warning("Only presence/absence survey sites are shown for this model's image")
  } else {
    # checks
    if (length(z) == 1) {
      if (typeof(z) == "character") {
        if (!(z %in% c(colnames(x$data), 'fitted', 'residuals'))) {
          stop(paste0(z, " must be one of 'fitted', 'residuals' or one of the columns of the model data frame.\nAlternatively supply a vector of z values directly"))
        } else {
          if (z %in% colnames(x$data)) {
            z.name <- z
            z <- x$data[x$pt.quad.id == 0, z.name]
          } else if (z == "fitted") {
            z.name <- "Fitted log-Intensity"
            z <- x$fitted.values[x$pt.quad.id == 0]
          } else {
            if (missing(residual.type)) {
              z.name <- "Residuals (raw)"
              tmp.z <- residuals.scampr(x, type = "raw")
              # z <- tmp.z[x$pt.quad.id ==0]
              # need to count the presence points in each quadrat for plotting
              pres <- x$data[x$pt.quad.id == 1, x$coord.names]
              quad <- x$data[x$pt.quad.id == 0, x$coord.names]
              dists <- fields::rdist(pres, quad)
              is.min <- matrix(data = F, nrow = nrow(dists), ncol = ncol(dists))
              for (i in 1:nrow(dists)) {
                is.min[i, which.min(dists[i, ])] <- T
              }
              pp.resid.in.quad <- apply(is.min, 2, sum)
              # add these into the residuals at the quadrature points
              z <- tmp.z[x$pt.quad.id ==0] + pp.resid.in.quad
              is.resid <- T
            } else {
              if (!residual.type %in% c("raw", "inverse", "pearson")) {
                stop("residual.type must be one of 'raw', 'inverse' or 'pearson'")
              } else {
                res.name <- residual.type
                substr(res.name, 1, 1) <- toupper(substr(res.name, 1, 1))
                z.name <- paste0("Residuals (", res.name, ")")
                tmp.z <- residuals.scampr(x, residual.type)
                # z <- tmp.z[x$pt.quad.id ==0]
                # need to count the presence points in each quadrat for plotting
                pres <- x$data[x$pt.quad.id == 1, x$coord.names]
                quad <- x$data[x$pt.quad.id == 0, x$coord.names]
                dists <- fields::rdist(pres, quad)
                is.min <- matrix(data = F, nrow = nrow(dists), ncol = ncol(dists))
                for (i in 1:nrow(dists)) {
                  is.min[i, which.min(dists[i, ])] <- T
                }
                # FOR NOW SEE IF THESE ARE NOT NEEDED #
                # # adjust for non-raw residuals
                # pp.resids <- matrix(0, nrow = nrow(dists), ncol = ncol(dists))
                # pp.resids[is.min] <- switch(residual.type,
                #        raw = 1,
                #        inverse = 1 / exp(x$fitted.values[x$pt.quad.id == 1]),
                #        pearson = 1 / sqrt(exp(x$fitted.values[x$pt.quad.id == 1]))
                # )
                # add.pp.resid.in.quad <- apply(pp.resids, 2, sum)
                add.pp.resid.in.quad <- apply(is.min, 2, sum)
                quad.resid.lambda <- switch(residual.type,
                                             raw = rep(1, sum(x$pt.quad.id == 0)),
                                             inverse = exp(x$fitted.values[x$pt.quad.id == 0]),
                                             pearson = sqrt(exp(x$fitted.values[x$pt.quad.id == 0]))
                                      )
                # add these into the residuals at the quadrature points
                z <- tmp.z[x$pt.quad.id ==0] + (add.pp.resid.in.quad / quad.resid.lambda)
                is.resid <- T
              }
            }
          }
        }
      }
    } else {
      z.name <- quote(z) # can't seem to get this to quote the input argument
    }
    if (length(z) != sum(x$pt.quad.id == 0)) {
      stop("z vector of field values must match the number (and order) of quadrature points in the model provided")
    }

    quad <- x$data[x$pt.quad.id == 0, ]
    # adjust for smoothing if plotting a residuals image
    if (is.resid) {
      im_field <- fields::as.image(z, x = as.matrix(cbind(quad$x, quad$y)))
      smooth.z <- fields::image.smooth(im_field, theta = residual.smoothing)
      smooth.z$z[is.na(im_field$z)] <- NA
      xtrargs$x <- smooth.z$x
      xtrargs$z <- smooth.z$z
      xtrargs$y <- smooth.z$y
      # xtrargs$residual.type <- NULL
    } else  {
      xs <- sort(unique(quad[ , x$coord.names[1]]))
      ys <- sort(unique(quad[ , x$coord.names[2]]))
      zs <- vec2mat(z, quad[ , x$coord.names[1]], quad[ , x$coord.names[2]])
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
    # fields::image.plot(x = xs, y = ys, z = zs, col = grDevices::topo.colors(100), asp = 1, xlab = x$coord.names[1], ylab = x$coord.names[2], main = z.name, bty = 'n')
    do.call(fields::image.plot, xtrargs)
  }
}
