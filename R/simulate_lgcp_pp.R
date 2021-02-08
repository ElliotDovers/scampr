#' Simulate a point pattern from log-Gaussian Cox Process - LEGACY FN: INTERRNAL FOR USE IN THESIS
#'
#' @description Simulate a point pattern from a LGCP with an intercept and slope on a single deterministic spatial covariate. Over a 1 hectare domain, includes quadrature. Control the correlation range and marginal variance of the latent field.
#'
#' @param Intercept Fixed effect intercept - used to control the expected number of points
#' @param Slope Fixed effect slope - controls the effect of the spatial covariate on the point pattern
#' @param latent.variance The marginal variance of the latent Gaussian field
#' @param latent.range The correlation range of the latent Gaussian field
#' @param rseed integer for setting the random seed
#'
#' @return the simulated data frame with attributes containing sim information/specifications.
#' @noRd
#'
#' @importFrom spatstat spatstat.options im owin rLGCP
#' @importFrom RandomFields RFoptions
#'
#' @examples
#' # Simulate a point pattern
#' dat <- simulate_lgcp_pp(Intercept = -3, Slope = 1.25, latent.range = 15, rseed = 1)
#'
#' # Fit an IPP model to the point pattern
#' m.ipp <- ipp(pres ~ X, data = dat)
#'
#' # Set up a simple 2D grid of basis functions to fit a LGCP model to the data
#' bfs <- simple_basis(nodes.on.long.edge = 10, data = dat)
#'
#' # Fit a LGCP model to the point pattern
#' m.lgcp_va <- lgcp(pres ~ X, data = dat, approx.with = "variational", simple.basis = bfs)
#'
#' summary(m.ipp)
#' summary(m.lgcp_va)
simulate_lgcp_pp <- function(Intercept = -3, Slope = 1.25, latent.variance = 1, latent.range = 15, rseed = NA) {

  # Checks #
  # if (!library(spatstat, logical.return = T)) {
  #   stop("Please install 'spatstat' package before using this function")
  # }
  # if (!library(fields, logical.return = T)) {
  #   stop("Please install 'fields' package before using this function")
  # }
  # if (!library(sp, logical.return = T)) {
  #   stop("Please install 'sp' package before using this function")
  # }

  # Need to set the number of pixels
  spatstat::spatstat.options(npixel=c(101, 101))

  ############################################################################################################################
  # Convert a vector to a spatstat image object via vector locations #####(mainly for plotting) ##############################
  vec2im <- function(vec, x.loc, y.loc){
    ux <- sort(unique(x.loc))
    uy <- sort(unique(y.loc))
    nx <- length(ux)
    ny <- length(uy)
    col.ref <- match(x.loc, ux)
    row.ref <- match(y.loc, uy)
    Vec <- rep(NA, max(row.ref)*max(col.ref))
    vec.ref <- (col.ref - 1)*max(row.ref) + row.ref
    Vec[vec.ref] <- vec
    grid.mask <- matrix(Vec, max(row.ref), max(col.ref),dimnames = list(uy, ux))
    vec.im <- spatstat::im(grid.mask, xcol = ux, yrow = uy)

    return(vec.im)
  }
  ############################################################################################################################
  # Interpolate some covariate at x, y locations #############################################################################
  f.env_covar <- function(x, y) {
      sin((x / 20) - 15) - cos((y / 20))
    }
  #plot(vec2im(f.env_covar(quad$x, quad$y), quad$x, quad$y))

  # This interpolation is for the realised log-lambda
  interp_fld <- function(x.loc, y.loc, covar.name, dframe){
    if (!covar.name %in% colnames(dframe)) {
      stop(paste0(covar.name, " is not a column of ", dframe))
    }
    sp.fld <- sp::SpatialPixelsDataFrame(points = dframe[,c("x", "y")], data = dframe[ , !colnames(dframe) %in% c("x", "y")])

    # turn coordinates into SpatialPoints object:
    spp = sp::SpatialPoints(data.frame(x = x.loc,y = y.loc))

    # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
    v <- sp::over(spp, sp.fld[ , covar.name])
    v[is.na(v)] = 0 # NAs are a problem! Remove them
    return(v[,1])
  }
  ############################################################################################################################

  # The domain is [0,100] x [0,100] on a 101 x 101 grid.
  grid <- list(x= seq(0, 100, by = 1), y= seq(0, 100, by = 1))
  # Place in long form in a data frame
  quad <- expand.grid(grid)
  quad$X <- f.env_covar(quad$x, quad$y)
  X.grid <- matrix(nrow = 101, ncol = 101)
  for (xs in 1:101) {
    for (ys in 1:101) {
      X.grid[ys, xs] <- quad$X[quad$x == xs - 1 & quad$y == ys - 1]
    }
  }

  ## Set up the quadrat sizes
  dx <- unique(diff(sort(unique(quad$x))))[1] # [1] required due to machine error
  dy <- unique(diff(sort(unique(quad$y))))[1] # [1] required due to machine error
  D.size <- diff(range(quad$x)) * diff(range(quad$y)) # Is one hectare
  ## Quadrature Sizes ##
  # rectangles centered at the quadrature
  quad$quad.size <- rep(dx * dy, nrow(quad))
  # vertical half of rectangle - along the edges
  quad$quad.size[quad$x == min(quad$x) | quad$x == max(quad$x)] <- (dx / 2) * dy
  # horizontal half of rectangle
  quad$quad.size[quad$y == min(quad$y) | quad$y == max(quad$y)] <- (dy / 2) * dx
  # small corner rectangles
  quad$quad.size[(quad$x == min(quad$x) | quad$x == max(quad$x)) & (quad$y == min(quad$y) | quad$y == max(quad$y))] <- (dx / 2) * (dy / 2)

  # Warning - if the linear predictor
  eta.fixed <- Intercept + (Slope * quad$X)
  if (sum(quad$quad.size * exp(eta.fixed)) < 100) {
    warning(paste0("Small linear predictor: Expected N = ", sum(quad$quad.size * exp(eta.fixed))))
  }
  # Set a seed as required
  if (!is.na(rseed)) {
    # RNGkind(sample.kind = "Rounding")
    set.seed(rseed)
    RandomFields::RFoptions(seed=rseed, printlevel = 0)
  }

  # Set the observation window
  wnd <- spatstat::owin(xrange = c(0, 100), yrange = c(0, 100))

  # Simulate the LGCP
  pp <- spatstat::rLGCP(model = "stable",
              mu =
                Intercept +
                (Slope * spatstat::im(X.grid, xcol = grid$x, yrow = grid$y)),
              var = latent.variance, scale = latent.range, alpha = 2,
              win = wnd,
              saveLambda = TRUE)
  lambda.im <- attr(pp, "Lambda")
  lambda.vec <- NULL
  xs <- NULL
  ys <- NULL
  for (i in 1:lambda.im$dim[1]) {
    lambda.vec[(1:lambda.im$dim[2]) + (i-1)*lambda.im$dim[2]] <- lambda.im[i, ]
    xs <- c(xs, grid$x)
    ys <- c(ys, rep(grid$y[i], length(grid$y)))
  }
  id <- match(paste(quad$x, quad$y, sep = ":"), paste(xs, ys, sep = ":"))
  quad$loglamb <- log(lambda.vec[id])
  quad$xi <- quad$loglamb - eta.fixed

  # # Check out the sim characteristics
  # par(mar = c(0, 0, 1, 0))
  # par(mfrow = c(2, 2))
  # plot(vec2im(quad$X, quad$x, quad$y), main = paste0(Slope, " x Env. Cov."))
  # plot(vec2im(quad$loglamb, quad$x, quad$y), main = "log-Intensity")
  # plot(vec2im(quad$xi, quad$x, quad$y), main = "Latent Field")
  # plot(pp, main = "Presence Pts.")
  # par(mar = c(5.1, 4.1, 4.1, 2.1))
  # par(mfrow = c(1, 1))

  # Set up the pres frame
  pres <- data.frame(cbind(x = pp$x, y = pp$y))
  # Calculate the covariate values at the presence points
  pres$X <- f.env_covar(pres$x, pres$y)
  # Interpolate the intensity and latent fields
  pres$loglamb <- interp_fld(x = pres$x, y = pres$y, covar.name = "loglamb", dframe = quad)
  pres$xi <- interp_fld(x = pres$x, y = pres$y, covar.name = "xi", dframe = quad)
  # Set quadrat size to zero
  pres$quad.size <- 0

  # Add the presence-quadrature identifier to the data frames and combine
  dat <- rbind(cbind(pres, pres = 1), cbind(quad, pres = 0))
  # include the simulation info
  attr(dat, "sim_info") <- c(Intercept = Intercept,
                             Slope = Slope,
                             latent.variance = latent.variance,
                             latent.range = latent.range,
                             rseed = rseed,
                             N = pp$n,
                             Expected_N = sum(quad$quad.size * exp(eta.fixed))
                             )
  return(dat)
}
