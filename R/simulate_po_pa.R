#' Simulate presence-only data (point pattern from log-Gaussian Cox Process) and presence/absence data (from the same underlying process)
#'
#' @param Intercept.env Fixed effect intercept
#' @param Intercept.bias Bias effect intercept
#' @param Beta length 2 vector of env. effects
#' @param Tau length 2 vector of bias effects
#' @param latent.variance variance for the latent field
#' @param latent.range correlation range for the latent field
#' @param sites.sampled number of survey sites to be sampled
#' @param rseed integer for setting the random seed
#'
#' @return
#' @export
#'
#' @examples
simulate_po_pa <- function(
  Intercept.env = -1.75,
  Intercept.bias = -2,
  Beta = c(-1.2, 0.75),
  Tau = c(1.3, -0.8),
  latent.variance = 1,
  latent.range = 30,
  sites.sampled = 1000,
  rseed = NA) {

  # Checks #
  if (!library(spatstat, logical.return = T)) {
    stop("Please install 'spatstat' package before using this function")
  }
  if (!library(fields, logical.return = T)) {
    stop("Please install 'fields' package before using this function")
  }
  if (!library(sp, logical.return = T)) {
    stop("Please install 'sp' package before using this function")
  }
  if (!library(RandomFields, logical.return = T)) {
    stop("Please install 'sp' package before using this function")
  }

  # Need to set the number of pixels
  spatstat.options(npixel=c(100, 100))

  ############################################################################################################################
  # Convert a vector to a spatstat image object via vector locations #####(mainly for plotting) ##############################
  vec2im <- function(vec, x.loc, y.loc){
    #### spatstat package is needed ####
    if (!library(spatstat, logical.return = T)) {
      stop("Please install 'spatstat' package before using this function")
    }
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
  # Interpolate some covariate at x, y locations ###########################################
  interp.covar <- function(x.loc, y.loc, covar.name, data = quad){
    #### sp package is needed ####
    if (!library(sp, logical.return = T)) {
      stop("Please install 'sp' package before using this function")
    }
    sp.quad <- SpatialPixelsDataFrame(points = quad[,c("x", "y")], data = quad[ , !colnames(quad) %in% c("x", "y", "quad.size")])

    # turn coordinates into SpatialPoints object:
    spp = SpatialPoints(data.frame(x = x.loc,y = y.loc))

    # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
    v <- over(spp, sp.quad[ , covar.name])
    v[is.na(v)] = 0 # NAs are a problem! Remove them
    return(v[,1])
  }
  ############################################################################################################################

  # The domain is [1,100] x [1,100] on a 100 x 100 grid.
  grid <- list(x= seq(1, 100, by = 1), y= seq(1, 100, by = 1))
  quad <- expand.grid(grid)
  if (is.na(rseed)) {
    add.seed <- 0
  } else {
    add.seed <- rseed
  }
  # Set the spatial covariates (to be used as either environmental or bias variables)
  RFoptions(seed=12345 + add.seed)
  cov.mod <- RMgauss(var = 1, scale = 20)
  # ec1 <- RFsimulate(model = cov.mod, x = grid$x, y = grid$y)
  tmp <- RFsimulate(model = cov.mod, x = grid$x, y = grid$y)
  ec1 <- tmp$variable1
  RFoptions(seed=12346 + add.seed)
  cov.mod <- RMgauss(var = 1, scale = 10)
  # ec2 <- RFsimulate(model = cov.mod, x = grid$x, y = grid$y)
  tmp <- RFsimulate(model = cov.mod, x = grid$x, y = grid$y)
  ec2 <- tmp$variable1
  RFoptions(seed=12347 + add.seed)
  cov.mod <- RMgauss(var = 1, scale = 20)
  # bc1 <- RFsimulate(model = cov.mod, x = grid$x, y = grid$y)
  tmp <- RFsimulate(model = cov.mod, x = grid$x, y = grid$y)
  bc1 <- tmp$variable1
  RFoptions(seed=12348 + add.seed)
  cov.mod <- RMgauss(var = 1, scale = 10)
  # bc2 <- RFsimulate(model = cov.mod, x = grid$x, y = grid$y)
  tmp <- RFsimulate(model = cov.mod, x = grid$x, y = grid$y)
  bc2 <- tmp$variable1
  RFoptions(seed=NA)

  # par(mfrow = c(2, 2))
  # plot(spatstat::im(ec1, xcol = grid$x, yrow = grid$y), main = "Env. Cov. 1")
  # plot(spatstat::im(ec2, xcol = grid$x, yrow = grid$y), main = "Env. Cov. 2")
  # plot(spatstat::im(bc1, xcol = grid$x, yrow = grid$y), main = "Bias Cov. 1")
  # plot(spatstat::im(bc2, xcol = grid$x, yrow = grid$y), main = "Bias Cov. 2")
  # par(mfrow = c(1, 1))

  # # Place in long form in a data frame
  # quad <- expand.grid(grid)
  # ec1.vec <- NULL
  # bc1.vec <- NULL
  # ec2.vec <- NULL
  # bc2.vec <- NULL
  # xs <- NULL
  # ys <- NULL
  # for (i in 1:length(grid$x)) {
  #   ec1.vec[(1:length(grid$y)) + (i-1)*length(grid$y)] <- ec1[i, ]
  #   bc1.vec[(1:length(grid$y)) + (i-1)*length(grid$y)] <- bc1[i, ]
  #   ec2.vec[(1:length(grid$y)) + (i-1)*length(grid$y)] <- ec2[i, ]
  #   bc2.vec[(1:length(grid$y)) + (i-1)*length(grid$y)] <- bc2[i, ]
  #   xs <- c(xs, grid$x)
  #   ys <- c(ys, rep(grid$y[i], length(grid$y)))
  # }
  # id <- match(paste(quad$x, quad$y, sep = ":"), paste(xs, ys, sep = ":"))
  # quad$ec1 <- ec1.vec[id]
  # quad$bc1 <- bc1.vec[id]
  # quad$ec2 <- ec2.vec[id]
  # quad$bc2 <- bc2.vec[id]

  quad$ec1 <- ec1
  quad$bc1 <- bc1
  quad$ec2 <- ec2
  quad$bc2 <- bc2

  ## Set up the quadrat sizes
  dx <- unique(diff(sort(unique(quad$x))))[1] # [1] required due to machine error
  dy <- unique(diff(sort(unique(quad$y))))[1] # [1] required due to machine error
  D.size <- diff(range(quad$x)) * diff(range(quad$y))
  ## Quadrature Sizes ##
  # rectangles centered at the quadrature
  quad$quad.size <- rep(dx * dy, nrow(quad))
  # vertical half of rectangle - along the edges
  quad$quad.size[quad$x == min(quad$x) | quad$x == max(quad$x)] <- (dx / 2) * dy
  # horizontal half of rectangle
  quad$quad.size[quad$y == min(quad$y) | quad$y == max(quad$y)] <- (dy / 2) * dx
  # small corner rectangles
  quad$quad.size[(quad$x == min(quad$x) | quad$x == max(quad$x)) & (quad$y == min(quad$y) | quad$y == max(quad$y))] <- (dx / 2) * (dy / 2)

  # Set a seed as required
  if (!is.na(rseed)) {
    set.seed(rseed)
    RFoptions(seed=rseed)
  }
  # Get the linear predictor
  eta.fixed <- Intercept.env + (Beta[1] * quad$ec1) + (Beta[2] * quad$ec2) + Intercept.bias + (Tau[1] * quad$bc1) + (Tau[2] * quad$bc2)

  # Set the observation window
  wnd <- owin(xrange = c(1, 100), yrange = c(1, 100))
  #
  # # Simulate the LGCP - FOR SOME REASON THIS WAY FLIPS THE COVARIATES AND FUCKS SHIT UP!
  # pp <- rLGCP(model = "stable",
  #             mu =
  #               Intercept.env +
  #               (Beta[1] * spatstat::im(ec1, xcol = grid$x, yrow = grid$y)) +
  #               (Beta[2] * spatstat::im(ec2, xcol = grid$x, yrow = grid$y)) +
  #               Intercept.bias +
  #               (Tau[1] * spatstat::im(bc1, xcol = grid$x, yrow = grid$y)) +
  #               (Tau[2] * spatstat::im(bc2, xcol = grid$x, yrow = grid$y)),
  #             var = latent.variance, scale = latent.range, alpha = 2, win = wnd,
  #             saveLambda = TRUE)
  # Simulate the LGCP
  pp <- rLGCP(model = "stable",
              mu =
                Intercept.env +
                (Beta[1] * vec2im(ec1, quad$x, quad$y)) +
                (Beta[2] * vec2im(ec2, quad$x, quad$y)) +
                Intercept.bias +
                (Tau[1] * vec2im(bc1, quad$x, quad$y)) +
                (Tau[2] * vec2im(bc2, quad$x, quad$y)),
              var = latent.variance, scale = latent.range, alpha = 2, win = wnd,
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
  quad$logabun <- (Intercept.env + (Beta[1] * quad$ec1) + (Beta[2] * quad$ec2)) + quad$xi

  # # Check out the field characteristics
  # par.default <- par()
  # par(mar = c(0, 0, 1, 0))
  # m <- matrix(c(1, 2, 3, 4, 5, 6, 7, 7), nrow = 4, ncol = 2, byrow = TRUE)
  # layout(mat = m)
  # plot(vec2im(quad$ec1, quad$x, quad$y), main = paste0(Beta[1], " x Env. Cov. 1"))
  # plot(vec2im(quad$bc1, quad$x, quad$y), main = paste0(Tau[1], " x Bias Cov. 1"))
  # plot(vec2im(quad$ec2, quad$x, quad$y), main = paste0(Beta[2], " x Env. Cov. 2"))
  # plot(vec2im(quad$bc2, quad$x, quad$y), main = paste0(Tau[2], " x Bias Cov. 2"))
  # plot(vec2im(quad$logabun, quad$x, quad$y), main = "log-Abundance")
  # plot(vec2im(quad$loglamb, quad$x, quad$y), main = "log-Intensity")
  # plot(vec2im(quad$xi, quad$x, quad$y), main = "Latent Field")
  # par(mar = c(5.1, 4.1, 4.1, 2.1))
  # layout(matrix(1, nrow = 1, ncol = 1))

  # Set up the presence points
  pres <- data.frame(cbind(x = pp$x, y = pp$y))
  # Calculate the covariate values at the presence points by assigning the value of the quadrature into which it falls
  pres$ec1 <- interp.covar(x = pres$x, y = pres$y, covar.name = "ec1")
  pres$bc1 <- interp.covar(x = pres$x, y = pres$y, covar.name = "bc1")
  pres$ec2 <- interp.covar(x = pres$x, y = pres$y, covar.name = "ec2")
  pres$bc2 <- interp.covar(x = pres$x, y = pres$y, covar.name = "bc2")
  pres$loglamb <- interp.covar(x = pres$x, y = pres$y, covar.name = "loglamb")
  pres$xi <- interp.covar(x = pres$x, y = pres$y, covar.name = "xi")
  pres$logabun <- interp.covar(x = pres$x, y = pres$y, covar.name = "logabun")
  pres$quad.size <- 0
  # Add the presence-quadrature identifier to the data frames and combine
  dat_po <- rbind(cbind(pres, pres = 1), cbind(quad, pres = 0))

  # Create the presence/absence data - sampling at 50 random sites within the domain
  dat_pa <- as.data.frame(cbind(x = runif(sites.sampled, 1, 100), y = runif(sites.sampled, 1, 100)))
  dat_pa$ec1 <- interp.covar(x = dat_pa$x, y = dat_pa$y, covar.name = "ec1")
  dat_pa$bc1 <- interp.covar(x = dat_pa$x, y = dat_pa$y, covar.name = "bc1")
  dat_pa$ec2 <- interp.covar(x = dat_pa$x, y = dat_pa$y, covar.name = "ec2")
  dat_pa$bc2 <- interp.covar(x = dat_pa$x, y = dat_pa$y, covar.name = "bc2")
  dat_pa$loglamb <- interp.covar(x = dat_pa$x, y = dat_pa$y, covar.name = "loglamb")
  dat_pa$xi <- interp.covar(x = dat_pa$x, y = dat_pa$y, covar.name = "xi")
  dat_pa$logabund <- interp.covar(x = dat_pa$x, y = dat_pa$y, covar.name = "logabun")
  dat_pa$pprob <- 1 -exp(-exp(dat_pa$logabun))
  dat_pa$pa <- rbinom(length(dat_pa$pprob), 1, dat_pa$pprob)

  # Compile both data sets and information about the simulaton into a result list
  temp.info = c(Intercept.env, Intercept.bias, Beta, Tau, latent.variance, latent.range, rseed, N_po = pp$n, N_pa_pres = sum(dat_pa$pa), N_pa_sites = sites.sampled)
  names(temp.info) <- c("Intercept_env", "Intercept_bias", "ec1", "ec2", "bc1", "bc2", "latent.variance", "latent.range", "rseed", "N_po", "N_pa_pres", "N_pa_sites")
  res.list <- list(
    PO = dat_po,
    PA = dat_pa,
    info =   temp.info
  )
  return(res.list)
}
