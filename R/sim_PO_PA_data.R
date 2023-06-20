sim_occurrence_data <- function(
    Intercept_po = -2,
    Intercept_pa = -2,
    Beta.env = 1,
    Beta.observer_bias = 1,
    sites.sampled = 500,
    rseed = NULL,
    env.covariate.type = c("gradient", "random_field"),
    env.covariate.range = 10,
    presence.only.observer.bias.covariate.type = c("none", "gradient", "random_field"),
    presence.only.observer.bias.covariate.range = 15,
    latent.field = F,
    latent.variance = 1,
    latent.range = 30,
    domain = list(x = 1:100, y = 1:100),
    plotting = TRUE,
    site.size = 1 # not yet implemented
    ) {
  
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
    stop("Please install 'RandomFields' package before using this function")
  }

  # match character arguments
  env.covariate.type <- match.arg(env.covariate.type)
  presence.only.observer.bias.covariate.type <- match.arg(presence.only.observer.bias.covariate.type)
  
  # store the arguments
  arg.info <- as.list(environment())
  
  # Need to set the number of pixels
  spatstat.options(npixel=c(length(unique(domain$x)), length(unique(domain$y))))
  # Set the seed if present (default is NULL => no set seed)
  set.seed(rseed)
  
  ##############################################################################
  # Interpolate some covariate at x, y locations ###############################
  interp.covar <- function(x.loc, y.loc, covar.name, data = quad){

    # turn the quadrature into a spatial pixels data frame
    sp.quad <- sp::SpatialPixelsDataFrame(points = quad[,c("x", "y")], data = quad[ , !colnames(quad) %in% c("x", "y", "quad.size")])
    
    # turn coordinates into SpatialPoints object:
    spp = sp::SpatialPoints(data.frame(x = x.loc,y = y.loc)) 
    
    # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
    v <- sp::over(spp, sp.quad[ , covar.name])
    v[is.na(v)] = 0 # NAs are a problem! Remove them
    return(v[,1])
  }
  ################################################################################
  # # Convert a vector with spatial locations to a matrix ########################
  # vec2mat <- function(vec, x.loc, y.loc){
  #   ux <- sort(unique(x.loc))
  #   uy <- sort(unique(y.loc))
  #   nx <- length(ux)
  #   ny <- length(uy)
  #   row.ref <- match(x.loc, ux)
  #   col.ref <- match(y.loc, uy)
  #   Vec <- rep(NA, max(row.ref)*max(col.ref))
  #   vec.ref <- (col.ref - 1)*max(row.ref) + row.ref
  #   Vec[vec.ref] <- vec
  #   mat <- matrix(Vec, max(row.ref), max(col.ref), dimnames = list(ux, uy))
  #   
  #   return(t(mat)) # needs to be transposed to play nice
  # }
  ###############################################################################
  # Convert a vector with spatial locations to a matrix ########################
  vec2im <- function(vec, x.loc, y.loc){
    ux <- sort(unique(x.loc))
    uy <- sort(unique(y.loc))
    nx <- length(ux)
    ny <- length(uy)
    row.ref <- match(x.loc, ux)
    col.ref <- match(y.loc, uy)
    Vec <- rep(NA, max(row.ref)*max(col.ref))
    vec.ref <- (col.ref - 1)*max(row.ref) + row.ref
    Vec[vec.ref] <- vec
    mat <- matrix(Vec, max(row.ref), max(col.ref), dimnames = list(ux, uy))
    
    return(spatstat.geom::im(t(mat), xcol = sort(unique(x.loc)), yrow = sort(unique(y.loc)))) # needs to be transposed to play nice
  }
  ##############################################################################
  # Convert matrix values to a vector to match spatial locations ###############
  mat2vec <- function(mat, x.locs, y.locs){
    grid.vec <- NULL
    xs <- NULL
    ys <- NULL
    for (i in 1:nrow(mat)) {
      grid.vec[(1:ncol(mat)) + (i-1)*ncol(mat)] <- mat[i, ]
      ys <- c(ys, 1:ncol(mat))
      xs <- c(xs, rep(i, length(1:ncol(mat))))
    }
    return(data.frame(row = xs, col = ys, value = grid.vec))
  }
  ##############################################################################
  
  # The domain is [1,100] x [1,100] on a 100 x 100 grid.  
  # grid <- list(x= seq(1, 100, by = 1), y= seq(1, 50, by = 1))
  grid <- domain
  quad <- expand.grid(grid)

  # Calculate the single environmental covariate over the grid
  if (env.covariate.type == "gradient") {
    env.covar <- as.numeric(scale(quad$x))
  } else {
    ## NOTE `fields` package GRF simulation fails too easily for ranges > 10 units on [1,100]x[1,100] 
    # try(assign("env.obj", fields::circulantEmbeddingSetup(grid,
    #                               cov.args=list(Covariance="Matern",
    #                                             aRange=env.covariate.range,
    #                                             smoothness=env.covariate.smoothness)
    # )))
    # if (!exists("env.obj")) {
    #   stop("Environmental covariate field simulation failed. Try a different combination of range and smoothness parameters")
    # }
    # env.grid <- fields::circulantEmbedding(env.obj)
    # env.covar <- mat2vec(env.grid)[ , "value"]
    ## USING RandomFields instead:
    
    # Need to block installation for RandomFields on Linux (since it is no longer on CRAN)
    RFoptions(install="no", seed = rseed + 1000) # ensures the GRF is entirely different from the latent and bias field (1000 > # simulations)
    env.cov.mod <- RMgauss(var = 1, scale = env.covariate.range)
    env.grid <- RFsimulate(model = env.cov.mod, x = grid$x, y = grid$y)
    env.covar <- env.grid$variable1
  }
  # add into the full quadrature
  quad$env <- env.covar
  
  # Calculate observer bias covariate and detection probability if required
  if (presence.only.observer.bias.covariate.type == "gradient") {
    bias.covar <- as.numeric(scale(quad$y)) # perpendicular to the env covariate so to be not correlated
  } else if (presence.only.observer.bias.covariate.type == "random_field") {
    ## NOTE `fields` package GRF simulation fails too easily for ranges > 10 units on [1,100]x[1,100]
    # try(assign("bias.obj", circulantEmbeddingSetup(grid,
    #                                cov.args=list(Covariance="Matern",
    #                                              aRange=observer.bias.covariate.range,
    #                                              smoothness=observer.bias.covariate.smoothness)
    # )))
    # if (!exists("bias.obj")) {
    #   stop("Presence-only biasing covariate field simulation failed. Try a different combination of range and smoothness parameters")
    # }
    # bias.grid <- circulantEmbedding(bias.obj)
    # bias.covar <- mat2vec(bias.grid)[ , "value"]
    ## USING RandomFields instead:
    
    # Need to block installation for RandomFields on Linux (since it is no longer on CRAN)
    RFoptions(install="no", seed = rseed + 2000) # ensures the GRF is entirely different from the latent and env. field (2000 > # simulations)
    bias.cov.mod <- RMgauss(var = 1, scale = presence.only.observer.bias.covariate.range)
    bias.grid <- RFsimulate(model = bias.cov.mod, x = grid$x, y = grid$y)
    bias.covar <- bias.grid$variable1
  } else {
    bias.covar <- rep(0, nrow(quad))
  }
  # add into the full quadrature
  quad$bias <- bias.covar
  
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
  
  # Get the linear predictor
  quad$eta.fixed <- Intercept_po + (Beta.env * quad$env) + (Beta.observer_bias * quad$bias)
  
  # Set the observation window
  wnd <- owin(xrange = range(quad$x), yrange = range(quad$y))
  
  ##############################################################################
  ## Simulate the PO data ######################################################
  ##############################################################################

  # reset the random seed for generating the latent field
  RFoptions(install="no", seed = rseed)
  
  # Simulate presence points
  if (latent.field) {
    pp <- spatstat.random::rLGCP(model = "stable",
                mu = vec2im(quad$eta.fixed, quad$x, quad$y),
                var = latent.variance, scale = latent.range, alpha = 2, win = wnd,
                saveLambda = TRUE)
    lambda_im <- attr(pp, "Lambda")
  } else {
    lambda_im <- vec2im(exp(quad$eta.fixed), quad$x, quad$y)
    pp <- spatstat.random::rpoispp(lambda = lambda_im)
  }
  
  # add the true PO intensity to the quadrature (lambda)
  quad$lambda <- mat2vec(lambda_im$v)[ , "value"]
  
  # calculate the true mean abundance rate (mu) - remove the PO bias and PO intercept, then add in the PA intercept
  tmp.xi <- log(quad$lambda) - quad$eta.fixed
  quad$mu <- exp(Intercept_pa + (Beta.env * quad$env) + tmp.xi)
  
  # calculate the true shared latent field (xi)
  quad$xi <- tmp.xi
  
  # calculate the presence probability across the domain
  quad$pres.prob <- 1 - exp(-quad$mu)
  
  # make a data frame for the presence points
  pp <- data.frame(x = pp$x, y = pp$y)
  
  # Interpolate the various fields onto the presence-only data #
  
  # interpolate the environmental variable
  if (env.covariate.type == "gradient") {
    pp$env <- (pp$x - attr(scale(quad$x), "scaled:center")) / attr(scale(quad$x), "scaled:scale")
  } else {
    pp$env <- interp.covar(x = pp$x, y = pp$y, covar.name = "env")
  }
  # interpolate the bias variable
  if (presence.only.observer.bias.covariate.type == "gradient") {
    pp$bias <- (pp$y - attr(scale(quad$y), "scaled:center")) / attr(scale(quad$y), "scaled:scale")
  } else {
    pp$bias <- interp.covar(x = pp$x, y = pp$y, covar.name = "bias")
  }
  
  # add in the quadrature size (zero at the presence points)
  pp$quad.size <- 0
  # interpolate the linear predictor
  pp$eta.fixed <- interp.covar(x = pp$x, y = pp$y, covar.name = "eta.fixed")
  # interpolate the true PO intensity
  pp$lambda <- interp.covar(x = pp$x, y = pp$y, covar.name = "lambda")
  # interpolate the true abundance rate
  pp$mu <- interp.covar(x = pp$x, y = pp$y, covar.name = "mu")
  # interpolate the true latent field
  pp$xi <- interp.covar(x = pp$x, y = pp$y, covar.name = "xi")
  # interpolate the true presence probability
  pp$pres.prob <- interp.covar(x = pp$x, y = pp$y, covar.name = "pres.prob")
  
  # add in a presence identifier
  pp$present <- 1
  
  ##############################################################################
  ## Abundance Counts and Presence/Absence Data ################################
  ##############################################################################
  
  # randomly sample survey sites
  survey <- data.frame(x = runif(sites.sampled, min(domain$x), max(domain$x)), y = runif(sites.sampled, min(domain$y), max(domain$y)))

  # Interpolate the various fields #
  
  # interpolate the environmental variable
  if (env.covariate.type == "gradient") {
    survey$env <- (survey$x - attr(scale(quad$x), "scaled:center")) / attr(scale(quad$x), "scaled:scale")
  } else {
    survey$env <- interp.covar(x = survey$x, y = survey$y, covar.name = "env")
  }
  # calculate the linear predictor
  survey$eta.fixed <- Intercept_pa + (Beta.env * survey$env)
  # calculate the true abundance rate (via interpolating the latent field)
  xi_pa <- interp.covar(x = survey$x, y = survey$y, covar.name = "xi")
  survey$mu <- exp(survey$eta.fixed + xi_pa)
  # add in the latent field
  survey$xi <- xi_pa
  # calculate the presence probability
  survey$pres.prob <-  1 - exp(-survey$mu)
  
  # add in a presence/absence
  survey$present <- rbinom(nrow(survey), 1, survey$pres.prob)

  ##############################################################################
  # Check out the field characteristics if required
  if (plotting) {
    par(mar = c(0, 0, 1, 0))
      if (latent.field & presence.only.observer.bias.covariate.type != "none") {
        m <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2, byrow = TRUE)
        layout(mat = m, heights = rep(1/3, 3), widths = c(0.5,0.5))
        plot(vec2im(quad$env, quad$x, quad$y), main = "Env. Covariate", col = terrain.colors(256))
        plot(vec2im(quad$bias, quad$x, quad$y), main = "PO Bias Covariate", col = cividis(256))
        plot(vec2im(quad$xi, quad$x, quad$y), main = "Latent Field", col = terrain.colors(256))
        plot(vec2im(log(quad$mu), quad$x, quad$y), main = "log-Abundance Rate", col = terrain.colors(256))
        plot(vec2im(quad$pres.prob, quad$x, quad$y), main = "Presence Prob. w. presence/absence")
        points(survey$x, survey$y, col = survey$present, pch = c(4, 1)[survey$present + 1])
        plot(vec2im(log(quad$lambda), quad$x, quad$y), main = "Presence record log-Intensity")
        points(pp[,c("x", "y")])
      } else if (latent.field & presence.only.observer.bias.covariate.type == "none") {
        m <- matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, ncol = 2, byrow = TRUE)
        layout(mat = m, heights = rep(1/3, 3), widths = c(0.5,0.5))
        plot(vec2im(quad$env, quad$x, quad$y), main = "Env. Covariate", col = terrain.colors(256))
        plot(vec2im(quad$xi, quad$x, quad$y), main = "Latent Field", col = terrain.colors(256))
        plot(vec2im(log(quad$mu), quad$x, quad$y), main = "log-Abundance Rate")
        plot(vec2im(quad$pres.prob, quad$x, quad$y), main = "Presence Prob. w. presence/absence")
        points(survey$x, survey$y, col = survey$present, pch = c(4, 1)[survey$present + 1])
        plot(vec2im(log(quad$lambda), quad$x, quad$y), main = "Presence record log-Intensity")
        points(pp[,c("x", "y")])
      } else if (!latent.field & presence.only.observer.bias.covariate.type != "none") {
        m <- matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, ncol = 2, byrow = TRUE)
        layout(mat = m, heights = rep(1/3, 3), widths = c(0.5,0.5))
        plot(vec2im(quad$env, quad$x, quad$y), main = "Env. Covariate", col = terrain.colors(256))
        plot(vec2im(quad$bias, quad$x, quad$y), main = "PO Bias Covariate", col = cividis(256))
        plot(vec2im(log(quad$mu), quad$x, quad$y), main = "log-Abundance Rate")
        plot(vec2im(quad$pres.prob, quad$x, quad$y), main = "Presence Prob. w. presence/absence")
        points(survey$x, survey$y, col = survey$present, pch = c(4, 1)[survey$present + 1])
        plot(vec2im(log(quad$lambda), quad$x, quad$y), main = "Presence record log-Intensity")
        points(pp[,c("x", "y")])
      } else {
        m <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
        layout(mat = m, heights = c(0.5,0.5), widths = c(0.5,0.5))
        plot(vec2im(quad$env, quad$x, quad$y), main = "Env. Covariate", col = terrain.colors(256))
        plot(vec2im(log(quad$mu), quad$x, quad$y), main = "log-Abundance Rate")
        plot(vec2im(quad$pres.prob, quad$x, quad$y), main = "Presence Prob. w. presence/absence")
        points(survey$x, survey$y, col = survey$present, pch = c(4, 1)[survey$present + 1])
        plot(vec2im(log(quad$lambda), quad$x, quad$y), main = "Presence record log-Intensity")
        points(pp[,c("x", "y")])
      }
    par(mar = c(5.1, 4.1, 4.1, 2.1))
    layout(matrix(1, nrow = 1, ncol = 1))
  }
  ##############################################################################
  
  # remove the random seed
  set.seed(NULL)
  
  res <- survey
  attr(res, "presence-only") <- pp
  attr(res, "truth.grid") <- quad
  
  # add information about the simulation
  attr(res, "sim.info") <- arg.info
  return(res)
}
