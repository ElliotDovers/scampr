#' Calculate the approximate boundary of 2D Euclidean data
#'
#' @param object Either a data frame containing the spatial locations of each point or a \code{scampr} model object. In the latter case the model's data will be used.
#' @param coord.names a vector of length 2 containing character strings describing the column names of the coordinates in any data provided. First coordinate name should refer to the horizontal axis. Will be over-ridden when \code{object} is a \code{scampr} model.
#' @param precision an integer value that determines the number of decimal places to search through. Controls the precision of the approximate boundary, 1 or 2 should be adequate. The higher (larger) the presicion the slower the result is obtained.
#' @param plotting a logical indicating whether the results should be plotted. Helpful for trialling presicion and/or encasing.
#' @param encase a numeric value indicating the amount of encasing of the data should occur. This will essentially add a buffer to the boundary.
#' @param smoothed a logical indicating whether or not the boundary should be smoothed using a Loess' smoother.
#'
#' @return a data frame of 2 columns for each coordinate.
#' @export
#'
#' @importFrom stats loess
#' @importFrom graphics points lines
#'
#' @examples
#' # Get the gorilla nesting data
#' dat <- gorillas
#'
#' # Calculate a boundary for the data
#' bnd <- bound(dat, coord.names = c("x", "y"))
bound <- function (object, coord.names = c("x", "y"), precision = 1, plotting = F, encase = NA, smoothed = F) {

  # Check the object type
  if (class(object) == "scampr") {
    dat <- object$data[object$pt.quad.id == 0, ]
    coord.names <- object$coord.names
  } else {
    dat <- object
  }

  # check coord.names are present in the data
  if (!all(coord.names %in% colnames(dat))) {
    stop(paste0("coord.names, ", coord.names, ", not found in data set provided"))
  }

  # Obtain functional dataframe
  A <- data.frame(cbind(X = dat[ , coord.names[1]], Y = dat[ , coord.names[2]]))

  # Find data "center"
  mid <- sapply(A, mean)

  # Initialise plot if necessary
  if (plotting) {
    plot.default(A$X, A$Y, pch = ".", xlab = coord.names[1], ylab = coord.names[2], asp = 1)
    graphics::points(mid["X"], mid["Y"], pch = 16, col = "blue")
  }

  # Calculate origin shifted coordinates
  A$x <- A$X - mid["X"]
  A$y <- A$Y - mid["Y"]

  # Remove the origin (if present)
  A <- A[!(A$x == 0 & A$y == 0), ]

  # Calculate polar coordinates
  A$theta <- atan2(A$y, A$x)
  A$r <- sqrt((A$x)^2 + (A$y)^2)

  # determine the 4 origin extremes
  ey <- range(A[A$x == A$x[which.min(abs(A$X - mid["X"]))], "y"]) # THESE ARE PROBLEMS (LOOKS LIKE I HAVEN"T ACCOUNTED FOR A SINGLE VALUE BEING FOUND)
  ex <- range(A[A$y == A$y[which.min(abs(A$Y - mid["Y"]))], "x"])

  # Determine the bandwidth of theta to loop through
  loop <- seq(-round(pi, precision), round(pi, precision), by = 1 / (10^precision))

  # Define the output
  ext <- NULL

  # Loop through theta and identify and retain the points furthest from the center of the data - axis cases handled separately
  for (i in 1:length(loop)) {
    if (loop[i] == round(-pi, precision)) {
      tempA <- data.frame(cbind(X = ex[1] + mid["X"], Y = 0 + mid["Y"], x = ex[1], y = 0, theta = round(-pi, precision), r = abs(ex[1])))
    } else if (loop[i] == round(-pi/2, precision)) {
      tempA <- data.frame(cbind(X = 0 + mid["X"], Y = ey[1] + mid["Y"], x = 0, y = ey[1], theta = round(-pi/2, precision), r = abs(ey[1])))
    } else if (loop[i] == round(0, precision)) {
      tempA <- data.frame(cbind(X = ex[2] + mid["X"], Y = 0 + mid["Y"], x = ex[2], y = 0, theta = round(0, precision), r = abs(ex[2])))
    } else if (loop[i] == round(pi/2, precision)) {
      tempA <- data.frame(cbind(X = 0 + mid["X"], Y = ey[2] + mid["Y"], x = 0, y = ey[2], theta = round(pi/2, precision), r = abs(ey[2])))
    } else if (loop[i] == round(pi, precision)) {
      tempA <- data.frame(cbind(X = ex[1] + mid["X"], Y = 0 + mid["Y"], x = ex[1], y = 0, theta = round(-pi, precision), r = abs(ex[1])))
    } else {
      tempA <- A[(A$theta >= (loop[i] - (1 / (10^precision)) / 2)) & (A$theta < (loop[i] + (1 / (10^precision)) / 2)), ]
    }

    tempid <- which.max(tempA$r) # no need to account for ties
    if (length(tempid) != 0) {
      ext$X[i] <- tempA$X[tempid]
      ext$Y[i] <- tempA$Y[tempid]
      ext$theta[i] <- loop[i]
      ext$r[i] <- tempA$r[tempid]
    }
    rm(tempid, tempA)
  }

  # Convert to data frame
  ext <- na.omit(as.data.frame(ext))

  # If encasing is defined then increase the boundary
  if (!is.na(encase)) {
    ext$X <- ext$X + (encase * cos(ext$theta))
    ext$Y <- ext$Y + (encase * sin(ext$theta))
  }

  # Plot the boundary
  if (plotting) {
    graphics::points(ext$X, ext$Y, col = "green", pch = 16)
  }

  # Remove the no-longer need shifted coordinates
  ext <- ext[ , c("X", "Y", "theta", "r")]

  # If including smoothing then do so via loess()
  if (smoothed) {

    # to avoid discontinuity at the start of the cyclical smoothing, extend the fit across an extra repeated quarter of data on either side
    temp1 <- seq(loop[length(loop)] + (1/(10^precision)), loop[length(loop)] + round(length(loop) / 4, 0) * (1/(10^precision)), by = (1/(10^precision)))
    temp2 <- seq(loop[length(1)] - round(length(loop) / 4, 0) * (1/(10^precision)), loop[1] - (1/(10^precision)), by = (1/(10^precision)))

    # 'smoothed' will determine the span parameter within loess()
    temp <- stats::loess(c(ext$r[(length(ext$r) - (length(temp2) - 1)):length(ext$r)], ext$r, ext$r[1:length(temp1)]) ~ c(temp2, ext$theta, temp1), span = smoothed)$fitted
    ext$rsmooth <- temp[(c(temp2, ext$theta, temp1) <= max(ext$theta)) & (c(temp2, ext$theta, temp1) >= min(ext$theta))]
    ext$Xsmooth <- (ext$rsmooth * cos(ext$theta)) + mid["X"]
    ext$Ysmooth <- (ext$rsmooth * sin(ext$theta)) + mid["Y"]

    # Adjust for encasing
    if (!is.na(encase)) {
      ext$Xsmooth <- ext$Xsmooth + (encase * cos(ext$theta))
      ext$Ysmooth <- ext$Ysmooth + (encase * sin(ext$theta))
    }

    if (plotting){
      graphics::lines(ext$Xsmooth, ext$Ysmooth, col = "red")
    }
  }

  bnd.res <- ext[ , c("X", "Y")]
  colnames(bnd.res) <- coord.names

  return(bnd.res)
}
