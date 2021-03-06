% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcov.scampr.R
\name{vcov.scampr}
\alias{vcov.scampr}
\title{vcov for objects of class 'scampr'}
\usage{
\method{vcov}{scampr}(object, ...)
}
\arguments{
\item{object}{A scampr model object}

\item{...}{NA}
}
\value{
A matrix of the estimated covariances between the parameter estimates in the linear or non-linear predictor of the model. This should have row and column names corresponding to the parameter names given by the coef method.
}
\description{
vcov for objects of class 'scampr'
}
\examples{
# Get the gorilla nesting data
dat <- gorillas

# Standardise the elevation covariate
dat$elev.std <- scale(dat$elevation)

# Fit a scampr model to the point pattern
m <- scampr(pres ~ elev.std, data = dat, model.type = "ipp")

vcov(m)
}
