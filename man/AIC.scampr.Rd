% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AIC.scampr.R
\name{AIC.scampr}
\alias{AIC.scampr}
\title{Akaike's Information Criteria for objects of class 'scampr'. Note that this currently calculates mAIC}
\usage{
\method{AIC}{scampr}(object, ..., k = 2)
}
\arguments{
\item{object}{a scampr model}

\item{...}{Optionally, additional scampr model objects}

\item{k}{a numeric describing the penalty per parameter - defaults to k = 2 i.e. classical AIC.}
}
\value{
a numeric value with the corresponding AIC (or BIC, or ..., depending on k).
}
\description{
Akaike's Information Criteria for objects of class 'scampr'. Note that this currently calculates mAIC
}
\examples{
# Get the gorilla nesting data
dat <- gorillas

# Standardise the elevation covariate
dat$elev.std <- scale(dat$elevation)

# Fit a scampr model to the point pattern
m <- scampr(pres ~ elev.std, data = dat, model.type = "PO")

AIC(m)
}
