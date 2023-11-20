# scampr
Spatially Correlated, Approximate Modelling of Point patterns in R

scampr is an R package that offers a regression-style framework for modelling spatial point patterns (referred to as presence-only data in ecology)
using a log-Gaussian Cox Process (LGCP). The package can also be used to fit a binomial model with spatial random effects to binary data
(*e.g.* presence/absence data in ecology), as well as an integrated model that combines the two while permitting shared spatially correlated latent field(s).

Unlike the inhomogeneous Poisson process, LGCPs offer a way to incorporate additional spatial clustering into models by including a Gaussian random field (GRF)
to induce additional spatial correlation between observations — effectively acting as a spatially correlated error term. LGCP models are particularly appropriate
in instances where clustering arises from missing or unmeasured environmental predictors/phenomena, as opposed to those in which clustering/dispersal is due to
interactions between point events.

Fitting LGCP models can be difficult and time consuming and, as a result, limits the ability of researchers to flexibly analyse spatial point pattern data.
scampr is a fast, approximate maximum-likelihood approach to fitting LGCP to 2D spatial point patterns, involving a combination of three innovations. First, variational
approximation (VA) permits a closed form approximation to the marginalised log-likelihood. Second, fixed rank kriging provides a rank reduced approximation to the
large spatial variance-covariance matrices that arise and are otherwise very computationally demanding. Finally, automatic differentiation is used to quickly obtain
gradient information for efficient optimization and inference.

Models are fit using `scampr()` that adopts a simple interface/syntax following the common regression modelling formats used on R, *e.g.* `lm()`, `glm()`. Many of the common S3
functions familiar to users (such as `summary`, `plot`, `simulate`, `predict`, `logLik`, `AIC`, `confint` and `vcov`) are also available to scampr models. The package is built upon the
advances of TMB (Kristensen et al., 2016) which enables coding the likelihoods in C++, as well as providing automatic differentiation for easy access to gradient
information — permitting fast optimisation, automated Laplace approximation, and automated estimation of the variance-covariance matrix of parameter estimates in scampr
models.

The package name stands for Spatially Correlated, Approximate Modelling of Point patterns in R however, the verb "scamper" — to run with quick, light steps — perfectly
captures the motivation of this package: to give researchers access to complex spatial models that fit quickly and require only a light touch.

## Installing scampr

As scampr is not on CRAN please install via:

`devtools::install_github("ElliotDovers/scampr", dependencies = "Imports", upgrade = "never")`

## Publications
- [Dovers et al. (2023) "Fast, Approximate Maximum Likelihood Estimation of Log-Gaussian Cox Processes", *Journal of Computational and Graphical Statistics*](https://www.tandfonline.com/doi/full/10.1080/10618600.2023.2182311)
- [Dovers et al. (2023) "A fast method for fitting integrated species distribution models", *Methods in Ecology and Evolution*](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.14252)
