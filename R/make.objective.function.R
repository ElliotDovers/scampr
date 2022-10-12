#' Internal scampr function that uses TMB to make the objective function of scampr models
#'
#' @param TMB.inputs a list created by the internal scampr function: \code{get.TMB.data.input}.
#' @param maxit a numeric indicating the maximum number of iterations for the optimizer. Default is 1000.
#'
#' @return a data.frame (sparse or dense depending on parameter bf.matrix.type)
#' @export
#'
#' @importFrom TMB MakeADFun
#'
#' @examples
#' # Get the gorilla nesting data
#' dat <- gorillas
#'
make.objective.function <- function(TMB.inputs, maxit = 1000) {
  objective.fn <- switch(TMB.inputs$args$bias.type,
                none = switch(TMB.inputs$args$model.type,
                              ipp = TMB::MakeADFun(data = TMB.inputs$tmb.data, parameters = TMB.inputs$tmb.pars, DLL = "scampr", map = list(random = factor(rep(NA, length(TMB.inputs$tmb.pars$random))), log_variance_component = factor(rep(NA, length(TMB.inputs$tmb.pars$log_variance_component)))), silent = T),
                              variational = TMB::MakeADFun(data = TMB.inputs$tmb.data, parameters = TMB.inputs$tmb.pars, DLL = "scampr", silent = T),
                              laplace = TMB::MakeADFun(data = TMB.inputs$tmb.data, parameters = TMB.inputs$tmb.pars, random = "random", DLL = "scampr", silent = T)
                ),
                covariates = switch(TMB.inputs$args$model.type,
                                    ipp = TMB::MakeADFun(data = TMB.inputs$tmb.data, parameters = TMB.inputs$tmb.pars, DLL = "scampr", map = list(random = factor(rep(NA, length(TMB.inputs$tmb.pars$random))), log_variance_component = factor(rep(NA, length(TMB.inputs$tmb.pars$log_variance_component)))), silent = T),
                                    variational = TMB::MakeADFun(data = TMB.inputs$tmb.data, parameters = TMB.inputs$tmb.pars, DLL = "scampr", silent = T),
                                    laplace = TMB::MakeADFun(data = TMB.inputs$tmb.data, parameters = TMB.inputs$tmb.pars, random = "random", DLL = "scampr", silent = T)
                ),
                latent = switch(TMB.inputs$args$model.type,
                                ipp = TMB::MakeADFun(data = TMB.inputs$tmb.data, parameters = TMB.inputs$tmb.pars, DLL = "scampr", map = list(random = factor(rep(NA, length(TMB.inputs$tmb.pars$random))), log_variance_component = factor(rep(NA, length(TMB.inputs$tmb.pars$log_variance_component)))), silent = T),
                                variational = TMB::MakeADFun(data = TMB.inputs$tmb.data, parameters = TMB.inputs$tmb.pars, random = c("random", "bias"), DLL = "scampr", silent = T),
                                laplace = TMB::MakeADFun(data = TMB.inputs$tmb.data, parameters = TMB.inputs$tmb.pars, random = c("random", "bias"), DLL = "scampr", silent = T)
                ),
                new_latent = switch(TMB.inputs$args$model.type,
                                    ipp = TMB::MakeADFun(data = TMB.inputs$tmb.data, parameters = TMB.inputs$tmb.pars, DLL = "scampr", map = list(random = factor(rep(NA, length(TMB.inputs$tmb.pars$random))), log_variance_component = factor(rep(NA, length(TMB.inputs$tmb.pars$log_variance_component)))), silent = T),
                                    variational = TMB::MakeADFun(data = TMB.inputs$tmb.data, parameters = TMB.inputs$tmb.pars, random = c("random", "bias"), DLL = "scampr", silent = T),
                                    laplace = TMB::MakeADFun(data = TMB.inputs$tmb.data, parameters = TMB.inputs$tmb.pars, random = c("random", "bias"), DLL = "scampr", silent = T)
                )
  )
  objective.fn$control = list(maxit = maxit)
  return(objective.fn)
}
