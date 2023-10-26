#' Internal scampr function that uses TMB to make the objective function of scampr models
#'
#' @param TMB.inputs a list created by the internal scampr function: \code{get.TMB.data.input}.
#' @param maxit a numeric indicating the maximum number of iterations for the optimizer. Default is 100 as this optimization uses gradient information.
#'
#' @return a data.frame (sparse or dense depending on parameter bf.matrix.type)
#' @export
#'
#' @importFrom TMB MakeADFun
#'
#' @examples
#' # Get the flora data for one of the species
#' dat_po <- flora$po$sp1
#' dat_pa <- flora$pa
#'
#' # obtain a sample of 10,000 quadrature points for the point process model
#' set.seed(1)
#' quad.pts <- flora$quad[sample(1:nrow(flora$quad), 10000, replace = F), ]
#' set.seed(NULL)
#'
#' # Attach the quadrature points to the presence-only data
#' dat_po <- rbind.data.frame(dat_po, quad.pts)
#'
#' # Ensure the "response" variable in each data set shares the same name
#' dat_po$presence <- dat_po$pres
#' dat_pa$presence <- dat_pa$sp1
#'
#' # Get the TMB data lists for a combined data model without latent field
#' tmb.input <- scampr:::get.TMB.data.input(presence ~ MNT, bias.formula ~ D.Main, po.data = dat_po, pa.data = dat_pa)
#'
#' # create the objective function
#' obj <- make.objective.function(tmb.input)
make.objective.function <- function(TMB.inputs, maxit = 100) {
  if (TMB.inputs$args$random.bias.type == "none") {
    objective.fn <-switch(TMB.inputs$args$approx.type,
           not_sre = TMB::MakeADFun(data = TMB.inputs$tmb.data, parameters = TMB.inputs$tmb.pars, DLL = "scampr", map = list(random = factor(rep(NA, length(TMB.inputs$tmb.pars$random))), log_variance_component = factor(rep(NA, length(TMB.inputs$tmb.pars$log_variance_component)))), silent = T),
           variational = TMB::MakeADFun(data = TMB.inputs$tmb.data, parameters = TMB.inputs$tmb.pars, DLL = "scampr", silent = T),
           laplace = TMB::MakeADFun(data = TMB.inputs$tmb.data, parameters = TMB.inputs$tmb.pars, random = "random", DLL = "scampr", silent = T)
    )
  } else {
    objective.fn <-switch(TMB.inputs$args$approx.type,
                          not_sre = stop("model does not contain SRE however, random biasing effects were found"),
                          variational = stop("model found to contain random biasing effects however, approx. type is variational"),
                          laplace = TMB::MakeADFun(data = TMB.inputs$tmb.data, parameters = TMB.inputs$tmb.pars, random = c("random", "random_bias"), DLL = "scampr", silent = T)
    )
  }
  objective.fn$control = list(maxit = maxit)
  return(objective.fn)
}
