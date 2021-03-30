#' Spatially correlated, approximate modelling of presences in R
#'
#' @description This is the main function for modelling presences within the \code{scampr} framework. The type of model will depend on the arguments provided. This can be used to model point patterns as a log-Gaussian Cox process (LGCP) or Inhomogeneous Poisson process (IPP), as well as, jointly fitting a model to presence-only and presence/absence data if both are provided. This function can also fit binary regression with a complimentary log-log link function (with optional spatial random effects) when only presence/absence data is provided.
#'
#' If only \code{formula} and \code{data} are provided, the function will fit either an IPP, or LGCP model to the point pattern (depending on argument \code{model.type}). This uses numerical quadrature (provided with the data, see e.g. scampr::gorillas) to approximate the spatial integral. If fitting a LGCP, uses one of either a Laplace or variational approximation to marginalise over the latent field.
#'
#' If only \code{pa.formula} and \code{pa.data} are provided, the function will fits a binary regression model to presence/absence data using a complimentary log-log link function. Can accomodate an approximate latent field as spatial random effects (depending on argument \code{model.type}).
#'
#' If both \code{formula} and \code{pa.formula}, as well as, \code{data} and \code{pa.data} are provided, the function jointly fits a model to presence-only and presence/absence data as linked by response to environmental predictors provided in each formula. The presence-only formula must also contain biasing predictors to account for opportunistic collection. If argument \code{model.type} is not equal to "ipp", then the two data sources will additionally share a latent Gaussian random field.
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the Presence-only data model to be fitted. The 'response' must be a binary that indicates whether a datum is a presence or quadrature point. See GLM function for further formula details.
#' @param data a data frame containing predictors at both point locations and quadrature as well as the formula 'response' found in \code{formula}.
#' @param pa.formula Optionally, an object of class "formula" (or one that can be coerced to that class): a symbolic description of the Presence/Absence data model to be fitted. All predictor terms must also be included in the \code{formula} argument. The 'response' must be a must be either the site abundance or a binary indicating whether there is a presence or not. See GLM function for further formula details.
#' @param pa.data Optionally, a data frame containing predictors and response described in \code{pa.formula}.
#' @param coord.names a vector of length 2 containing character strings describing the column names of the coordinates in any data provided. First coordinate name should refer to the horizontal axis.
#' @param quad.weights.name a charater string of the column name of quadrature weights in \code{data}.
#' @param FRK.basis.functions an optional object of class "Basis" from FRK package. If neither \code{FRK.basis.functions} nor \code{simple.basis} are specified will use default \code{FRK::auto_basis()} with a single spatial resolutions - if the model type is not set to an IPP.
#' @param simple.basis an alternative to \code{FRK.basis.functions}: a data frame of basis functions information created by '\code{scampr::simple_basis()}.
#' @param model.type a character string indicating the type of model to be used. May be one of 'laplace' or 'variational' for Cox Processes involving spatially correlated errors or 'ipp' for a model that follows an inhomgeneous Poisson process.
#' @param bf.matrix.type a character string, one of 'sparse' or 'dense' indicating whether to use sparse or dense matrix computations for the basis functions created.
#' @param se a logical indicating whether standard errors should be calculated.
#' @param starting.pars an optional named list or previously fit scampr model object that gives warm starting values for the parameters of the model.
#' @param subset an optional vector describing a subset of the data to be used.
#'
#' @return a scampr model object
#' @export
#'
#' @importFrom methods as
#' @importFrom stats as.formula optim
#' @importFrom sp coordinates
#' @importFrom FRK auto_basis eval_basis
#' @importFrom TMB MakeADFun sdreport
#'
#' @examples
#' # Get the Eucalypt data
#' dat_po <- eucalypt[["po"]]
#' dat_pa <- eucalypt[["pa"]]
#'
#' # Fit models without a latent effects (IPP) #
#' # Point Process Model
#' m.ipp <- scampr(pres ~ TMP_MIN + D_MAIN_RDS, dat_po, model.type = "ipp")
#' # Binary Regression
#' m.bin <- scampr(pa.formula = Y ~ TMP_MIN, pa.data = dat_pa, model.type = "ipp")
#' # Combined Data Model
#' m.comb <- scampr(pres ~ TMP_MIN + D_MAIN_RDS, dat_po, Y ~ TMP_MIN,
#' dat_pa, model.type = "ipp")
#'
#' # Set up a simple 2D grid of basis functions to fit a LGCP model to the data
#' bfs <- simple_basis(nodes.on.long.edge = 9, data = dat_po)
#'
#' \dontrun{
#' # Fit with a shared latent field (LGCP) #
#' # Point Process Model
#' m.lgcp <- scampr(pres ~ TMP_MIN + D_MAIN_RDS, dat_po)
#' # Binary Regression with spatial random effects
#' m.bin_w_sre <- scampr(pa.formula = Y ~ TMP_MIN, pa.data = dat_pa)
#' # Combined Data Model with spatial random effects
#' m.comb_w_sre <- scampr(pres ~ TMP_MIN + D_MAIN_RDS, dat_po, Y ~ TMP_MIN,
#' dat_pa)
#' }
scampr <- function(formula, data, pa.formula, pa.data, coord.names = c("x", "y"), quad.weights.name = "quad.size", FRK.basis.functions, simple.basis, model.type = c("variational", "laplace", "ipp"), bf.matrix.type = c("sparse", "dense"), se = TRUE, starting.pars, subset) {

  # Determine the data model type to be used based off provided arguments
  data.type <- NULL
  if((missing(pa.formula) & missing(pa.data)) & !missing(formula) & !missing(data)) {
    data.type <- "po"
  } else if (!missing(pa.formula) & !missing(pa.data) & missing(formula) & missing(data)) {
    data.type <- "pa"
  } else if (!missing(pa.formula) & !missing(pa.data) & !missing(formula) & !missing(data)) {
    data.type <- "popa"
  } else {
    stop("Please supply either 'formula' and 'data' or 'pa.formula' and 'pa.data' or all of the above for a combined data model!")
  }

  # Obtain the arguments to be passed to the appropriate function
  argpass <- as.list(match.call())
  argpass[[1]] <- NULL
  names(argpass)[names(argpass) %in% "formula"] <- "po.formula"
  names(argpass)[names(argpass) %in% "data"] <- "po.data"

  # Try using eval(parse( to correct the scoping issues involved with do.call()
  # tmp <- paste(paste(names(argpass), unname(argpass), sep = " = "), collapse = ", ")
  # THIS HAS PROBLEMS WITH CHARACTER ARGUMENTS! JUST GOING TO DITCH TESTING ON THIS

  if (data.type == "po") { # Model the presence-only type
    #######################
    #### PO data model ####
    #######################
    cpu <- system.time(assign("mod", do.call(po, argpass)))
    mod$cpu <- cpu
    # fn_as_char <- paste0("po(", tmp, ")")

  } else if (data.type == "pa") {
    #######################
    #### PA data model ####
    #######################
    cpu <- system.time(assign("mod", do.call(pa, argpass)))
    mod$cpu <- cpu
    # fn_as_char <- paste0("pa(", tmp, ")")

  } else if (data.type == "popa") {
    #########################
    #### POPA data model ####
    #########################
    cpu <- system.time(assign("mod", do.call(popa, argpass)))
    mod$cpu <- cpu
    # fn_as_char <- paste0("popa(", tmp, ")")

  }

  # cpu <- system.time(assign("mod", eval(parse(text = fn_as_char))))
  # mod$cpu <- cpu

  return(mod)

}
