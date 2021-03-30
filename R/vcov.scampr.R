#' vcov for objects of class 'scampr'
#'
#' @param object A scampr model object
#' @param ... NA
#'
#' @return A matrix of the estimated covariances between the parameter estimates in the linear or non-linear predictor of the model. This should have row and column names corresponding to the parameter names given by the coef method.
#' @export
#'
#' @importFrom TMB MakeADFun sdreport
#'
#' @examples
#' # Get the gorilla nesting data
#' dat <- gorillas
#'
#' # Standardise the elevation covariate
#' dat$elev.std <- scale(dat$elevation)
#'
#' # Fit a scampr model to the point pattern
#' m <- ipp(pres ~ elev.std, data = dat)
#'
#' vcov(m)
vcov.scampr <- function(object, ...) {

  # get an identifier for the model type
  mod.id <- NULL
  if (is.na(object$approx.type)) {
    mod.id <- "ipp"
  } else {
    mod.id <- object$approx.type
  }

  # Need separate cases for each data type #

  if (object$data.model.type == "po") {

    # Check the type of basis functions used. Currently supports simple or FRK basis
    if (is.null(object$FRK.basis.functions)) {
      if (is.null(object$basis.fn.info)) {
        if (!is.na(object$approx.type)) {
          stop("The model's basis functions are incompatible with obtaining a Hessian matrix")
        } else {
          inputs <- get.TMB.data.input(formula = object$formula, data = object$data,
                                                coord.names = object$coord.names,
                                                quad.weights.name = object$quad.weights.name,
                                                model.type = mod.id,
                                                data.type = object$data.model.type, starting.pars = object)
        }
      } else {
        inputs <- get.TMB.data.input(formula = object$formula, data = object$data,
                                              coord.names = object$coord.names,
                                              quad.weights.name = object$quad.weights.name,
                                              simple.basis = object$basis.fn.info,
                                              model.type = mod.id, bf.matrix.type = object$bf.matrix.type,
                                              data.type = object$data.model.type, starting.pars = object)
      }
    } else {
      inputs <- get.TMB.data.input(formula = object$formula, data = object$data,
                                  coord.names = object$coord.names,
                                  quad.weights.name = object$quad.weights.name,
                                  FRK.basis.functions = object$FRK.basis.functions,
                                  model.type = mod.id, bf.matrix.type = object$bf.matrix.type,
                                  data.type = object$data.model.type, starting.pars = object)
    }
    # TMB required data setup
    dat.list <- inputs$tmb.data
    # initilise starting parameters at existing ML estimates to speed things up
    start.pars <- inputs$tmb.pars
    # set up the objective function w.r.t. mod.id
    obj <- switch(mod.id,
                  ipp = TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", map = list(bias = as.factor(rep(NA, ncol(dat.list$B_PO_pres))), random = factor(rep(NA, ncol(dat.list$Z_PA))), log_variance_component = factor(rep(NA, length(dat.list$bf_per_res)))), silent = T),
                  variational = TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", map = list(bias = as.factor(rep(NA, ncol(dat.list$B_PO_pres)))), silent = T),
                  laplace = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", map = list(bias = as.factor(rep(NA, ncol(dat.list$B_PO_pres)))), silent = T)
    )
    # get the sdreport that produces the full precision matrix
    tmp <- TMB::sdreport(obj, getJointPrecision = T) # This will change pre to post optim but ok with ML starting pars
    if (is.null(tmp$jointPrecision)) {
      vcov_mat <- tmp$cov.fixed
    } else {
      Hess <- as.matrix(tmp$jointPrecision)
      vcov_mat <- solve(Hess)
    }
    # re-name the rows and columns
    if (mod.id == "laplace") {
      rownames(vcov_mat)[rownames(vcov_mat) != "random"] <- names(object$coefficients)
      colnames(vcov_mat)[colnames(vcov_mat) != "random"] <- names(object$coefficients)
    } else {
      rownames(vcov_mat) <- names(object$coefficients)
      colnames(vcov_mat) <- names(object$coefficients)
    }
  } else if (object$data.model.type == "pa") {

    # Check the type of basis functions used. Currently supports simple or FRK basis
    if (is.null(object$FRK.basis.functions)) {
      if (is.null(object$basis.fn.info)) {
        if (!is.na(object$approx.type)) {
          stop("The model's basis functions are incompatible with obtaining a Hessian matrix")
        } else {
          inputs <- get.TMB.data.input(pa.formula = object$formula, pa.data = object$data,
                                                coord.names = object$coord.names,
                                                quad.weights.name = object$quad.weights.name,
                                                model.type = mod.id,
                                                data.type = object$data.model.type, starting.pars = object)
        }
      } else {
        inputs <- get.TMB.data.input(pa.formula = object$formula, pa.data = object$data,
                                              coord.names = object$coord.names,
                                              quad.weights.name = object$quad.weights.name,
                                              simple.basis = object$basis.fn.info,
                                              model.type = mod.id, bf.matrix.type = object$bf.matrix.type,
                                              data.type = object$data.model.type, starting.pars = object)
      }
    } else {
      inputs <- get.TMB.data.input(pa.formula = object$formula, pa.data = object$data,
                                            coord.names = object$coord.names,
                                            quad.weights.name = object$quad.weights.name,
                                            FRK.basis.functions = object$FRK.basis.functions,
                                            model.type = mod.id, bf.matrix.type = object$bf.matrix.type,
                                            data.type = object$data.model.type, starting.pars = object)
    }
    # TMB required data setup
    dat.list <- inputs$tmb.data
    # initilise starting parameters at existing ML estimates to speed things up
    start.pars <- inputs$tmb.pars
    # # AT THIS STAGE CAN ONLY PERFORM LAPLAC APPROX.
    # set up the objective function w.r.t. mod.id
    obj <- switch(mod.id,
                  ipp = TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", map = list(bias = as.factor(rep(NA, ncol(dat.list$B_PO_pres))), random = factor(rep(NA, ncol(dat.list$Z_PA))), log_variance_component = factor(rep(NA, length(dat.list$bf_per_res)))), silent = T),
                  variational = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", map = list(bias = as.factor(rep(NA, ncol(dat.list$B_PO_pres)))), silent = T),
                  laplace = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", map = list(bias = as.factor(rep(NA, ncol(dat.list$B_PO_pres)))), silent = T)
    )
    # get the sdreport that produces the full precision matrix
    tmp <- TMB::sdreport(obj, getJointPrecision = T) # This will change pre to post optim but ok with ML starting pars
    if (is.null(tmp$jointPrecision)) {
      vcov_mat <- tmp$cov.fixed
    } else {
      Hess <- as.matrix(tmp$jointPrecision)
      vcov_mat <- solve(Hess)
    }
    # re-name the rows and columns
    if (mod.id == "laplace") {
      rownames(vcov_mat)[rownames(vcov_mat) != "random"] <- names(object$coefficients)
      colnames(vcov_mat)[colnames(vcov_mat) != "random"] <- names(object$coefficients)
    } else {
      rownames(vcov_mat) <- names(object$coefficients)
      colnames(vcov_mat) <- names(object$coefficients)
    }

  } else if (object$data.model.type == "popa") {

    data.po <- object$data
    data.pa <- attr(object$data, "pa")
    form.po <- object$formula
    form.pa <- attr(object$formula, "pa")

    # Check the type of basis functions used. Currently supports simple or FRK basis
    if (is.null(object$FRK.basis.functions)) {
      if (is.null(object$basis.fn.info)) {
        if (!is.na(object$approx.type)) {
          stop("The model's basis functions are incompatible with obtaining a Hessian matrix")
        } else {
          inputs <- get.TMB.data.input(formula = form.po, pa.formula = form.pa, data = data.po,
                                                pa.data = data.pa,
                                                coord.names = object$coord.names,
                                                quad.weights.name = object$quad.weights.name,
                                                model.type = mod.id,
                                                data.type = object$data.model.type, starting.pars = object)
        }
      } else {
        inputs <- get.TMB.data.input(formula = form.po, pa.formula = form.pa, data = data.po,
                                              pa.data = data.pa,
                                              coord.names = object$coord.names,
                                              quad.weights.name = object$quad.weights.name,
                                              simple.basis = object$basis.fn.info,
                                              model.type = mod.id, bf.matrix.type = object$bf.matrix.type,
                                              data.type = object$data.model.type, starting.pars = object)
      }
    } else {
      inputs <- get.TMB.data.input(formula = form.po, pa.formula = form.pa, data = data.po,
                                            pa.data = data.pa,
                                            coord.names = object$coord.names,
                                            quad.weights.name = object$quad.weights.name,
                                            FRK.basis.functions = object$FRK.basis.functions,
                                            model.type = mod.id, bf.matrix.type = object$bf.matrix.type,
                                            data.type = object$data.model.type, starting.pars = object)
    }
    # TMB required data setup
    dat.list <- inputs$tmb.data
    # initilise starting parameters at existing ML estimates to speed things up
    start.pars <- inputs$tmb.pars
    # # AT THIS STAGE CAN ONLY PERFORM LAPLAC APPROX.
    # set up the objective function w.r.t. mod.id
    obj <- switch(mod.id,
                  ipp = TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", map = list(random = factor(rep(NA, ncol(dat.list$Z_PA))), log_variance_component = factor(rep(NA, length(dat.list$bf_per_res)))), silent = T),
                  variational = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", silent = T),
                  laplace = TMB::MakeADFun(data = dat.list, parameters = start.pars, random = "random", DLL = "scampr", silent = T)
    )
    # get the sdreport that produces the full precision matrix
    tmp <- TMB::sdreport(obj, getJointPrecision = T) # This will change pre to post optim but ok with ML starting pars
    if (is.null(tmp$jointPrecision)) {
      vcov_mat <- tmp$cov.fixed
    } else {
      Hess <- as.matrix(tmp$jointPrecision)
      vcov_mat <- solve(Hess)
    }
    # re-name the rows and columns
    if (mod.id == "laplace") {
      rownames(vcov_mat)[rownames(vcov_mat) != "random"] <- names(object$coefficients)
      colnames(vcov_mat)[colnames(vcov_mat) != "random"] <- names(object$coefficients)
    } else {
      rownames(vcov_mat) <- names(object$coefficients)
      colnames(vcov_mat) <- names(object$coefficients)
    }

  } else {
    stop("Incorrect 'data.type' provided. Must be one of 'po', 'pa' or 'popa'")
  }

  return(vcov_mat)
}
