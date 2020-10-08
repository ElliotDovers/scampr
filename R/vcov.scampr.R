#' vcov for objects of class 'scampr'
#'
#' @param object A scampr model object
#'
#' @return
#' @export
#'
#' @examples
vcov.scampr <- function(object) {
  # get an identifier for the model type
  mod.id <- NULL
  if (is.na(object$approx.type)) {
    mod.id <- "ipp"
  } else {
    mod.id <- object$approx.type
  }
  des.mat <- scampr:::get.desgin.matrix(object$formula, object$data)
  pt.quad.id <- object$pt.quad.id
  fixed.names <- colnames(des.mat)
  # Check the type of basis functions used. Currently supports simple or FRK basis
  if (is.null(object$FRK.basis.functions)) {
    if (is.null(object$basis.fn.info)) {
      if (!is.na(object$approx.type)) {
        stop("The model's basis functions are incompatible with obtaining a Hessian matrix")
      } else {
        bf.matrix <- matrix(rep(0, 2*nrow(des.mat)), nrow = nrow(des.mat))
      }
    } else {
      bf.matrix <- scampr:::get.bf.matrix(object$basis.fn.info, object$data[ , object$coord.names])
    }
  } else {
    bf.matrix <- FRK::eval_basis(basis = object$FRK.basis.functions, as.matrix(object$data[ , object$coord.names]))
  }
  # TMB required data setup
  dat.list <- list(
    X_pres = des.mat[pt.quad.id == 1, ],
    Z_pres = as(bf.matrix[pt.quad.id == 1, ], "sparseMatrix"),
    X_quad = des.mat[pt.quad.id == 0, ],
    Z_quad = as(bf.matrix[pt.quad.id == 0, ], "sparseMatrix"),
    quad_size = object$data[ , object$quad.weights.name][pt.quad.id == 0],
    bf_per_res = switch(mod.id,
                        ipp = 2,
                        variational = object$basis.per.res,
                        laplace = object$basis.per.res),
    mod_type = as.integer(which(mod.id == c("ipp", "variational", "laplace")) - 1)
  )
  # initilise starting parameters at existing ML estimates to speed things up
  start.pars <- lapply(split(object$par, names(object$par)), unname)
  if (mod.id == "ipp") {
    start.pars$random <- rep(0, ncol(dat.list$Z_pres))
    start.pars$log_variance_component <- rep(0, ncol(dat.list$Z_pres))
  } else if (mod.id == "laplace") {
    start.pars$random <- unname(object$random.effects[grepl(" Mean ", row.names(object$random.effects), fixed = T), 1L])
  }
  # set up the objective function w.r.t. model identifier
  obj <- switch(mod.id,
                ipp = TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", map = list(random = factor(rep(NA, ncol(dat.list$Z_pres))), log_variance_component = factor(rep(NA, ncol(dat.list$Z_pres)))), silent = T),
                variational = TMB::MakeADFun(data = dat.list, parameters = start.pars, DLL = "scampr", silent = T),
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
  return(vcov_mat)
}
