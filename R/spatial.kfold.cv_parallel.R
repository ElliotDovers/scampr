#' Performs a spatial hold-one-out K fold cross validation on a scampr model (using parallelisation via \code{foreach}) #WIP#
#'
#' @param object a scampr model object
#' @param po.fold.id An integer or factor vector the same length as the presence-only data used in the model \code{object} - describes the CV fold that each location falls into.
#' @param pa.fold.id An integer or factor vector the same length as the presence/absence data used in the model \code{object} - describes the CV fold that each location falls into.
#' @param trunc.pa.prob A small positive number by which the predicted probability of presence is truncated. This can be used to ensure infinite values are avoided within the cross-validation.
#'
#' @return a list including- 'predicted.cll.pa': the predicted logLikelihood, conditional on the model predictors and latent field, on the held out presence/absence data. 'predicted.cll.po': the predicted logLikelihood, conditional on the model predictors and latent field, on the held out presence-only data. 'auc': the area under the ROC curve on the held out presence/absence sites. 'model': the model upon which the cross validation was performed.
#' @noRd
#'
#' @importFrom pROC roc auc
#' @importFrom foreach foreach "%dopar%"
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#'
#' @examples
#' #' # Get the Eucalypt data
#' dat_po <- eucalypt[["po"]]
#' dat_pa <- eucalypt[["pa"]]
#'
#' \dontrun{
#' # Fit with a shared latent field
#' res <- simple_basis_search_popa(pres ~ TMP_MIN + D_MAIN_RDS, Y ~ TMP_MIN,
#' po.data = dat_po, pa.data = dat_pa)
#' }
spatial.kfold.cv_parallel <- function(object, po.fold.id, pa.fold.id, trunc.pa.prob = 1e-7) {

  if (object$model.type == "IDM") {

    # check the appropriate spatial folds are present
    if (missing(po.fold.id) & missing(pa.fold.id)) {
      stop(paste0("'po.fold.id' and 'pa.fold.id' must be supplied for a model of type = ", object$model.type))
    }

    # Convert CV folds to numerics in case factor is provided
    po.fold.id <- as.numeric(po.fold.id)
    pa.fold.id <- as.numeric(pa.fold.id)

    # Extract the data (these elements need to be evaluated here to be accessible to parallel cores)
    po.data <- object$data
    pa.data <- attr(object$data, "PA")
    po.resp <- po.data[, all.vars(object$formula[[2]])]
    pa.resp <- pa.data[, all.vars(object$formula[[2]])]
    quad.sizes <- po.data[, object$quad.weights.name]
    bfs_shared <- object$basis.functions
    bfs_bias <- object$po.biasing.basis.functions
    folds <- sort(unique(po.fold.id)) # fold levels
    po.row.id <- 1:nrow(po.data) # index reference for returning original PO dataset order
    pa.row.id <- 1:nrow(pa.data) # index reference for returning original PA dataset order

    # Determine the number of cores and set number of instances to initialise in parallel
    ncores <- parallel::detectCores() # partition number will be based on the cores available
    if (length(folds) < ncores) {
      k <- length(folds)
    } else {
      k <- ncores - 1
    }

    # Start the clusters
    sock <- parallel::makeCluster(rep("localhost", k), type = "SOCK")
    doParallel::registerDoParallel(sock)
    # print("Computing CV folds in parallel.", quote=FALSE)
    comb_res <- foreach::foreach(fold = folds, .combine = rbind, .packages='scampr') %dopar% {
      call.list <- as.list(object$call)
      # adjust the training data
      call.list$data <- po.data[po.fold.id != fold, ]
      call.list$IDM.presence.absence.df <- pa.data[pa.fold.id != fold, ]
      call.list$basis.functions <- bfs_shared
      call.list$po.biasing.basis.functions <- bfs_bias
      # train the model
      call.list[[1]] <- NULL
      train_mod <- do.call("scampr", call.list)
      # test on held out fold
      test.pa <- predict.scampr(train_mod, newdata = pa.data[pa.fold.id == fold, ])
      test.po <- predict.scampr(train_mod, newdata = po.data[po.fold.id == fold, ], include.bias.accounting = TRUE)
      # calculate the row indices of the new test datasets
      po.test.rows <- po.row.id[po.fold.id == fold]
      pa.test.rows <- pa.row.id[pa.fold.id == fold]
      # set the return object
      ret.obj <- data.frame(preds = c(test.pa, test.po),
                       row_id = c(pa.test.rows, po.test.rows),
                       # fold_id = rep(fold, length(c(pa.test.rows, po.test.rows))),
                       data_id = c(rep(0, length(pa.test.rows)), rep(1, length(po.test.rows)))
      )
      return(ret.obj)
    }
    parallel::stopCluster(sock)
    # print("Stopping parallel clusters.", quote=FALSE)

    # separate the results by data sets
    res_pa <- data.frame(comb_res[comb_res$data_id == 0, ])
    res_po <- data.frame(comb_res[comb_res$data_id == 1, ])

    # re-order the test datasets to match the raw datasets
    pred.po <- res_po$preds[order(res_po$row_id)]
    pred.pa <- res_pa$preds[order(res_pa$row_id)]

    # predicted presence probability
    pres_prob <- 1 - exp(-exp(pred.pa))
    # NEED TO TRUNCATE PRESENCE PROBABILITY BY 'trunc.pa.prob'
    pres_prob[pres_prob <= trunc.pa.prob] <- trunc.pa.prob
    pres_prob[1 - pres_prob <= trunc.pa.prob] <- 1 - trunc.pa.prob
    # predicted conditional log-likelihood on the PA data
    pll_pa <- sum((pa.resp * log(pres_prob)) + ((1 - pa.resp) * log(1 - pres_prob)))
    # predicted conditional log-likelihood on the PO data
    pll_po <- sum(pred.po[po.resp == 1]) - sum(quad.sizes[po.resp == 0] * exp(pred.po[po.resp == 0]))
    # area under ROC curve
    roccurve <- pROC::roc(pa.resp, as.vector(pres_prob), quiet = T)
    auc_val <- as.numeric(pROC::auc(roccurve))

  } else  if (object$model.type == "PO") {

    # check the appropriate spatial fold is present
    if (missing(po.fold.id)) {
      stop(paste0("'po.fold.id' must be supplied for a model of type = ", object$model.type))
    }

    # Convert CV fold into a numeric in case factor is provided
    po.fold.id <- as.numeric(po.fold.id)

    # Extract the data (these elements need to be evaluated here to be accessible to parallel cores)
    po.data <- object$data
    po.resp <- po.data[, all.vars(object$formula[[2]])]
    quad.sizes <- po.data[, object$quad.weights.name]
    bfs <- object$basis.functions
    folds <- sort(unique(po.fold.id)) # fold levels
    po.row.id <- 1:nrow(po.data) # index reference for returning original PO dataset order

    # Determine the number of cores and set number of instances to initialise in parallel
    ncores <- parallel::detectCores() # partition number will be based on the cores available
    if (length(folds) < ncores) {
      k <- length(folds)
    } else {
      k <- ncores - 1
    }

    # Start the clusters
    sock <- parallel::makeCluster(rep("localhost", k), type = "SOCK")
    doParallel::registerDoParallel(sock)
    print("Computing CV folds in parallel.", quote=FALSE)
    comb_res <- foreach::foreach(fold = folds, .combine = rbind, .packages='scampr') %dopar% {
      call.list <- as.list(object$call)
      # adjust the training data
      call.list$data <- po.data[po.fold.id != fold, ]
      call.list$basis.functions <- bfs
      # train the model
      call.list[[1]] <- NULL
      train_mod <- do.call("scampr", call.list)
      # test on held out fold
      test.po <- predict.scampr(train_mod, newdata = po.data[po.fold.id == fold, ], include.bias.accounting = TRUE)
      # calculate the row indices of the new test datasets
      po.test.rows <- po.row.id[po.fold.id == fold]
      # set the return object
      ret.obj <- data.frame(preds = test.po,
                       row_id = po.test.rows
      )
      return(ret.obj)
    }
    parallel::stopCluster(sock)
    print("Stopping parallel clusters.", quote=FALSE)

    # re-order the test datasets to match the raw dataset
    pred.po <- comb_res$preds[order(comb_res$row_id)]

    # predicted conditional log-likelihood on the PA data
    pll_pa <- NA
    # predicted conditional log-likelihood on the PO data
    pll_po <- sum(pred.po[po.resp == 1]) - sum(quad.sizes[po.resp == 0] * exp(pred.po[po.resp == 0]))
    # area under ROC curve
    auc_val <- NA

  } else if (object$model.type == "PA") {

    # check the appropriate spatial folds are present
    if (missing(pa.fold.id)) {
      stop(paste0("'pa.fold.id' must be supplied for a model of type = ", object$model.type))
    }

    # Convert CV fold to a numeric in case factor is provided
    pa.fold.id <- as.numeric(pa.fold.id)

    # Extract the data (these elements need to be evaluated here to be accessible to parallel cores)
    pa.data <- object$data
    pa.resp <- pa.data[, all.vars(object$formula[[2]])]
    bfs <- object$basis.functions
    folds <- sort(unique(pa.fold.id)) # fold levels
    pa.row.id <- 1:nrow(pa.data) # index reference for returning original PA dataset order

    # Determine the number of cores and set number of instances to initialise in parallel
    ncores <- parallel::detectCores() # partition number will be based on the cores available
    if (length(folds) < ncores) {
      k <- length(folds)
    } else {
      k <- ncores - 1
    }

    # Start the clusters
    sock <- parallel::makeCluster(rep("localhost", k), type = "SOCK")
    doParallel::registerDoParallel(sock)
    print("Computing CV folds in parallel.", quote=FALSE)
    comb_res <- foreach::foreach(fold = folds, .combine = rbind, .packages='scampr') %dopar% {
      call.list <- as.list(object$call)
      # adjust the training data
      call.list$data <- pa.data[pa.fold.id != fold, ]
      call.list$basis.functions <- bfs
      # train the model
      call.list[[1]] <- NULL
      train_mod <- do.call("scampr", call.list)
      # test on held out fold
      test.pa <- predict.scampr(train_mod, newdata = pa.data[pa.fold.id == fold, ])
      # calculate the row indices of the new test dataset
      pa.test.rows <- pa.row.id[pa.fold.id == fold]
      # set the return object
      ret.obj <- data.frame(preds = test.pa,
                            row_id = pa.test.rows
      )
      return(ret.obj)
    }
    parallel::stopCluster(sock)
    print("Stopping parallel clusters.", quote=FALSE)

    # re-order the test datasets to match the raw dataset
    pred.pa <- comb_res$preds[order(comb_res$row_id)]

    # predicted presence probability
    pres_prob <- 1 - exp(-exp(pred.pa))
    # NEED TO TRUNCATE PRESENCE PROBABILITY BY 'trunc.pa.prob'
    pres_prob[pres_prob <= trunc.pa.prob] <- trunc.pa.prob
    pres_prob[1 - pres_prob <= trunc.pa.prob] <- 1 - trunc.pa.prob
    # predicted conditional log-likelihood on the PA data
    pll_pa <- sum((pa.resp * log(pres_prob)) + ((1 - pa.resp) * log(1 - pres_prob)))
    # predicted conditional log-likelihood on the PO data
    pll_po <- NA
    # area under ROC curve
    roccurve <- pROC::roc(pa.resp, as.vector(pres_prob), quiet = T)
    auc_val <- as.numeric(pROC::auc(roccurve))

  } else {
    stop(paste0("model type not reognised: ", object$model.type))
  }

  # return the CV predictions and model
  return(list(predicted.cll.pa = pll_pa, predicted.cll.po = pll_po, auc = auc_val, model = object))
}
