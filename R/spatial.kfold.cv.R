#' Performs a spatial hold-one-out K fold cross validation on a scampr model
#'
#' @param object a scampr model object
#' @param po.fold.id An integer or factor vector the same length as the presence-only data used in the model \code{object} - describes the CV fold that each location falls into.
#' @param pa.fold.id An integer or factor vector the same length as the presence/absence data used in the model \code{object} - describes the CV fold that each location falls into.
#' @param trunc.pa.prob A small positive number by which the predicted probability of presence is truncated. This can be used to ensure infinite values are avoided within the cross-validation.
#'
#' @return a list including- 'predicted.cll.pa': the predicted logLikelihood, conditional on the model predictors and latent field, on the held out presence/absence data. 'predicted.cll.po': the predicted logLikelihood, conditional on the model predictors and latent field, on the held out presence-only data. 'auc': the area under the ROC curve on the held out presence/absence sites. 'model': the model upon which the cross validation was performed.
#' @export
#'
#' @importFrom pROC roc auc
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
spatial.kfold.cv <- function(object, po.fold.id, pa.fold.id, trunc.pa.prob = 1e-7) {

  if (object$model.type == "IDM") {

    # check the appropriate spatial folds are present
    if (missing(po.fold.id) & missing(pa.fold.id)) {
      stop(paste0("'po.fold.id' and 'pa.fold.id' must be supplied for a model of type = ", object$model.type))
    }

    # Convert CV folds to numerics in case factor is provided
    po.fold.id <- as.numeric(po.fold.id)
    pa.fold.id <- as.numeric(pa.fold.id)

    # Extract the data
    po.data <- object$data
    pa.data <- attr(object$data, "PA")
    po.resp <- po.data[, all.vars(object$formula[[2]])]
    pa.resp <- pa.data[, all.vars(object$formula[[2]])]
    quad.sizes <- po.data[, object$quad.weights.name]

    # Initialise Objects #
    train_mods <- list() # training models
    test.pa <- NULL # predicted abundance rate on PA data
    test.po <- NULL # predicted intensity rate on PO data
    po.test.rows <- NULL # storage for row indices
    pa.test.rows <- NULL # storage for row indices
    po.row.id <- 1:nrow(po.data) # index reference for returning original PO dataset order
    pa.row.id <- 1:nrow(pa.data) # index reference for returning original PA dataset order

    # Loop through CV folds
    for (i in sort(unique(po.fold.id))) {
      # get the model's call
      call.list <- as.list(object$call)
      # adjust the training data
      call.list$data <- po.data[po.fold.id != i, ]
      call.list$IDM.presence.absence.df <- pa.data[pa.fold.id != i, ]
      # train the model
      call.list[[1]] <- NULL
      train_mods[[i]] <- do.call("scampr", call.list)
      # test on held out fold
      test.pa <- c(test.pa, predict.scampr(train_mods[[i]], newdata = pa.data[pa.fold.id == i, ]))
      test.po <- c(test.po, predict.scampr(train_mods[[i]], newdata = po.data[po.fold.id == i, ], include.bias.accounting = TRUE))
      # calculate the row indices of the new test datasets
      po.test.rows <- c(po.test.rows, po.row.id[po.fold.id == i])
      pa.test.rows <- c(pa.test.rows, pa.row.id[pa.fold.id == i])
    }
    # re-order the test datasets to match the raw datasets
    pred.po <- test.po[order(po.test.rows)]
    pred.pa <- test.pa[order(pa.test.rows)]

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

    # Extract the data
    po.data <- object$data
    po.resp <- po.data[, all.vars(object$formula[[2]])]
    quad.sizes <- po.data[, object$quad.weights.name]

    # Initialise Objects #
    train_mods <- list() # training models
    test.po <- NULL # predicted intensity rate on PO data
    po.test.rows <- NULL # storage for row indices
    po.row.id <- 1:nrow(po.data) # index reference for returning original PO dataset order

    # Loop through CV folds
    for (i in sort(unique(po.fold.id))) {
      # get the model's call
      call.list <- as.list(object$call)
      # adjust the training data
      call.list$data <- po.data[po.fold.id != i, ]
      # train the model
      call.list[[1]] <- NULL
      train_mods[[i]] <- do.call("scampr", call.list)
      # test on held out fold
      test.po <- c(test.po, predict.scampr(train_mods[[i]], newdata = po.data[po.fold.id == i, ], include.bias.accounting = TRUE))
      # calculate the row indices of the new test dataset
      po.test.rows <- c(po.test.rows, po.row.id[po.fold.id == i])
    }
    # re-order the test datasets to match the raw dataset
    pred.po <- test.po[order(po.test.rows)]

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

    # Extract the data
    pa.data <- object$data
    pa.resp <- pa.data[, all.vars(object$formula[[2]])]

    # Initialise Objects #
    train_mods <- list() # training models
    test.pa <- NULL # predicted abundance rate on PA data
    pa.test.rows <- NULL # storage for row indices
    pa.row.id <- 1:nrow(pa.data) # index reference for returning original PA dataset order

    # Loop through CV folds
    for (i in sort(unique(pa.fold.id))) {
      # get the model's call
      call.list <- as.list(object$call)
      # adjust the training data
      call.list$data <- pa.data[pa.fold.id != i, ]
      # train the model
      call.list[[1]] <- NULL
      train_mods[[i]] <- do.call("scampr", call.list)
      # test on held out fold
      test.pa <- c(test.pa, predict.scampr(train_mods[[i]], newdata = pa.data[pa.fold.id == i, ]))
      # calculate the row indices of the new test dataset
      pa.test.rows <- c(pa.test.rows, pa.row.id[pa.fold.id == i])
    }
    # re-order the test datasets to match the raw dataset
    pred.pa <- test.pa[order(pa.test.rows)]

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
