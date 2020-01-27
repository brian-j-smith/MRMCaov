#' ROC Performance Metrics
#'
#' Calculation of TPR and FPR pairs for all values of a numeric rating of a
#' true binary response.
#'
#' @param truth vector of true binary statuses.
#' @param rating vector of numeric ratings.
#' @param ... interaction factors within which to calculate metrics.
#'
#' @rdname roc
#'
#' @seealso \code{\link{plot}}
#'
#' @examples
#' with(VanDyke, roc(truth, rating, treatment, reader))
#'
roc <- function(truth, rating, ...) {
  groups <- list(...)
  if (length(groups) == 0) groups <- list(rep(1, length(truth)))
  if (is.null(names(groups))) {
    names(groups) <- make.unique(rep("Group", length(groups)))
  }
  groups <- lapply(groups, as.character)

  df <- data.frame(truth, rating, groups)
  perf_list <- by(df, groups, function(data) {
    metrics <- .roc(data$truth, data$rating)
    data.frame(data[rep(1, nrow(metrics)), -(1:2), drop = FALSE],
               metrics, row.names = NULL)
  }, simplify = FALSE)
  names(perf_list) <- NULL

  structure(
    do.call(rbind, perf_list),
    class = c("roc_frame", "data.frame")
  )
}


.roc <- function(truth, rating) {
  roc <- pROC::roc(truth, rating, auc = FALSE, quiet = TRUE)
  cbind(FPR = 1 - roc$specificities, TPR = roc$sensitivities)
}


proproc_params <- function(truth, rating) {
  truth <- as.factor(truth)
  is_pos <- truth == levels(truth)[2]

  pred_pos <- as.double(rating[is_pos])
  pred_neg <- as.double(rating[!is_pos])

  auc <- .Fortran("pbmroc",
                  length(pred_neg), length(pred_pos), pred_neg, pred_pos,
                  d_a = double(1), c = double(1),
                  est = double(1), var = double(1))

  structure(list(d_a = auc$d_a, c = auc$c), class = "proproc_params")
}


#' @rdname roc
#'
proproc <- function(truth, rating, ...) {
  groups <- list(...)
  if (length(groups) == 0) groups <- list(rep(1, length(truth)))
  if (is.null(names(groups))) {
    names(groups) <- make.unique(rep("Group", length(groups)))
  }
  groups <- lapply(groups, as.character)

  all_data <- data.frame(truth, rating)
  proproc_by <- by(all_data, groups, function(data) {
    proproc_params(data$truth, data$rating)
  }, simplify = FALSE)
  params <- data.frame(row.names = seq(proproc_by))
  params$coef <- as.data.frame(t(sapply(proproc_by, unlist)))
  params$group <- expand.grid(dimnames(proproc_by))

  roc <- roc(truth, rating, ...)
  groups <- roc[!(names(roc) %in% c("FPR", "TPR"))]
  roc_split <- split(roc, groups)
  for (i in seq(roc_split)) {
    proproc <- structure(as.list(params$coef[i, ]), class = "proproc_params")
    spec <- 1 - roc_split[[i]]$FPR
    roc_split[[i]]$TPR <- sensitivity(proproc, spec)
  }
  roc <- unsplit(roc_split, groups)

  structure(roc, params = params, class = c("proproc_frame", "data.frame"))
}


sensitivity <- function(x, ...) {
  UseMethod("sensitivity")
}


sensitivity.proproc_params <- function(x, specificity, ...) {
  sapply(specificity, function(spec) {
    fpf <- 1.0 - spec
    .Fortran("pbmrocfpf2tpf",
             x$d_a, x$c,
             fpf = fpf, tpf = double(1), double(1))$tpf
  })
}


specificity <- function(x, ...) {
  UseMethod("specificity")
}


specificity.proproc_params <- function(x, sensitivity, ...) {
  sapply(sensitivity, function(sens) {
    1 - .Fortran("pbmroctpf2fpf",
                 x$d_a, x$c,
                 tpf = as.double(sens), fpf = as.double(1), double(1))$fpf
  })
}
