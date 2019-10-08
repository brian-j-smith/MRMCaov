#' ROC Performance Metrics
#' 
#' Calculation of TPR and FPR pairs for all values of a numeric rating of a
#' true binary response.
#' 
#' @param truth vector of true binary statuses.
#' @param rating vector of numeric ratings.
#' @param ... interaction factors within which to calculate metrics.
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
    class = c("roc", "data.frame")
  )
}


.roc <- function(truth, rating) {
  perf <- ROCR::prediction(rating, truth) %>%
    ROCR::performance("tpr", "fpr")
  cbind(FPR = perf@x.values[[1]], TPR = perf@y.values[[1]])
}


proproc <- function(truth, rating) {
  truth <- as.factor(truth)
  is_pos <- truth == levels(truth)[2]
  
  pred_pos <- as.double(rating[is_pos])
  pred_neg <- as.double(rating[!is_pos])
  
  auc <- .Fortran("pbmroc",
                  length(pred_neg), length(pred_pos), pred_neg, pred_pos,
                  d_a = double(1), c = double(1),
                  est = double(1), var = double(1))
  
  structure(list(d_a = auc$d_a, c = auc$c), class = "proproc")
}


sensitivity <- function(x, ...) {
  UseMethod("sensitivity")
}


sensitivity.proproc <- function(x, specificity, ...) {
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


specificity.proproc <- function(x, sensitivity, ...) {
  sapply(sensitivity, function(sens) {
    1 - .Fortran("pbmroctpf2fpf",
                 x$d_a, x$c,
                 tpf = as.double(sens), fpf = as.double(1), double(1))$fpf
  })
}
