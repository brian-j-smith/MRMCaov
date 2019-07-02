psi <- function(x_pos, x_neg) {
  (x_pos > x_neg) + 0.5 * (x_pos == x_neg)
}


#' Performance Metrics
#' 
#' @name metrics
#' @rdname metrics
#' 
#' @param truth vector of true binary statuses.
#' @param rating vector of numeric ratings.
#' 
#' @details
#' Performance metrics measure the degree to which higher case ratings are
#' associated with positive case statuses, where positive status is taken to be
#' the highest level of \code{truth}.  Available metrics include the following.
#' \describe{
#' \item{binormroc_auc:}{Area under the binormal ROC curve.}
#' \item{proproc_auc:}{Area under a proper ROC curve.}
#' \item{roc_auc:}{Area under the empirical, or trapezoidal, ROC curve.}
#' }
#' 
binormroc_auc <- function(truth, rating) {
  truth <- as.factor(truth)
  is_pos <- truth == levels(truth)[2]

  pred_pos <- as.double(rating[is_pos])
  pred_neg <- as.double(rating[!is_pos])
  
  fit <- .Fortran("cvbmroc",
                  length(pred_neg), length(pred_pos), pred_neg, pred_pos,
                  double(1), double(1), est = double(1), var = double(1))
  fit$est
}


#' @rdname metrics
#' 
proproc_auc <- function(truth, rating) {
  truth <- as.factor(truth)
  is_pos <- truth == levels(truth)[2]

  pred_pos <- as.double(rating[is_pos])
  pred_neg <- as.double(rating[!is_pos])
  
  fit <- .Fortran("pbmroc",
                  length(pred_neg), length(pred_pos), pred_neg, pred_pos,
                  double(1), double(1), est = double(1), var = double(1))
  fit$est
}


#' @rdname metrics
#' 
roc_auc <- function(truth, rating) {
  truth <- as.factor(truth)
  
  R <- rank(rating)
  is_pos <- truth == levels(truth)[2]
  num_pos <- sum(is_pos)
  num_neg <- length(truth) - num_pos
  (sum(R[is_pos]) - num_pos * (num_pos + 1) / 2) / (num_pos * num_neg)
}


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
