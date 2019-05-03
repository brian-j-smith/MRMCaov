psi <- function(x_pos, x_neg) {
  (x_pos > x_neg) + 0.5 * (x_pos == x_neg)
}


#' Performance Metrics
#' 
#' @name metrics
#' @rdname metrics
#' 
#' @param observed vector of observed binary responses.
#' @param predicted vector of predcited responses.
#' 
proproc_auc <- function(observed, predicted) {
  observed <- as.factor(observed)
  is_pos <- observed == levels(observed)[2]

  pred_pos <- as.double(predicted[is_pos])
  pred_neg <- as.double(predicted[!is_pos])
  
  fit <- .Fortran("pbmroc",
                  length(pred_neg), length(pred_pos), pred_neg, pred_pos,
                  double(1), double(1), est = double(1), var = double(1))
  fit$est
}


#' @rdname metrics
#' 
roc_auc <- function(observed, predicted) {
  observed <- as.factor(observed)
  
  R <- rank(predicted)
  is_pos <- observed == levels(observed)[2]
  num_pos <- sum(is_pos)
  num_neg <- length(observed) - num_pos
  (sum(R[is_pos]) - num_pos * (num_pos + 1) / 2) / (num_pos * num_neg)
}


#' ROC Performance Metrics
#' 
#' Calculation of TPR and FPR pairs for all values of a numeric predictor
#' of an observed binary response.
#' 
#' @param observed vector of observed binary responses.
#' @param predicted vector of predcited responses.
#' @param ... interaction factors within which to calculate metrics.
#' 
#' @seealso \code{\link{plot}}
#' 
#' @examples
#' with(VanDyke, roc(truth, rating, treatment, reader))
#' 
roc <- function(observed, predicted, ...) {
  groups <- list(...)
  if (length(groups) == 0) groups <- list(rep(1, length(observed)))
  if (is.null(names(groups))) {
    names(groups) <- make.unique(rep("Group", length(groups)))
  }
  groups <- lapply(groups, as.character)
  
  df <- data.frame(observed, predicted, groups)
  perf_list <- by(df, groups, function(data) {
    metrics <- .roc(data$observed, data$predicted)
    data.frame(data[rep(1, nrow(metrics)), -(1:2), drop = FALSE],
               metrics, row.names = NULL)
  }, simplify = FALSE)
  names(perf_list) <- NULL

  structure(
    do.call(rbind, perf_list),
    class = c("roc", "data.frame")
  )
}


.roc <- function(observed, predicted) {
  perf <- ROCR::prediction(predicted, observed) %>%
    ROCR::performance("tpr", "fpr")
  cbind(FPR = perf@x.values[[1]], TPR = perf@y.values[[1]])
}
