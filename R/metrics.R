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
#' \item{empirical_auc, trapezoidal_auc:}{Area under the empirical, or
#' trapezoidal, ROC curve.}
#' \item{proproc_auc:}{Area under a proper ROC curve.}
#' }
#' 
NULL


binormal_auc <- function(truth, rating) {
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
empirical_auc <- function(truth, rating) {
  truth <- as.factor(truth)
  
  R <- rank(rating)
  is_pos <- truth == levels(truth)[2]
  num_pos <- sum(is_pos)
  num_neg <- length(truth) - num_pos
  (sum(R[is_pos]) - num_pos * (num_pos + 1) / 2) / (num_pos * num_neg)
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
trapezoidal_auc <- function(truth, rating) {
  empirical_auc(truth, rating)
}


psi <- function(x_pos, x_neg) {
  (x_pos > x_neg) + 0.5 * (x_pos == x_neg)
}
