#' Performance Metrics
#'
#' Estimated performance metrics from ROC curves.
#'
#' @name metrics
#' @rdname metrics
#'
#' @param truth vector of true binary statuses.
#' @param rating vector of 0-1 binary ratings for the binary metrics and ranges
#'   of numeric ratings for the others.
#' @param partial character string \code{"sensitivity"} or \code{"specificity"}
#'   for calculation of partial AUC, or \code{FALSE} for full AUC.  Partial
#'   matching of the character strings is allowed.
#' @param min,max minimum and maximum sensitivity or specificity values over
#'   which to calculate partial AUC.
#' @param sens,spec numeric sensitivity/specificity at which to calculate
#'   specificity/sensitivity.
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


#' @rdname metrics
#'
binary_sens <- function(truth, rating) {
  binary_metric(truth, rating, function(truth, rating) {
    events <- truth == levels(truth)[2]
    pos <- rating == levels(truth)[2]
    mean(pos[events])
  })
}


#' @rdname metrics
#'
binary_spec <- function(truth, rating) {
  binary_metric(truth, rating, function(truth, rating) {
    nonevents <- truth != levels(truth)[2]
    neg <- rating != levels(truth)[2]
    mean(neg[nonevents])
  })
}


binary_metric <- function(truth, rating, f) {
  truth <- as.factor(truth)
  stopifnot(nlevels(truth) == 2)
  if (!is.factor(rating)) {
    rating <- factor(as.numeric(rating), levels = 0:1, labels = levels(truth))
  }
  stopifnot(nlevels(rating) == 2)
  f(truth, rating)
}


#' @rdname metrics
#'
binormal_auc <- function(truth, rating, partial = FALSE, min = 0, max = 1) {
  curve <- roc_curves(truth, rating, method = "binormal")
  auc(curve, partial = partial, min = min, max = max)
}


#' @rdname metrics
#'
binormal_sens <- function(truth, rating, spec) {
  curve <- roc_curves(truth, rating, method = "binormal")
  sensitivity(curve, specificity = spec)
}


#' @rdname metrics
#'
binormal_spec <- function(truth, rating, sens) {
  curve <- roc_curves(truth, rating, method = "binormal")
  specificity(curve, sensitivity = sens)
}


#' @rdname metrics
#'
empirical_auc <- function(truth, rating, partial = FALSE, min = 0, max = 1) {
  curve <- roc_curves(truth, rating, method = "empirical")
  auc(curve, partial = partial, min = min, max = max)
}


#' @rdname metrics
#'
empirical_sens <- function(truth, rating, spec) {
  curve <- roc_curves(truth, rating, method = "empirical")
  sensitivity(curve, specificity = spec)
}


#' @rdname metrics
#'
empirical_spec <- function(truth, rating, sens) {
  curve <- roc_curves(truth, rating, method = "empirical")
  specificity(curve, sensitivity = sens)
}


#' @rdname metrics
#'
proproc_auc <- function(truth, rating, partial = FALSE, min = 0, max = 1) {
  curve <- roc_curves(truth, rating, method = "proproc")
  auc(curve, partial = partial, min = min, max = max)
}


#' @rdname metrics
#'
proproc_sens <- function(truth, rating, spec) {
  curve <- roc_curves(truth, rating, method = "empirical")
  sensitivity(curve, specificity = spec)
}


#' @rdname metrics
#'
proproc_spec <- function(truth, rating, sens) {
  curve <- roc_curves(truth, rating, method = "empirical")
  specificity(curve, sensitivity = sens)
}


#' @rdname metrics
#'
trapezoidal_auc <- function(truth, rating, partial = FALSE, min = 0, max = 1) {
  empirical_auc(truth, rating, partial = partial, min = min, max = max)
}


#' @rdname metrics
#'
trapezoidal_sens <- function(truth, rating, spec) {
  empirical_sens(truth, rating, spec)
}


#' @rdname metrics
#'
trapezoidal_spec <- function(truth, rating, sens) {
  empirical_spec(truth, rating, sens)
}


psi <- function(x_pos, x_neg) {
  (x_pos > x_neg) + 0.5 * (x_pos == x_neg)
}
