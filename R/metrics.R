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
#'   matching of the character strings is allowed.  \code{"specificity"}
#'   results in area under the ROC curve between the given \code{min} and
#'   \code{max} specificity values, whereas \code{"sensitivity"} results in area to
#'   the right of the curve between the given sensitivity values.
#' @param min,max minimum and maximum sensitivity or specificity values over
#'   which to calculate partial AUC.
#' @param normalize logical indicating whether partial AUC is divided by the
#'   interval width (\code{max - min}) over which it is calculated.
#' @param sens,spec numeric sensitivity/specificity at which to calculate
#'   specificity/sensitivity.
#' @param slope slope of the iso-utility line at which to compute expected
#'   utility of the ROC curve.
#'
#' @details
#' Performance metrics measure the degree to which higher case ratings are
#' associated with positive case statuses, where positive status is taken to be
#' the highest level of \code{truth}.  Available metrics include area under the
#' ROC curve (auc), expected utility of the ROC curve (eu) at a given
#' iso-utility line (Abbey, 2013), sensitivity (sens) at a given specificity,
#' and specificity (spec) at a given sensitivity.
#'
#' @return
#' Returns a numeric value.
#'
#' @seealso \code{\link{mrmc}}, \code{\link{srmc}}, \code{\link{stmc}}
#'
#' @references
#' Abbey CK, Samuelson FW and Gallas BD (2013). Statistical power considerations
#' for a utility endpoint in observer performance studies. Academic Radiology,
#' 20: 798-806.
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
binormal_auc <- function(truth, rating, partial = FALSE, min = 0, max = 1,
                         normalize = FALSE) {
  curve <- roc_curves(truth, rating, method = "binormal")
  auc(curve, partial = partial, min = min, max = max, normalize = normalize)
}


#' @rdname metrics
#'
binormal_eu <- function(truth, rating, slope = 1) {
  curve <- roc_curves(truth, rating, method = "binormal")
  roc_eu(curve, slope = slope)
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
binormalLR_auc <- function(truth, rating, partial = FALSE, min = 0, max = 1,
                           normalize = FALSE) {
  curve <- roc_curves(truth, rating, method = "binormalLR")
  auc(curve, partial = partial, min = min, max = max, normalize = normalize)
}


#' @rdname metrics
#'
binormalLR_eu <- function(truth, rating, slope = 1) {
  curve <- roc_curves(truth, rating, method = "binormalLR")
  roc_eu(curve, slope = slope)
}


#' @rdname metrics
#'
binormalLR_sens <- function(truth, rating, spec) {
  curve <- roc_curves(truth, rating, method = "binormalLR")
  sensitivity(curve, specificity = spec)
}


#' @rdname metrics
#'
binormalLR_spec <- function(truth, rating, sens) {
  curve <- roc_curves(truth, rating, method = "binormalLR")
  specificity(curve, sensitivity = sens)
}


#' @rdname metrics
#'
empirical_auc <- function(truth, rating, partial = FALSE, min = 0, max = 1,
                          normalize = FALSE) {
  curve <- roc_curves(truth, rating, method = "empirical")
  auc(curve, partial = partial, min = min, max = max, normalize = normalize)
}


#' @rdname metrics
#'
empirical_eu <- function(truth, rating, slope = 1) {
  curve <- roc_curves(truth, rating, method = "empirical")
  roc_eu(curve, slope = slope)
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
trapezoidal_auc <- function(truth, rating, partial = FALSE, min = 0, max = 1,
                            normalize = FALSE) {
  do.call(empirical_auc, as.list(environment()))
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
  (sign(x_pos - x_neg) + 1) / 2
}
