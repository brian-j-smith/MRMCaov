#' Performance Metrics
#'
#' Estimated performance metrics from ROC curves.
#'
#' @name metrics
#' @rdname metrics
#'
#' @param truth vector of true binary statuses.
#' @param rating vector of numeric ratings.
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
binormal_auc <- function(truth, rating, partial = FALSE, min = 0, max = 1) {
  fit <- if (isFALSE(partial)) {
    truth <- as.factor(truth)
    is_pos <- truth == levels(truth)[2]
    pred_pos <- as.double(rating[is_pos])
    pred_neg <- as.double(rating[!is_pos])
    .Fortran("cvbmroc",
             length(pred_neg), length(pred_pos), pred_neg, pred_pos,
             double(1), double(1), est = double(1), var = double(1))
  } else {
    partial <- match.arg(partial, c("sensitivity", "specificity"))
    min <- as.double(min)
    max <- as.double(max)
    if (partial == "specificity") {
      flag <- 1L
      min <- 1 - max
      max <- 1 - min
    } else {
      flag <- 2L
    }
    params <- binormal_params(truth, rating)
    .Fortran("cvbmrocpartial",
             params$a, params$b, min, max, flag, est = double(1),
             err = integer(1))
  }

  fit$est
}


#' @rdname metrics
#'
binormal_sens <- function(truth, rating, spec) {
  x <- proproc_params(truth, rating)
  sensitivity(x, specificity = spec)
}


#' @rdname metrics
#'
binormal_spec <- function(truth, rating, sens) {
  x <- proproc_params(truth, rating)
  specificity(x, sensitivity = sens)
}


#' @rdname metrics
#'
empirical_auc <- function(truth, rating, partial = FALSE, min = 0, max = 1) {
  args <- list(pROC::roc(truth, rating, quiet = TRUE))
  if (!isFALSE(partial)) {
    partial <- match.arg(partial, c("sensitivity", "specificity"))
    args$partial.auc <- c(min, max)
    args$partial.auc.focus <- partial
  }
  as.numeric(do.call(pROC::auc, args))
}


#' @rdname metrics
#'
empirical_sens <- function(truth, rating, spec) {
  roc <- pROC::roc(truth, rating, quiet = TRUE)
  pROC::coords(roc, spec, input = "specificity", ret = "sensitivity",
               transpose = TRUE, drop = TRUE)
}


#' @rdname metrics
#'
empirical_spec <- function(truth, rating, sens) {
  roc <- pROC::roc(truth, rating, quiet = TRUE)
  pROC::coords(roc, sens, input = "sensitivity", ret = "specificity",
               transpose = TRUE, drop = TRUE)
}


#' @rdname metrics
#'
proproc_auc <- function(truth, rating, partial = FALSE, min = 0, max = 1) {
  fit <- if (isFALSE(partial)) {
    truth <- as.factor(truth)
    is_pos <- truth == levels(truth)[2]
    pred_pos <- as.double(rating[is_pos])
    pred_neg <- as.double(rating[!is_pos])
    .Fortran("pbmroc",
             length(pred_neg), length(pred_pos), pred_neg, pred_pos,
             double(1), double(1), est = double(1), var = double(1))
  } else {
    partial <- match.arg(partial, c("sensitivity", "specificity"))
    min <- as.double(min)
    max <- as.double(max)
    if (partial == "specificity") {
      flag <- 1L
      min <- 1 - max
      max <- 1 - min
    } else {
      flag <- 2L
    }
    params <- proproc_params(truth, rating)
    .Fortran("pbmrocpartial",
             params$d_a, params$c, min, max, flag, est = double(1),
             err = integer(1))
  }

  fit$est
}


#' @rdname metrics
#'
proproc_sens <- function(truth, rating, spec) {
  x <- proproc_params(truth, rating)
  sensitivity(x, specificity = spec)
}


#' @rdname metrics
#'
proproc_spec <- function(truth, rating, sens) {
  x <- proproc_params(truth, rating)
  specificity(x, sensitivity = sens)
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
