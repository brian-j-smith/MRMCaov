#' Deprecated Functions
#'
#' Functions that have been deprecated and will be removed in the future.
#'
#' @name deprecated
#' @rdname deprecated
#'
#' @param truth,rating,... arguments passed to non-deprecated equivalent.
#'
NULL


#' @rdname deprecated
#'
#' @details
#' Use \code{\link[=roc_curves]{roc_curves()}} instead of \code{roc()}.
#'
roc <- function(truth, rating, ...) {
  warning("roc() is deprecated; use roc_curves() instead")
  roc_curves(truth, rating, list(...), method = "empirical")
}


#' @rdname deprecated
#'
#' @details
#' Use \code{\link[=roc_curves]{roc_curves()}} instead of \code{proproc()}.
#'
proproc <- function(truth, rating, ...) {
  warning("proproc() is deprecated; use roc_curves() instead")
  roc_curves(truth, rating, list(...), method = "proproc")
}
