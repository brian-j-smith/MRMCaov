#' Deprecated Functions
#'
#' Functions that have been deprecated and will be removed in the future.
#'
#' @name deprecated
#' @rdname deprecated
#'
#' @param ... arguments passed to non-deprecated equivalent.
#'
NULL


#' @rdname deprecated
#'
#' @details
#' Use \code{\link[=binormalLR_auc]{binormalLR_auc()}} instead of
#' \code{proproc_auc()}.
#'
proproc_auc <- function(...) {
  warning("proproc_auc() is deprecated; use binormalLR_auc() instead")
  binormalLR_auc(...)
}


#' @rdname deprecated
#'
#' @details
#' Use \code{\link[=binormalLR_sens]{binormalLR_sens()}} instead of
#' \code{proproc_sens()}.
#'
proproc_sens <- function(...) {
  warning("proproc_sens() is deprecated; use binormalLR_sens() instead")
  binormalLR_sens(...)
}


#' @rdname deprecated
#'
#' @details
#' Use \code{\link[=binormalLR_spec]{binormalLR_spec()}} instead of
#' \code{proproc_spec()}.
#'
proproc_spec <- function(...) {
  warning("proproc_spec() is deprecated; use binormalLR_spec() instead")
  binormalLR_spec(...)
}
