#' Deprecated Functions
#'
#' Functions that have been deprecated and will be removed in the future.
#'
#' @name deprecated
#' @rdname deprecated
#'
#' @param ... arguments passed to non-deprecated equivalent.
#'
#' @noRd
NULL


dep_methodarg <- function(method) {
  if (!missing(method)) {
    warning("Argument 'method' is deprecated; use 'cov' instead.",
            call. = FALSE)
  }
}