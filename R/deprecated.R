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


dep_methodarg <- function(missing) {
  if (!missing) {
    warning("Argument 'method' is deprecated; use 'cov' instead.",
            call. = FALSE)
  }
}