#' Performance Plots
#'
#' @name plot
#' @rdname plot-methods
#'
#' @param x object to plot.
#' @param n number of equally spaced false-positive rate points at which to
#'   calculate true-positive rates and interpolate through for display of the
#'   curve.
#' @param coord_fixed logical indicating whether to fix the scales of x and y
#'   axes.
#' @param ... arguments passed to other methods.
#'
#' @seealso \code{\link{roc_curves}}
#'
#' @examples
#' curves <- with(VanDyke,
#'   roc_curves(truth, rating, groups = list(Test = treatment, Reader = reader))
#' )
#' plot(curves)
#'
NULL


#' @rdname plot-methods
#'
plot.roc_curve <- function(x, n = 100, ...) {
  plot(points(x, values = seq(0, 1, length = n)))
}


#' @rdname plot-methods
#'
plot.roc_curves <- function(x, n = 100, ...) {
  plot(points(x, values = seq(0, 1, length = n)))
}


#' @rdname plot-methods
#'
plot.empirical_curve <- function(x, ...) {
  plot(points(x))
}


#' @rdname plot-methods
#'
plot.empirical_curves <- function(x, ...) {
  plot(points(x))
}


#' @rdname plot-methods
#'
plot.roc_points <- function(x, coord_fixed = TRUE, ...) {
  df <- data.frame(
    FPR = x$FPR,
    TPR = x$TPR
  )
  aes_args <- list(~ FPR, ~ TPR)

  group_vars <- x[["Group"]]
  n_group_vars <- length(group_vars)
  if (n_group_vars) {
    df$Group <- group_vars[[n_group_vars]]
    aes_args$color <- ~ Group
    group_vars[[n_group_vars]] <- NULL
    if (length(group_vars)) {
      prefix <- paste(names(group_vars), collapse = ".")
      df$facet <- paste0(prefix, ": ", interaction(group_vars))
    }
  }

  p <- ggplot(df, do.call(aes_, aes_args)) +
    geom_path() +
    labs(x = "False Positive Rate", y = "True Positive Rate",
         color = names(x[["Group"]])[n_group_vars])

  if (coord_fixed) p <- p + coord_fixed()

  if (!is.null(df$facet)) p <- p + facet_wrap(~ facet)

  p
}


#' @rdname plot-methods
#'
plot.mrmc <- function(x, n = 100, ...) {
  plot(roc_curves(x), n = n)
}
