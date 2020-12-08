#' Performance Plots
#'
#' @name plot
#' @rdname plot-methods
#'
#' @param x object to plot.
#' @param n number of equally spaced false-positive rate points at which to
#'   calculate true-positive rates and interpolate through for display of the
#'   curve.
#' @param emp_points logical indicating whether to overlay empirical ROC points
#'   on parametric curves.
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
plot.roc_curve <- function(x, n = 100, emp_points = FALSE, ...) {
  emp_points <- if (emp_points) {
    points(roc_curves(x, method = "empirical"), which = "observed")
  }
  x <- points(x, values = seq(0, 1, length = n))
  plot_roc_points(x, emp_points = emp_points, ...)
}


#' @rdname plot-methods
#'
plot.roc_curves <- function(x, n = 100, emp_points = FALSE, ...) {
  emp_points <- if (emp_points) {
    points(roc_curves(x, method = "empirical"), which = "observed")
  }
  x <- points(x, values = seq(0, 1, length = n))
  plot_roc_points(x, emp_points = emp_points, ...)
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
  plot_roc_points(x, coord_fixed = coord_fixed, ...)
}


plot_roc_points <- function(x, coord_fixed = TRUE, emp_points = NULL, ...) {
  x_df <- make_plot_df(x)
  emp_df <- make_plot_df(emp_points)

  aes_args <- list(~ FPR, ~ TPR)
  if (!is.null(x_df[["Group"]])) aes_args$color <- ~ Group

  p <- ggplot(x_df, do.call(aes_, aes_args)) +
    geom_path() +
    labs(x = "False Positive Rate", y = "True Positive Rate",
         color = tail(names(x[["Group"]]), 1))

  if (coord_fixed) p <- p + coord_fixed()
  if (!is.null(emp_points)) {
    p <- p + geom_point(data = emp_df, inherit.aes = TRUE)
  }
  if (!is.null(x_df$facet)) p <- p + facet_wrap(~ facet)

  p
}


make_plot_df <- function(x) {
  df <- data.frame(
    FPR = x$FPR,
    TPR = x$TPR
  )

  group_vars <- x[["Group"]]
  n_group_vars <- length(group_vars)
  if (n_group_vars) {
    df$Group <- group_vars[[n_group_vars]]
    group_vars[[n_group_vars]] <- NULL
    if (length(group_vars)) {
      prefix <- paste(names(group_vars), collapse = ".")
      df$facet <- paste0(prefix, ": ", interaction(group_vars))
    }
  }

  df
}


#' @rdname plot-methods
#'
plot.mrmc <- function(x, n = 100, ...) {
  plot(roc_curves(x), n = n)
}


#' @rdname plot-methods
#'
plot.stmc <- function(x, n = 100, ...) {
  plot(roc_curves(x), n = n)
}
