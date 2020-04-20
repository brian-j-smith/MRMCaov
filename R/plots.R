#' Performance Plots
#'
#' @name plot
#' @rdname plot-methods
#'
#' @param x object to plot.
#' @param n number of equally spaced false-positive rate points at which to
#'   calculate true-positive rates and interpolate through for display of the
#'   curve.
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
plot.empirical_curves <- function(x, ...) {
  plot(points(x))
}


#' @rdname plot-methods
#'
plot.param_curves <- function(x, n = 100, ...) {
  plot(points(x, values = seq(0, 1, length = n)))
}


#' @rdname plot-methods
#'
plot.roc_points <- function(x, ...) {
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
    coord_fixed() +
    labs(x = "False Positive Rate", y = "True Positive Rate",
         color = names(x[["Group"]])[n_group_vars])

  if (!is.null(df$facet)) p <- p + facet_wrap(~ facet)

  p
}


#' @rdname plot-methods
#'
plot.mrmc <- function(x, n = 100, ...) {
  plot(x$roc, n = n)
}
