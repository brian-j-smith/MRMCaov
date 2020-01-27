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
#' @seealso \code{\link{roc}}
#'
#' @examples
#' perf <- with(VanDyke, roc(truth, rating, treatment, reader))
#' plot(perf)
#'
plot.roc_frame <- function(x, ...) {
  df <- data.frame(
    x = x$FPR,
    y = x$TPR
  )
  aes_args <- list(~ x, ~ y)

  x[c("FPR", "TPR")] <- NULL
  n <- ncol(x)

  if (n > 0) {
    df$group <- x[[n]]
    aes_args$color <- ~ group
  }
  if (n > 1) {
    prefix <- paste(names(x)[-n], collapse = ".")
    df$facet <- paste0(prefix, ": ", interaction(x[-n]))
  }

  p <- ggplot(df, do.call(aes_, aes_args)) +
    geom_path() +
    coord_fixed() +
    labs(x = "False Positive Rate", y = "True Positive Rate",
         color = names(x)[n])

  if (!is.null(df$facet)) p <- p + facet_wrap(~ facet)

  p
}


#' @rdname plot-methods
#'
plot.proproc_frame <- function(x, n = 101, ...) {
  x <- attr(x, "params")
  params <- split(x$coef, seq_len(nrow(x)))
  roc_list <- lapply(params, function(param) {
    proproc <- structure(as.list(param), class = "proproc_params")
    fpr <- seq(0, 1, length = n)
    tpr <- sensitivity(proproc, specificity = 1 - fpr)
    cbind(FPR = fpr, TPR = tpr)
  })

  data <- cbind(
    do.call(rbind, roc_list),
    x$group[rep(seq_len(nrow(x)), each = n), , drop = FALSE]
  )

  df <- data.frame(
    x = data$FPR,
    y = data$TPR
  )
  aes_args <- list(~ x, ~ y)

  data[c("FPR", "TPR")] <- NULL
  n <- ncol(data)

  if (n > 0) {
    df$group <- data[[n]]
    aes_args$color <- ~ group
  }
  if (n > 1) {
    prefix <- paste(names(data)[-n], collapse = ".")
    df$facet <- paste0(prefix, ": ", interaction(data[-n]))
  }

  p <- ggplot(df, do.call(aes_, aes_args)) +
    geom_path() +
    coord_fixed() +
    labs(x = "False Positive Rate", y = "True Positive Rate",
         color = names(data)[n])

  if (!is.null(df$facet)) p <- p + facet_wrap(~ facet)

  p
}


#' @rdname plot-methods
#'
plot.mrmc <- function(x, ...) {
  plot(x$roc)
}
