#' Performance Plots
#' 
#' @name plot
#' @rdname plot-methods
#' 
#' @param x object to plot.
#' @param ... arguments passed to other methods.
#' 
#' @seealso \code{\link{roc}}
#' 
#' @examples
#' perf <- with(VanDyke, roc(truth, rating, treatment, reader))
#' plot(perf)
#'
plot.roc <- function(x, ...) {
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
plot.mrmc <- function(x, ...) {
  plot(x$roc)
}
