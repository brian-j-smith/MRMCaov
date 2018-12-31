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
  
  x[c("FPR", "TPR")] <- NULL
  if (ncol(x) > 1) {
    df$group <- interaction(x)
    group_name <- "Group"
  } else {
    df$group <- x[[1]]
    group_name <- names(x)
  }

  aes_perf <- if (nlevels(df$group) > 1) {
    aes(x, y, color = group)
  } else {
    aes(x, y)
  }
  
  ggplot(df, aes_perf) +
    geom_step() +
    coord_fixed() +
    labs(x = "False Positive Rate", y = "True Positive Rate",
         color = group_name)
}


#' @rdname plot-methods
#' 
plot.mrmc <- function(x, ...) {
  plot(x$roc)
}
