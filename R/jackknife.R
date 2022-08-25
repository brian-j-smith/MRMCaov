#' @rdname cov_methods
#'
#' @references
#' Efron B (1982). The Jackknife, the Bootstrap and Other Resampling Plans.
#' Philadelphia: SIAM.
#'
jackknife <- function() {
  structure(
    function(data, ...) {
      cases <- data$case
      groups <- interaction(data$test, data$reader)
      df <- data[c("truth", "rating")]

      metrics <- matrix(NA, nlevels(cases), nlevels(groups))
      for (i in 1:nlevels(cases)) {
        keep <- cases != levels(cases)[i]
        metrics[i, ] <- by(df[keep, ], groups[keep], function(split) {
          eval(attr(data, "metric_call"), split)
        })
      }

      n <- nrow(metrics)
      structure(
        ((n - 1) / n) * crossprod(scale(metrics, scale = FALSE)),
        dimnames = list(levels(groups), levels(groups)),
        class = c("cov_jackknife", "cov_matrix")
      )
    },
    class = c("cov_method", "function")
  )
}
