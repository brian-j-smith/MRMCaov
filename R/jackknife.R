#' @rdname cov_methods
#'
#' @references
#' Efron B (1982). The Jackknife, the Bootstrap and Other Resampling Plans.
#' Philadelphia: SIAM.
#'
jackknife <- function() {
  structure(
    function(data, ...) {
      metric_call <- attr(data, "metric_call")
      groups <- interaction(data$test, data$reader)

      pb <- progress::progress_bar$new(
        format = "jackknife [:bar] :percent | :eta",
        total = nlevels(data$case) * nlevels(groups),
        show_after = 1
      )

      select <- c("truth", "rating", "case")
      metrics <- sapply(split(data[select], groups), function(x) {
        case_splits <- split(x[-3], x[[3]])
        lookup <- tibble(
          case = names(case_splits),
          index = match(case_splits, case_splits)
        )
        res <- numeric(nrow(lookup))
        inds <- which(lookup$index == seq_along(res))
        for (ind in inds) {
          res[ind] <- suppressWarnings(
            eval(metric_call, x[x$case != lookup$case[ind], -3])
          )
          pb$tick()
        }
        inds <- -inds
        res[inds] <- res[lookup$index[inds]]
        pb$tick(length(res) - length(inds))
        res
      })

      pb$terminate()

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
