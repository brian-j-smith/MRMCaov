#' @rdname cov_methods
#' 
jackknife <- function() {
  structure(
    function(formula, data, ...) {
      vars <- extract_vars(formula)
      cases <- data[["(cases)"]]
      groups <- interaction(data[[vars["tests"]]], data[[vars["readers"]]])
      df <- data[vars[c("observed", "predicted")]]
      
      metrics <- matrix(NA, nlevels(cases), nlevels(groups))
      for (i in 1:nlevels(cases)) {
        include <- cases != levels(cases)[i]
        metrics[i, ] <- by(df[include, ], groups[include], function(split) {
          eval(formula[[2]], split)
        })
      }
      
      n <- nrow(metrics)
      structure(
        ((n - 1) / n) * crossprod(scale(metrics, scale = FALSE)),
        dimnames = list(levels(groups), levels(groups))
      )
    },
    class = c("cov_method", "function")
  )
}
