#' Performance Metric Covariance Estimation
#'
#' @name cov_methods
#' @rdname cov_methods
#'
#' @seealso \code{\link{mrmc}}
#'
DeLong <- function() {
  structure(
    function(data, ...) {

      metric_call <- attr(data, "metric_call")
      metric_name <- as.character(metric_call)[1]
      if (!(metric_name %in% c("empirical_auc", "trapezoidal_auc"))) {
        stop("response metric must be 'empirical_auc' or 'trapezoidal_auc' for",
             " DeLong covariance method")
      }

      partial <- as.list(metric_call)$partial
      if (!(is.null(partial) || isFALSE(partial))) {
        stop("DeLong covariance method not available for partial AUC")
      }

      if (any(table(data[c("test", "reader", "case")]) != 1)) {
        stop("balanced design required for DeLong covariance method")
      }

      truths <- data$truth
      ratings <- data$rating
      groups <- interaction(data$test, data$reader, drop = TRUE)

      varcomps <- lapply(levels(groups), function(group) {
        indices <- groups == group
        truths <- truths[indices]
        ratings <- ratings[indices]
        varcomp <- varcomp_Sen(truths, ratings)
        auc <- empirical_auc(truths, ratings)
        list(varcomp10 = varcomp$v10 - auc, varcomp01 = varcomp$v01 - auc,
             auc = auc)
      })

      varcomp10_mat <- sapply(varcomps, getElement, name = "varcomp10")
      n_pos <- nrow(varcomp10_mat)
      s10_mat <- crossprod(varcomp10_mat) / (n_pos - 1)

      varcomp01_mat <- sapply(varcomps, getElement, name = "varcomp01")
      n_neg <- nrow(varcomp01_mat)
      s01_mat <- crossprod(varcomp01_mat) / (n_neg - 1)

      structure(
        s10_mat / n_pos + s01_mat / n_neg,
        dimnames = list(levels(groups), levels(groups))
      )
    },
    class = c("cov_method", "function")
  )
}


varcomp_Sen <- function(truths, ratings) {
  is_pos <- truths == levels(truths)[2]
  ratings_pos <- ratings[is_pos]
  ratings_neg <- ratings[!is_pos]

  indices <- expand.grid(pos = seq_along(ratings_pos),
                         neg = seq_along(ratings_neg))

  psi_all <- psi(ratings_pos[indices$pos], ratings_neg[indices$neg])

  v10 <- tapply(psi_all, indices$pos, mean)
  v01 <- tapply(psi_all, indices$neg, mean)

  list(v10 = v10, v01 = v01)
}
