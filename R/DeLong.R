#' Covariance Methods
#'
#' Reader performance metric covariance estimation methods to be used with
#' package-supplied multi-case statistical analysis functions.
#'
#' @name cov_methods
#' @rdname cov_methods
#'
#' @return
#' Returns a function of class \code{cov_method} specifying a covariance
#' method for \code{\link{mrmc}}, \code{\link{srmc}}, or \code{\link{stmc}}.
#'
#' @seealso \code{\link{mrmc}}, \code{\link{srmc}}, \code{\link{stmc}}
#'
#' @references
#' DeLong ER, DeLong DM, and Clarke-Pearson DL (1988). Comparing the areas under
#' two or more correlated receiver operating characteristic curves: a
#' nonparametric approach. Biometrics, 44: 837–45.
#'
DeLong <- function() {
  structure(
    function(data, ...) {

      metric_call <- attr(data, "metric_call")
      metric_name <- as.character(metric_call)[1]
      if (!(metric_name %in% c("empirical_auc", "trapezoidal_auc"))) {
        stop(
          "response metric must be 'empirical_auc' or 'trapezoidal_auc' for",
          " DeLong covariance method"
        )
      }

      partial <- as.list(metric_call)$partial
      if (!(is.null(partial) || isFALSE(partial))) {
        stop("DeLong covariance method not available for partial AUC")
      }

      if (!is_fully_paired(data)) {
        stop("balanced design required for DeLong covariance method")
      }

      truth <- data$truth
      rating <- data$rating
      group <- interaction(data$test, data$reader)

      varcomps <- lapply(levels(group), function(level) {
        keep <- group == level
        truth <- truth[keep]
        rating <- rating[keep]
        varcomp <- varcomp_Sen(truth, rating)
        auc <- empirical_auc(truth, rating)
        list(
          varcomp10 = varcomp$v10 - auc, varcomp01 = varcomp$v01 - auc,
          auc = auc
        )
      })

      varcomp10_mat <- sapply(varcomps, getElement, name = "varcomp10")
      n_pos <- nrow(varcomp10_mat)
      s10_mat <- crossprod(varcomp10_mat) / (n_pos - 1)

      varcomp01_mat <- sapply(varcomps, getElement, name = "varcomp01")
      n_neg <- nrow(varcomp01_mat)
      s01_mat <- crossprod(varcomp01_mat) / (n_neg - 1)

      structure(
        s10_mat / n_pos + s01_mat / n_neg,
        dimnames = list(levels(group), levels(group)),
        class = c("cov_DeLong", "cov_matrix")
      )
    },
    class = c("cov_method", "function")
  )
}


varcomp_Sen <- function(truth, rating) {
  is_pos <- is_reference(truth)
  rating_pos <- rating[is_pos]
  rating_neg <- rating[!is_pos]

  indices <- expand.grid(
    pos = seq_along(rating_pos), neg = seq_along(rating_neg)
  )

  psi_all <- psi(rating_pos[indices$pos], rating_neg[indices$neg])

  v10 <- tapply(psi_all, indices$pos, mean)
  v01 <- tapply(psi_all, indices$neg, mean)

  list(v10 = v10, v01 = v01)
}
