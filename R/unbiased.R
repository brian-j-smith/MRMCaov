#' @rdname cov_methods
#'
#' @param abar logical indicating whether to compute mean covariance components.
#'
unbiased <- function(abar = FALSE) {
  structure(
    function(data, ...) {

      metric_call <- attr(data, "metric_call")
      metric_name <- as.character(metric_call)[1]
      if (!(metric_name %in% c("empirical_auc", "trapezoidal_auc"))) {
        stop("response metric must be 'empirical_auc' or 'trapezoidal_auc' for",
             " for unbiased covariance method")
      }

      partial <- as.list(metric_call)$partial
      if (!(is.null(partial) || isFALSE(partial))) {
        stop("unbiased covariance method not available for partial AUC")
      }

      data$group <- interaction(data$test, data$reader)

      f <- if (is_fully_paired(data)) unbiased_balanced else unbiased_default
      covmat <- f(data$truth, data$rating, data$group, data$case)
      dimnames(covmat) <- list(levels(data$group), levels(data$group))

      a <- attr(covmat, "a")
      class(covmat) <- c("cov_unbiased", "cov_matrix")

      if (abar) {
        group_levels <- expand.grid(test = levels(data$test),
                                    reader = levels(data$reader))
        same_test <- outer(group_levels$test, group_levels$test, "==")
        same_reader <- outer(group_levels$reader, group_levels$reader, "==")
        abar_sigma2 <- sapply(a, function(x) mean(x[same_test & same_reader]))
        abar_cov1 <- sapply(a, function(x) mean(x[!same_test & same_reader]))
        abar_cov2 <- sapply(a, function(x) mean(x[same_test & !same_reader]))
        abar_cov3 <- sapply(a, function(x) mean(x[!(same_test | same_reader)]))
        attr(covmat, "abar") <- cbind(abar_sigma2, abar_cov1, abar_cov2,
                                      abar_cov3)
      }

      covmat

    },
    class = c("cov_method", "function")
  )
}


unbiased_default <- function(truth, rating, group, case) {

  n <- nlevels(group)
  is_pos <- truth == levels(truth)[2]

  V <- matrix(0, n, n)
  x <- list()
  for (i in 1:n) {
    x[[i]] <- get_scores(rating, group == levels(group)[i], is_pos, case)
    for (j in 1:i) {
      V[i, j] <- get_cov(x[[i]], x[[j]])
      V[j, i] <- V[i, j]
    }
  }
  V

}


get_scores <- function(rating, is_group, is_pos, case) {
  res <- list()
  is_group_pos <- is_group & is_pos
  is_group_neg <- is_group & !is_pos
  res$id_pos <- case[is_group_pos]
  res$id_neg <- case[is_group_neg]
  res$scores <- outer(rating[is_group_pos], rating[is_group_neg], psi)
  res$sum_scores_pos <- rowSums(res$scores)
  res$sum_scores_neg <- colSums(res$scores)
  res$sum_scores <- sum(res$sum_scores_neg)
  res
}


get_cov <- function(x, y) {

  id_pos <- intersect(x$id_pos, y$id_pos)
  id_neg <- intersect(x$id_neg, y$id_neg)
  nxny <- length(x$scores) * length(y$scores)
  n_pos <- length(id_pos)
  n_neg <- length(id_neg)

  delta <- nxny / (
    nxny -
    nrow(x$scores) * n_neg * nrow(y$scores) -
    ncol(x$scores) * n_pos * ncol(y$scores) +
    n_pos * n_neg
  )

  inds_pos <- match(id_pos, x$id_pos)
  inds_neg <- match(id_neg, x$id_neg)
  x$sum_scores_pos <- x$sum_scores_pos[inds_pos]
  x$sum_scores_neg <- x$sum_scores_neg[inds_neg]
  x$scores <- x$scores[inds_pos, inds_neg]

  inds_pos <- match(id_pos, y$id_pos)
  inds_neg <- match(id_neg, y$id_neg)
  y$sum_scores_pos <- y$sum_scores_pos[inds_pos]
  y$sum_scores_neg <- y$sum_scores_neg[inds_neg]
  y$scores <- y$scores[inds_pos, inds_neg]

  (
    (1 - delta) * (x$sum_scores * y$sum_scores) +
    delta * (sum(x$sum_scores_pos * y$sum_scores_pos) +
             sum(x$sum_scores_neg * y$sum_scores_neg) -
             sum(x$scores * y$scores))
  ) / nxny

}


unbiased_balanced <- function(truth, rating, group, case) {

  sort_inds <- order(group, case)

  m <- nlevels(case)
  n <- nlevels(group)
  rating_mat <- matrix(rating[sort_inds], m, n)
  truth <- truth[sort_inds[1:m]]

  is_pos <- truth == levels(truth)[2]
  rating_mat_pos <- rating_mat[is_pos, ]
  rating_mat_neg <- rating_mat[!is_pos, ]
  n_pos <- nrow(rating_mat_pos)
  n_neg <- nrow(rating_mat_neg)

  scores <- matrix(0, n_pos * n_neg, n)
  sum_scores_pos <- matrix(0, n_pos, n)
  sum_scores_neg <- matrix(0, n_neg, n)

  for (i in 1:n) {
    x <- outer(rating_mat_pos[, i], rating_mat_neg[, i], psi)
    scores[, i] <- c(x)
    sum_scores_pos[, i] <- rowSums(x)
    sum_scores_neg[, i] <- colSums(x)
  }

  delta <- (n_neg * n_pos)^2 / (n_neg * (n_neg - 1) * n_pos * (n_pos - 1))
  sum_scores <- rbind(colSums(scores))
  (
    (1 - delta) * (t(sum_scores) %*% sum_scores) +
    delta * (t(sum_scores_pos) %*% sum_scores_pos +
             t(sum_scores_neg) %*% sum_scores_neg -
             t(scores) %*% scores)
  ) / (n_neg * n_pos)^2

}
