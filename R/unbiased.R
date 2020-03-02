#' @rdname cov_methods
#'
unbiased <- function() {
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

      df <- data[c("truth", "rating", "case")]
      df$group <- interaction(data$test, data$reader)

      tbl <- table(df$group, df$case)
      balanced <- all(tbl == tbl[1])

      covmat <- ifelse(balanced, .unbiased_balanced, .unbiased_unbalanced)(df)
      dimnames(covmat) <- list(levels(df$group), levels(df$group))

      a <- attr(covmat, "a")
      group_levels <- expand.grid(test = levels(data$test),
                                  reader = levels(data$reader))
      same_test <- outer(group_levels$test, group_levels$test, "==")
      same_reader <- outer(group_levels$reader, group_levels$reader, "==")
      abar_sigma2 <- sapply(a, function(x) mean(x[same_test & same_reader]))
      abar_cov1 <- sapply(a, function(x) mean(x[!same_test & same_reader]))
      abar_cov2 <- sapply(a, function(x) mean(x[same_test & !same_reader]))
      abar_cov3 <- sapply(a, function(x) mean(x[!(same_test | same_reader)]))

      structure(covmat,
                a = NULL,
                abar = cbind(abar_sigma2, abar_cov1, abar_cov2, abar_cov3),
                class = "cov_unbiased")

    },
    class = c("cov_method", "function")
  )
}


.unbiased_balanced <- function(data) {

  df <- reshape(data, idvar = "case", v.names = "rating",
                timevar = "group", direction = "wide")
  subset_indices <- -match(c("case", "truth"), names(df))

  is_pos <- df$truth == levels(df$truth)[2]
  df_pos <- df[is_pos, subset_indices, drop = FALSE]
  df_neg <- df[!is_pos, subset_indices, drop = FALSE]

  n_pos <- nrow(df_pos)
  n_neg <- nrow(df_neg)

  indices <- expand.grid(pos = seq_len(n_pos), neg = seq_len(n_neg))

  same_case_pos <- outer(indices$pos, indices$pos, "==")
  same_case_neg <- outer(indices$neg, indices$neg, "==")

  in_a1 <- same_case_pos & same_case_neg
  in_a2 <- same_case_pos & !same_case_neg
  in_a3 <- !same_case_pos & same_case_neg
  in_a4 <- !same_case_pos & !same_case_neg

  psi_list <- lapply(1:nlevels(data$group), function(i) {
    psi(df_pos[indices$pos, i], df_neg[indices$neg, i])
  })

  a1 <- a2 <- a3 <- a4 <- matrix(NA, nlevels(data$group), nlevels(data$group))
  for (i in 1:nlevels(data$group)) {
    psi_i <- psi_list[[i]]
    for (j in 1:i) {
      psi_j <- psi_list[[j]]

      psi_cross <- tcrossprod(psi_i, psi_j)

      a1[i, j] <- a1[j, i] <- mean(psi_cross[in_a1])
      a2[i, j] <- a2[j, i] <- mean(psi_cross[in_a2])
      a3[i, j] <- a3[j, i] <- mean(psi_cross[in_a3])
      a4[i, j] <- a4[j, i] <- mean(psi_cross[in_a4])

    }
  }

  covmat <- (a1 + (n_neg - 1) * a2 + (n_pos - 1) * a3 +
               (1 - n_pos - n_neg) * a4) / (n_pos * n_neg)

  structure(covmat, a = list(a1, a2, a3, a4))

}


.unbiased_unbalanced <- function(data) {

  subset_indices <- -match("truth", names(data))

  is_pos <- data$truth == levels(data$truth)[2]
  df_pos <- data[is_pos, subset_indices, drop = FALSE]
  df_neg <- data[!is_pos, subset_indices, drop = FALSE]

  df_list <- lapply(levels(data$group), function(group) {
    .unbiased_psi(df_pos, df_neg, group)
  })

  a1 <- a2 <- a3 <- a4 <- matrix(NA, nlevels(data$group), nlevels(data$group))
  for(i in 1:nlevels(data$group)) {
    df_i <- df_list[[i]]
    for(j in 1:i) {
      df_j <- df_list[[j]]

      psi_cross <- tcrossprod(df_i$psi, df_j$psi)
      same_case_pos <- outer(df_i$case_pos, df_j$case_pos, "==")
      same_case_neg <- outer(df_i$case_neg, df_j$case_neg, "==")

      mean_of <- function(include) mean(psi_cross[include])
      a1[i, j] <- a1[j, i] <- mean_of(same_case_pos & same_case_neg)
      a2[i, j] <- a2[j, i] <- mean_of(same_case_pos & !same_case_neg)
      a3[i, j] <- a3[j, i] <- mean_of(!same_case_pos & same_case_neg)
      a4[i, j] <- a4[j, i] <- mean_of(!same_case_pos & !same_case_neg)
    }
  }

  n_pos <- sum(!duplicated(df_pos$case))
  n_neg <- sum(!duplicated(df_neg$case))
  covmat <- (a1 + (n_neg - 1) * a2 + (n_pos - 1) * a3 +
               (1 - n_pos - n_neg) * a4) / (n_pos * n_neg)

  structure(covmat, a = list(a1, a2, a3, a4))

}


.unbiased_psi <- function(df_pos, df_neg, group) {

  is_group_pos <- df_pos$group == group
  ratings_pos <- df_pos$rating[is_group_pos]
  cases_pos <- df_pos$case[is_group_pos]

  is_group_neg <- df_neg$group == group
  ratings_neg <- df_neg$rating[is_group_neg]
  cases_neg <- df_neg$case[is_group_neg]

  indices <- expand.grid(pos = seq_along(ratings_pos),
                         neg = seq_along(ratings_neg))

  data.frame(psi = psi(ratings_pos[indices$pos], ratings_neg[indices$neg]),
             case_pos = cases_pos[indices$pos],
             case_neg = cases_neg[indices$neg])

}
