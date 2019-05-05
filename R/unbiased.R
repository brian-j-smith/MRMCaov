#' @rdname cov_methods
#' 
unbiased <- function() {
  structure(
    function(formula, data, ...) {
      vars <- extract_vars(formula)
      
      if (vars["metric"] != "roc_auc") {
        stop("response metric must be 'roc_auc' for unbiased covariance method")
      }
      
      df <- data.frame(
        observed = data[[vars["observed"]]],
        case = data[["(cases)"]],
        predicted = data[[vars["predicted"]]],
        group = interaction(data[[vars["tests"]]], data[[vars["readers"]]],
                            drop = TRUE)
      )
      
      tbl <- table(df$group, df$case)
      balanced <- all(tbl == tbl[1])

      structure(
        ifelse(balanced, .unbiased_balanced, .unbiased_unbalanced)(df),
        dimnames = list(levels(data$group), levels(data$group))
      )
    },
    class = c("cov_method", "function")
  )
}


.unbiased_balanced <- function(data) {
  
  df <- reshape(data, idvar = "case", v.names = "predicted",
                timevar = "group", direction = "wide")
  subset_indices <- -match(c("case", "observed"), names(df))
  
  is_pos <- df$observed == levels(df$observed)[2]
  df_pos <- df[is_pos, subset_indices, drop = FALSE]
  df_neg <- df[!is_pos, subset_indices, drop = FALSE]

  n_pos <- nrow(df_pos)
  n_neg <- nrow(df_neg)
  
  indices <- expand.grid(pos = seq_len(n_pos), neg = seq_len(n_neg))
  
  same_case_pos <- outer(indices$pos, indices$pos, "==")
  same_case_neg <- outer(indices$neg, indices$neg, "==")
  
  in_a1 <- same_case_pos & same_case_neg
  in_a2 <- !same_case_pos & same_case_neg
  in_a3 <- same_case_pos & !same_case_neg
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
  
  (a1 + (n_pos - 1) * a2 + (n_neg - 1) * a3 + (1 - n_pos - n_neg) * a4) /
    (n_pos * n_neg)
  
}


.unbiased_unbalanced <- function(data) {

  subset_indices <- -match("observed", names(data))

  is_pos <- data$observed == levels(data$observed)[2]
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
      a2[i, j] <- a2[j, i] <- mean_of(!same_case_pos & same_case_neg)
      a3[i, j] <- a3[j, i] <- mean_of(same_case_pos & !same_case_neg)
      a4[i, j] <- a4[j, i] <- mean_of(!same_case_pos & !same_case_neg)
    }
  }

  n_pos <- sum(!duplicated(df_pos$case))
  n_neg <- sum(!duplicated(df_neg$case))
  (a1 + (n_pos - 1) * a2 + (n_neg - 1) * a3 + (1 - n_pos - n_neg) * a4) /
    (n_pos * n_neg)
  
}


.unbiased_psi <- function(df_pos, df_neg, group) {
  
  is_group_pos <- df_pos$group == group
  predicted_pos <- df_pos$predicted[is_group_pos]
  cases_pos <- df_pos$case[is_group_pos]
  
  is_group_neg <- df_neg$group == group
  predicted_neg <- df_neg$predicted[is_group_neg]
  cases_neg <- df_neg$case[is_group_neg]
  
  indices <- expand.grid(pos = seq_along(predicted_pos),
                         neg = seq_along(predicted_neg))
  
  data.frame(psi = psi(predicted_pos[indices$pos], predicted_neg[indices$neg]),
             case_pos = cases_pos[indices$pos],
             case_neg = cases_neg[indices$neg])
  
}