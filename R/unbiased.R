#' @rdname cov_methods
#' 
unbiased <- function() {
  structure(
    function(formula, data) {
      vars <- extract_vars(formula)
      
      if (vars["metric"] != "roc_auc") {
        stop("response metric must be 'roc_auc' for unbiased covariance method")
      }
      
      observed <- data[[vars["observed"]]]
      df <- data.frame(
        group = interaction(data[[vars["tests"]]], data[[vars["readers"]]]),
        predicted = data[[vars["predicted"]]],
        case = data[["(cases)"]]
      )

      is_pos <- observed == levels(observed)[2]
      df_pos <- df[is_pos, , drop = FALSE]
      df_neg <- df[!is_pos, , drop = FALSE]
      
      df_list <- lapply(levels(df$group), function(group) {
        .unbiased_psi(df_pos, df_neg, group)
      })

      a1 <- a2 <- a3 <- a4 <- matrix(NA, nlevels(df$group), nlevels(df$group))
      for(i in 1:nlevels(df$group)) {
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
      structure(
        (a1 + (n_pos - 1) * a2 + (n_neg - 1) * a3 + (1 - n_pos - n_neg) * a4) /
          (n_pos * n_neg),
        dimnames = list(levels(df$group), levels(df$group))
      )
    },
    class = c("cov_method", "function")
  )
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

