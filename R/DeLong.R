#' Performance Metric Covariance Estimation
#' 
#' @name cov_methods
#' @rdname cov_methods
#' 
#' @seealso \code{\link{mrmc}}
#' 
DeLong <- function() {
  structure(
    function(formula, data, ...) {
      vars <- extract_vars(formula)
      
      if (vars["metric"] != "roc_auc") {
        stop("response metric must be 'roc_auc' for DeLong covariance method")
      }
      
      if (any(table(data[c(vars[c("tests", "readers")], "(cases)")]) != 1)) {
        stop("balanced design required for DeLong covariance method")
      }
      
      observed <- data[[vars["observed"]]]
      predicted <- data[[vars["predicted"]]]
      groups <- interaction(data[[vars["tests"]]], data[[vars["readers"]]],
                            drop = TRUE)
    
      varcomps <- lapply(levels(groups), function(group) {
        indices <- groups == group
        observed <- observed[indices]
        predicted <- predicted[indices]
        varcomp <- varcomp_Sen(observed, predicted)
        auc <- roc_auc(observed, predicted)
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


varcomp_Sen <- function(observed, predicted) {
  is_pos <- observed == levels(observed)[2]
  predicted_pos <- predicted[is_pos]
  predicted_neg <- predicted[!is_pos]
  
  indices <- expand.grid(pos = seq_along(predicted_pos),
                         neg = seq_along(predicted_neg))
  
  psi_all <- psi(predicted_pos[indices$pos], predicted_neg[indices$neg])
  
  v10 <- tapply(psi_all, indices$pos, mean)
  v01 <- tapply(psi_all, indices$neg, mean)
  
  list(v10 = v10, v01 = v01)
}
