aov_mrmc <- function(data) {
  fo <- update(formula(data), . ~ .^2)
  aov(fo, data = data)
}


chol2det <- function(x, log = FALSE) {
  d <- diag(x)
  if (log) 2 * sum(log(d)) else prod(d)^2
}


curves2tibble <- function(x, groups) {
  groups_ind <- rep(seq(nrow(groups)), times = sapply(x, nrow))
  x <- do.call(rbind, x)
  tibble(Group = groups[groups_ind, , drop = FALSE], FPR = x$FPR, TPR = x$TPR)
}


dim.mrmc <- function(x) {
  sapply(levels(x), length)
}


dim.mrmc_tests <- function(x) {
  sapply(levels(x), length)
}


get_design <- function(data) {
  crosstab <- function(...) table(data[c(...)]) > 0

  test_x_reader <- crosstab("test", "reader") > 0
  if (all(test_x_reader)) {
    if (all(colSums(crosstab("reader", "case")) == 1)) {
      2
    } else if (all(colSums(crosstab("test", "case")) == 1)) {
      3
    } else {
      1
    }
  } else if (all(colSums(test_x_reader) == 1)) {
    4
  }
}


get_method <- function(x) {
  if (is(x, "character")) x <- get(x, mode = "function")
  if (!is(x, "cov_method")) x <- x()
  if (!is(x, "cov_method")) stop("invalid covariance method")
  x
}


is_balanced <- function(data) {
  all(table(data[c("test", "reader", "case")]) == 1)
}


levels.mrmc <- function(x) {
  structure(x$aov$xlevels, names = c("test", "reader"))
}


levels.mrmc_tests <- function(x) {
  structure(x$aov$xlevels, names = c("reader"))
}


mean.mrmc_tests <- function(x, ...) {
  mean(x$aov$model[[1]])
}


meansq <- function(x, ...) {
  UseMethod("meansq")
}


meansq.mrmc <- function(x, ...) {
  res <- summary(x$aov)[[1]][["Mean Sq"]]
  structure(res, names = head(c("T", "R", "T:R"), length(res)))
}


meansq.mrmc_tests <- function(x, ...) {
  c("T" = 0, "R" = summary(x$aov)[[1]][["Mean Sq"]], "T:R" = 0)
}


preprocess <- function(data) {
  metric_name <- as.character(attr(data, "metric_call"))[1]
  level <- switch(metric_name,
                  "binary_sens" = 2,
                  "binary_spec" = 1)
  if (!is.null(level)) {
    keep <- data$truth == levels(data$truth)[level]
    droplevels(data[keep, , drop = FALSE],
               except = which(names(data) == "truth"))
  } else data
}


trunc_ci <- function(object, x, ...) {
  UseMethod("trunc_ci")
}


trunc_ci.character <- function(object, x, ...) {
  unit_metrics <- c("_auc", "_sens", "_spec")
  if (any(endsWith(object, unit_metrics))) {
    x[, 1] <- pmax(0, x[, 1])
    x[, 2] <- pmin(x[, 2], 1)
  }
  x
}


trunc_ci.mrmc <- function(object, x, ...) {
  metric <- as.character(attr(object$mrmc_data, "metric_call"))[1]
  trunc_ci(metric, x)
}


trunc_ci.srmc <- function(object, x, ...) {
  c(trunc_ci(object$metric, rbind(x)))
}


vcov_comps <- function(object, ...) {
  UseMethod("vcov_comps")
}


vcov_comps.mrmc <- function(object, design = object$design, test = NULL,
                            reader = NULL, ...) {
  data <- object$data

  tests <- data[[2]]
  readers <- data[[3]]

  same_test <- outer(tests, tests, "==")
  same_reader <- outer(readers, readers, "==")

  is_group <- rep(TRUE, nrow(data))
  if (!is.null(test)) is_group <- is_group & (tests == test)
  if (!is.null(reader)) is_group <- is_group & (readers == reader)
  in_group <- as.logical(is_group %o% is_group)

  cov <- c(mean(object$cov[!same_test & same_reader & in_group]),
           mean(object$cov[same_test & !same_reader & in_group]),
           mean(object$cov[!same_test & !same_reader & in_group]))

  if (design == 2) {
    cov[2:3] <- 0
  } else if (design == 3) {
    cov[c(1, 3)] <- 0
  } else if (design == 4) {
    cov[1] <- 0
  }

  n_mat <- table(data[2:3])
  names(dimnames(n_mat)) <- c("test", "reader")

  structure(
    list(vars = object$vars,
         n = dim(object),
         n_mat = n_mat,
         MS = meansq(object),
         var = mean(diag(object$cov[is_group, is_group])),
         cov = cov),
    class ="vcov_comps"
  )
}


vcov_comps.mrmc_tests <- function(object, design = object$design, ...) {
  data <- object$aov$model

  cov2 <- if (design == 2) 0 else {
    readers <- data[[2]]
    same_reader <- outer(readers, readers, "==")
    mean(object$cov[!same_reader])
  }

  structure(
    list(vars = object$vars,
         n = dim(object),
         MS = meansq(object),
         var = mean(diag(object$cov)),
         cov = c(0, cov2, 0)),
    class = "vcov_comps"
  )

}


summary.vcov_comps <- function(object, ...) {
  n <- object$n
  MS <- object$MS
  var_error <- object$var
  cov <- object$cov

  if ("T:R" %in% names(MS)) {
    est <- c((MS[["R"]] - MS[["T:R"]]) / n[["test"]] - cov[1] + cov[3],
             MS[["T:R"]] - var_error + cov[1] + (cov[2] - cov[3]))
    names(est) <- c(object$vars["reader"],
                    paste0(object$vars[c("test", "reader")], collapse = ":"))
  } else {
    cat(str(object))
    est <- MS[["R"]] - var_error + cov[2]
    names(est) <- object$vars["reader"]
  }
  est["Error"] <- var_error
  est[paste0("Cov", 1:3)] <- cov

  vcov_comps <- data.frame(Estimate = est)
  vcov_comps$Correlation <- vcov_comps$Estimate /
    vcov_comps["Error", "Estimate"]
  vcov_comps$Correlation[1:3] <- NA

  vcov_comps
}
