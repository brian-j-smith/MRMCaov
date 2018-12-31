#' Multi-Reader Multi-Case ROC Analysis
#' 
#' Estimation and comparison of ROC performance metrics for multi-reader
#' multi-case studies.
#' 
#' @param formula formula defining the response metric and variables on the left
#' hand side and the test and reader variables on the right hand side.
#' The response should be expressed in terms of a package-supplied performance
#' \code{\link[mrROC:metrics]{metric}}.  Test and reader variables should be specified in that
#' order and separated by the \code{+} operator.
#' @param cases variable containing the case identifiers.
#' @param data data frame containing the \code{formula} and \code{cases}
#' variables.
#' @param method function, function call, or character string naming the
#' \code{\link[mrROC:cov_methods]{method}} to use in calculating performance
#' metric covariances.
#' 
#' @seealso \code{\link{metrics}}, \code{\link{cov_methods}},
#' \code{\link{summary}}, \code{\link{plot}}
#' 
#' @examples
#' (est <- mrmc(roc_auc(truth, rating) ~ treatment + reader, cases = case,
#'              data = VanDyke))
#' summary(est)
#' plot(est)
#' 
mrmc <- function(formula, cases, data, method = DeLong) {
  formula <- update(formula, ~ .^2)
  vars <- extract_vars(formula)
  factor_vars <- c("observed", "tests", "readers")
  data[vars[factor_vars]] <- lapply(data[vars[factor_vars]], factor)
  data[["(cases)"]] <- factor(eval(substitute(cases), data))

  method <- get_method(method)

  fo <- formula
  df_by <- by(data, data[vars[3:4]], function(split) {
    structure(
      c(nrow(split), eval(fo[[2]], split)),
      names = c("N", vars["metric"])
    )
  })
  df <- cbind(expand.grid(dimnames(df_by)), do.call(rbind, df_by))
  fo[[2]] <- as.name(names(df)[ncol(df)])

  structure(
    list(aov = aov(fo, data = df),
         aov_data = df,
         vars = vars,
         cov = method(formula, data),
         mrmc_tests = mrmc_tests(formula, data = data, method = method),
         roc = do.call(roc, c(list(data[[vars["observed"]]],
                                   data[[vars["predicted"]]]),
                              data[vars[c("tests", "readers")]]))),
    class = "mrmc"
  )
}


get_cov_comps <- function(cov, tests, readers) {
  same_test <- outer(tests, tests, "==")
  same_reader <- outer(readers, readers, "==")
  
  c(mean(cov[!same_test & same_reader]),
    mean(cov[same_test & !same_reader]),
    mean(cov[!same_test & !same_reader]))
}


print.mrmc <- function(x, n = 20, ...) {
  model <- x$aov$model
  vars <- x$vars

  cat(vars["tests"], "by", vars["readers"], "Estimates\n\n")
  print(x$aov_data)
  cat("\nANOVA Table\n\n")
  print(summary(x$aov))
  
  cat("\n\nObuchowski-Rockette error variance and covariance estimates\n\n")
  var_comps <- data.frame(
    Estimate = c(mean(diag(x$cov)),
                 get_cov_comps(x$cov,
                               model[[vars["tests"]]],
                               model[[vars["readers"]]])),
    row.names = c("Error", "Cov1", "Cov2", "Cov3")
  )
  var_comps$Correlation = var_comps$Estimate / var_comps$Estimate[1]
  var_comps$Correlation[1] <- NA
  print(var_comps)
  
  cat("\n\nFirst", n, "ROC performance metrics\n\n")
  print(head(x$roc, n = n))
  
  invisible()
}



#' Summary Estimates and Statistical Tests
#' 
#' @name summary
#' @rdname summary-methods
#' 
#' @param object object to summarize.
#' @param design one of the following study designs: 1 = factorial, 2 = cases
#' nested within readers, or 3 = cases nested within tests.
#' @param conf.level confidence level for confidence intervals.
#' @param ... additional arguments affecting the summary.
#' 
#' @seealso \code{\link{mrmc}}
#' 
summary.mrmc <- function(object, design = 1, conf.level = 0.95, ...) {
  MS <- structure(
    summary(object$aov)[[1]][["Mean Sq"]],
    names = c("T", "R", "T:R")
  )
  
  model <- object$aov$model
  vars <- object$vars
  n_tests <- nlevels(model[[vars["tests"]]])
  n_readers <- nlevels(model[[vars["readers"]]])
  var_error <- mean(diag(object$cov))
  cov <- get_cov_comps(object$cov,
                       model[[vars["tests"]]],
                       model[[vars["readers"]]])
  
  if (design == 2) {
    cov[2:3] <- 0
  } else if (design == 3) {
    cov[c(1, 3)] <- 0
  }

  var_comps <- data.frame(
    Estimate = c(
      (MS[["R"]] - MS[["T:R"]]) / n_tests - cov[1] + cov[3],
      MS[["T:R"]] - var_error + cov[1] + (cov[2] - cov[3]),
      var_error,
      cov
    ),
    row.names = c(labels(object$aov$terms)[c(2, 3)], "Error", "Cov1", "Cov2",
                  "Cov3")
  )
  var_comps$Correlation <- var_comps$Estimate / var_comps["Error", "Estimate"]
  var_comps$Correlation[1:3] <- NA

  denominator <- MS[["T:R"]] + n_readers * max(cov[2] - cov[3], 0)
  global_test <- data.frame(
    "MS(T)" = MS[["T"]],
    "MS(T:R)" = MS[["T:R"]],
    "Cov2" = cov[2],
    "Cov3" = cov[3],
    "Denominator" = denominator,
    "F" = MS[["T"]] / denominator,
    "df1" = n_tests - 1,
    "df2" = (MS[["T:R"]] + n_readers * max(cov[2] - cov[3], 0))^2 /
              (MS[["T:R"]]^2 / ((n_tests - 1) * (n_readers - 1))),
    check.names = FALSE
  )
  global_test["p-value"] = 1 - pf(global_test$F, global_test$df1,
                                  global_test$df2)

  metric_means <- c(tapply(model[[1]], model[[vars["tests"]]], mean))
  combs <- combinations(n_tests, 2)
  metric_diffs <- data.frame(
    "Comparison" = paste(combs[, 1], "-", combs[, 2]),
    "Estimate" = metric_means[combs[, 1]] - metric_means[combs[, 2]],
    "StdErr" = sqrt(2 / n_readers * denominator),
    "df" = global_test$df2
  )
  metric_diffs$CI <- metric_diffs$Estimate +
    qt((1 + conf.level) / 2, metric_diffs$df) *
    metric_diffs$StdErr %o% c(Lower = -1, Upper = 1)
  metric_diffs$t <- abs(metric_diffs$Estimate / metric_diffs$StdErr)
  metric_diffs$`p-value` <- 2 * (1 - pt(metric_diffs$t, metric_diffs$df))
  
  metrics_list <- lapply(object$mrmc_tests, summary, design = design,
                         conf.level = conf.level)
  metrics <- do.call(rbind, metrics_list)
  rownames(metrics) <- names(metrics_list)

  structure(
    list(var_comps = var_comps, global_test = global_test,
         metric_diffs = metric_diffs, metrics = metrics,
         design = design, conf.level = conf.level, vars = vars),
    class = "summary.mrmc"
  )
}


print.summary.mrmc <- function(x, ...) {
  cat("Experimental design:",
      switch(x$design,
             "factorial",
             paste("cases nested within", x$vars["readers"]),
             paste("cases nested within", x$vars["tests"])),
      "\n\n")
  
  cat("Obuchowski-Rockette variance component and covariance estimates\n\n")
  print(x$var_comps)
  
  tests_metric <- paste(x$vars["tests"], x$vars["metric"])
  
  cat("\n\nANOVA global test of equal", tests_metric, "\n\n")
  print(x$global_test)
  
  cat("\n\n", 100 * x$conf.level, "% CIs and tests for ", tests_metric,
      " pairwise differences\n\n", sep = "")
  print(x$metric_diffs)
  
  cat("\n\n", 100 * x$conf.level, "% ", tests_metric, " CIs\n\n", sep = "")
  print(x$metrics)
  
  invisible()
}
