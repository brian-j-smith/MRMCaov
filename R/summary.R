#' Summary Estimates and Statistical Tests
#'
#' @name summary
#' @rdname summary-methods
#'
#' @param object object to summarize.
#' @param conf.level confidence level for confidence intervals.
#' @param ... additional arguments affecting the summary.
#'
#' @seealso \code{\link{mrmc}}
#'
summary.mrmc <- function(object, conf.level = 0.95, ...) {
  f <- ifelse(object$design == 4, .summary_nested, .summary)
  f(object, conf.level = conf.level, ...)
}


#' @rdname summary-methods
#'
summary.srmc <- function(object, conf.level = 0.95, ...) {
  z <- qnorm((1 + conf.level) / 2)
  structure(
    c(object$est, object$se, object$est + c(-1, 1) * z * object$se),
    names = c(object$metric, "StdErr", "CI.Lower", "CI.Upper")
  )
}


new_summary_mrmc <- function(object, conf.level, vcov_comps,
                             test_equality = NULL, test_diffs = NULL,
                             test_means = NULL, reader_means = NULL) {
  cov_method <- class(object$cov)[1]
  structure(
    list(data_name = as.character(object$call$data),
         cov_method = substring(cov_method,
                                4 * startsWith(cov_method, "cov_") + 1),
         design = object$design,
         vars = object$vars,
         conf.level = conf.level,
         vcov_comps = vcov_comps,
         test_equality = test_equality,
         test_diffs = test_diffs,
         test_means = test_means,
         reader_means = reader_means),
    class = "summary.mrmc"
  )
}


.summary <- function(object, ...) {
  UseMethod(".summary")
}

.summary_nested <- function(object, ...) {
  UseMethod(".summary_nested")
}


.summary.mrmc_rrrc <- function(object, conf.level, ...) {
  comps <- vcov_comps(object)
  n <- comps$n
  MS <- comps$MS
  cov <- comps$cov

  test_levels <- levels(object)$test

  denominator <- MS[["T:R"]] + n[["reader"]] * max(cov[2] - cov[3], 0)
  test_equality <- data.frame(
    `MS(T)` = MS[["T"]],
    `MS(T:R)` = MS[["T:R"]],
    Cov2 = cov[2],
    Cov3 = cov[3],
    Denominator = denominator,
    F = MS[["T"]] / denominator,
    df1 = n[["test"]] - 1,
    df2 = (MS[["T:R"]] + n[["reader"]] * max(cov[2] - cov[3], 0))^2 /
      (MS[["T:R"]]^2 / ((n[["test"]] - 1) * (n[["reader"]] - 1))),
    check.names = FALSE
  )
  test_equality$`p-value` <- with(test_equality, 1 - pf(F, df1, df2))

  test_means_list <- lapply(object$mrmc_tests, summary.mrmc_tests_rrrc,
                            conf.level = conf.level)
  test_means <- do.call(rbind, test_means_list)
  rownames(test_means) <- names(test_means_list)

  estimates <- test_means$Estimate
  combs <- combinations(length(estimates), 2)
  test_diffs <- data.frame(
    Comparison = paste(test_levels[combs[, 1]], "-", test_levels[combs[, 2]]),
    Estimate = estimates[combs[, 1]] - estimates[combs[, 2]],
    StdErr = sqrt(2 / n[["reader"]] * denominator),
    df = test_equality$df2
  )
  test_diffs$CI <- with(test_diffs, {
    Estimate + qt((1 + conf.level) / 2, df) * StdErr %o% c(Lower = -1, Upper = 1)
  })
  test_diffs$t <- with(test_diffs, Estimate / StdErr)
  test_diffs$`p-value` <- with(test_diffs, 2 * (1 - pt(abs(t), df)))

  res <- new_summary_mrmc(object,
                          conf.level = conf.level,
                          vcov_comps = summary(comps),
                          test_equality = test_equality,
                          test_diffs = test_diffs,
                          test_means = test_means)
  structure(res, class = c("summary.mrmc_rrrc", class(res)))
}


.summary_nested.mrmc_rrrc <- function(object, conf.level, ...) {
  comps <- vcov_comps(object)
  n_test <- nrow(comps$n_mat)
  n_readers <- rowSums(comps$n_mat)
  n_reader <- sum(n_readers)
  MS <- comps$MS
  cov <- comps$cov

  test_levels <- levels(object)$test

  n_formula <- (n_reader - sum(n_readers^2)) / (n_reader * (n_test - 1))
  denominator <- MS[["R"]] + n_formula * max(cov[2] - cov[3], 0)
  test_equality <- data.frame(
    `MS[T]` = MS[["T"]],
    `MS[R(T)]` = MS[["R"]],
    Cov2 = cov[2],
    Cov3 = cov[3],
    Denominator = denominator,
    F = MS[["T"]] / denominator,
    df1 = n_test - 1,
    df2 = (MS[["R"]] + n_formula * max(cov[2] - cov[3], 0))^2 /
      (MS[["R"]]^2 / (n_reader - n_test)),
    check.names = FALSE
  )
  test_equality$`p-value` <- with(test_equality, 1 - pf(F, df1, df2))

  test_means_list <- lapply(object$mrmc_tests, summary.mrmc_tests_rrrc,
                            conf.level = conf.level)
  test_means <- do.call(rbind, test_means_list)
  rownames(test_means) <- names(test_means_list)

  estimates <- test_means$Estimate
  combs <- combinations(length(estimates), 2)
  test_diffs <- data.frame(
    Comparison = paste(test_levels[combs[, 1]], "-", test_levels[combs[, 2]]),
    Estimate = estimates[combs[, 1]] - estimates[combs[, 2]],
    StdErr = sqrt((1 / n_readers[combs[, 1]] + 1 / n_readers[combs[, 2]]) *
                    denominator),
    df = test_equality$df2,
    row.names = NULL
  )
  test_diffs$CI <- with(test_diffs, {
    Estimate + qt((1 + conf.level) / 2, df) * StdErr %o% c(Lower = -1, Upper = 1)
  })
  test_diffs$t <- with(test_diffs, Estimate / StdErr)
  test_diffs$`p-value` <- with(test_diffs, 2 * (1 - pt(abs(t), df)))

  res <- new_summary_mrmc(object,
                          conf.level = conf.level,
                          vcov_comps = summary(comps),
                          test_equality = test_equality,
                          test_diffs = test_diffs,
                          test_means = test_means)
  structure(res, class = c("summary.mrmc_rrrc", class(res)))
}


.summary.mrmc_frrc <- function(object, conf.level, ...) {
  comps <- vcov_comps(object)
  n <- comps$n
  MS <- comps$MS
  cov <- comps$cov

  test_levels <- levels(object)$test
  reader_levels <- levels(object)$reader

  denominator <- comps$var - cov[1] + (n[["reader"]] - 1) * (cov[2] - cov[3])
  test_equality <- data.frame(
    `MS(T)` = MS[["T"]],
    Cov1 = cov[1],
    Cov2 = cov[2],
    Cov3 = cov[3],
    Denominator = denominator,
    X2 = (n[["test"]] - 1) * MS[["T"]] / denominator,
    df = n[["test"]] - 1,
    check.names = FALSE
  )
  test_equality$`p-value` <- with(test_equality, 1 - pchisq(X2, df))

  test_means_list <- lapply(object$mrmc_tests, summary.mrmc_tests_frrc,
                            conf.level = conf.level)
  test_means <- do.call(rbind, test_means_list)
  rownames(test_means) <- names(test_means_list)

  estimates <- test_means$Estimate
  combs <- combinations(length(estimates), 2)
  test_diffs <- data.frame(
    Comparison = paste(test_levels[combs[, 1]], "-", test_levels[combs[, 2]]),
    Estimate = estimates[combs[, 1]] - estimates[combs[, 2]],
    StdErr = sqrt(2 / n[["reader"]] * denominator)
  )
  test_diffs$CI <- with(test_diffs, {
    Estimate + qnorm((1 + conf.level) / 2) * StdErr %o% c(Lower = -1, Upper = 1)
  })
  test_diffs$z <- with(test_diffs, Estimate / StdErr)
  test_diffs$`p-value` <- with(test_diffs, 2 * (1 - pnorm(abs(z))))


  reader_means <- object$data
  reader_means$StdErr <- sqrt(diag(object$cov))
  reader_means$CI <- reader_means[[3]] + qnorm((1 + conf.level) / 2) *
    reader_means$StdErr %o% c(Lower = -1, Upper = 1)

  res <- new_summary_mrmc(object,
                          conf.level = conf.level,
                          vcov_comps = summary(comps),
                          test_equality = test_equality,
                          test_diffs = test_diffs,
                          test_means = test_means,
                          reader_means = reader_means)
  res$reader_test_diffs = reader_test_diffs(object, conf.level)
  structure(res, class = c("summary.mrmc_frrc", class(res)))
}


reader_test_diffs <- function(object, conf.level) {

  n <- dim(object)
  test_levels <- levels(object)$test
  reader_levels <- levels(object)$reader

  estimates <- matrix(object$aov$model[[1]], ncol = n["test"], byrow = TRUE)

  stderrs <- sapply(reader_levels, function(reader) {
    comps <- vcov_comps(object, reader = reader)
    sqrt(2 * (comps$var - comps$cov[1]))
  })

  combs <- combinations(n["test"], 2)
  combs_reader <- rep(1, n["reader"]) %x% combs
  reader_indices <- rep(seq(reader_levels), each = nrow(combs))

  df <- data.frame(
    Reader = reader_levels[reader_indices],
    Comparison = paste(test_levels[combs_reader[, 1]],
                       "-",
                       test_levels[combs_reader[, 2]]),
    Estimate = estimates[cbind(reader_indices, combs_reader[, 1])] -
      estimates[cbind(reader_indices, combs_reader[, 2])],
    StdErr = stderrs
  )
  df$CI <- with(df, {
    Estimate + qnorm((1 + conf.level) / 2) * StdErr %o% c(Lower = -1, Upper = 1)
  })
  df$z <- with(df, Estimate / StdErr)
  df$`p-value` <- with(df, 2 * (1 - pnorm(abs(z))))

  df
}


.summary_nested.mrmc_frrc <- function(object, conf.level, ...) {
  comps <- vcov_comps(object)
  n_test <- nrow(comps$n_mat)
  n_readers <- rowSums(comps$n_mat)
  n_reader <- sum(n_readers)
  MS <- comps$MS
  cov <- comps$cov

  test_levels <- levels(object)$test

  n_formula <- (n_reader - sum(n_readers^2)) / (n_reader * (n_test - 1))
  denominator <- comps$var - cov[2] + n_formula * max(cov[2] - cov[3], 0)
  test_equality <- data.frame(
    `MS(T)` = MS[["T"]],
    Cov2 = cov[2],
    Cov3 = cov[3],
    Denominator = denominator,
    X2 = (n_test - 1) * MS[["T"]] / denominator,
    df = n_test - 1,
    check.names = FALSE
  )
  test_equality$`p-value` <- with(test_equality, 1 - pchisq(X2, df))

  test_means_list <- lapply(object$mrmc_tests, summary.mrmc_tests_frrc,
                            conf.level = conf.level)
  test_means <- do.call(rbind, test_means_list)
  rownames(test_means) <- names(test_means_list)

  estimates <- test_means$Estimate
  combs <- combinations(length(estimates), 2)
  test_diffs <- data.frame(
    Comparison = paste(test_levels[combs[, 1]], "-", test_levels[combs[, 2]]),
    Estimate = estimates[combs[, 1]] - estimates[combs[, 2]],
    StdErr = sqrt((1 / n_readers[combs[, 1]] + 1 / n_readers[combs[, 2]]) *
                    denominator),
    check.names = FALSE
  )
  test_diffs$CI <- with(test_diffs, {
    Estimate + qnorm((1 + conf.level) / 2) * StdErr %o% c(Lower = -1, Upper = 1)
  })
  test_diffs$z <- with(test_diffs, Estimate / StdErr)
  test_diffs$`p-value` <- with(test_diffs, 2 * (1 - pnorm(abs(z))))

  reader_means <- object$data
  reader_means$StdErr <- sqrt(diag(object$cov))
  reader_means$CI <- reader_means[[3]] + qnorm((1 + conf.level) / 2) *
    reader_means$StdErr %o% c(Lower = -1, Upper = 1)

  res <- new_summary_mrmc(object,
                          conf.level = conf.level,
                          vcov_comps = summary(comps),
                          test_equality = test_equality,
                          test_diffs = test_diffs,
                          reader_means = reader_means)
  structure(res, class = c("summary.mrmc_frrc", class(res)))
}


.summary.mrmc_rrfc <- function(object, conf.level, ...) {
  comps <- vcov_comps(object)
  n <- comps$n
  MS <- comps$MS
  cov <- comps$cov

  test_levels <- levels(object)$test

  test_equality <- data.frame(
    `MS(T)` = MS[["T"]],
    `MS(T:R)` = MS[["T:R"]],
    F = MS[["T"]] / MS[["T:R"]],
    df1 = n[["test"]] - 1,
    df2 = (n[["test"]] - 1) * (n[["reader"]] - 1),
    check.names = FALSE
  )
  test_equality$`p-value` <- with(test_equality, 1 - pf(F, df1, df2))

  test_means_list <- lapply(object$mrmc_tests, summary.mrmc_tests_rrfc,
                            conf.level = conf.level)
  test_means <- do.call(rbind, test_means_list)
  rownames(test_means) <- names(test_means_list)

  estimates <- test_means$Estimate
  combs <- combinations(length(estimates), 2)
  test_diffs <- data.frame(
    Comparison = paste(test_levels[combs[, 1]], "-", test_levels[combs[, 2]]),
    Estimate = estimates[combs[, 1]] - estimates[combs[, 2]],
    df = (n[["test"]] - 1) * (n[["reader"]] - 1),
    StdErr = sqrt(2 / n[["reader"]] * MS[["T:R"]])
  )
  test_diffs$CI <- with(test_diffs, {
    Estimate + qt((1 + conf.level) / 2, df) * StdErr %o% c(Lower = -1, Upper = 1)
  })
  test_diffs$t <- with(test_diffs, Estimate / StdErr)
  test_diffs$`p-value` <- with(test_diffs, 2 * (1 - pt(abs(t), df)))

  res <- new_summary_mrmc(object,
                          conf.level = conf.level,
                          vcov_comps = summary(comps),
                          test_equality = test_equality,
                          test_diffs = test_diffs,
                          test_means = test_means)
  structure(res, class = c("summary.mrmc_rrfc", class(res)))
}


.summary_nested.mrmc_rrfc <- function(object, conf.level, ...) {
  comps <- vcov_comps(object)
  n_test <- nrow(comps$n_mat)
  n_readers <- rowSums(comps$n_mat)
  n_reader <- sum(n_readers)
  MS <- comps$MS
  cov <- comps$cov

  test_levels <- levels(object)$test

  test_equality <- data.frame(
    `MS[T]` = MS[["T"]],
    `MS[R(T)]` = MS[["R"]],
    F = MS[["T"]] / MS[["R"]],
    df1 = n_test - 1,
    df2 = n_reader - n_test,
    check.names = FALSE
  )
  test_equality$`p-value` <- with(test_equality, 1 - pf(F, df1, df2))

  test_means_list <- lapply(object$mrmc_tests, summary.mrmc_tests_rrfc,
                            conf.level = conf.level)
  test_means <- do.call(rbind, test_means_list)
  rownames(test_means) <- names(test_means_list)

  estimates <- test_means$Estimate
  combs <- combinations(length(estimates), 2)
  test_diffs <- data.frame(
    Comparison = paste(test_levels[combs[, 1]], "-", test_levels[combs[, 2]]),
    Estimate = estimates[combs[, 1]] - estimates[combs[, 2]],
    df = n_reader - n_test,
    StdErr = sqrt((1 / n_readers[combs[, 1]] + 1 / n_readers[combs[, 2]]) *
                    MS[["R"]]),
    check.names = FALSE
  )
  test_diffs$CI <- with(test_diffs, {
    Estimate + qt((1 + conf.level) / 2, df) * StdErr %o% c(Lower = -1, Upper = 1)
  })
  test_diffs$t <- with(test_diffs, Estimate / StdErr)
  test_diffs$`p-value` <- with(test_diffs, 2 * (1 - pt(abs(t), df)))

  res <- new_summary_mrmc(object,
                          conf.level = conf.level,
                          vcov_comps = summary(comps),
                          test_equality = test_equality,
                          test_diffs = test_diffs,
                          test_means = test_means)
  structure(res, class = c("summary.mrmc_frrc", class(res)))
}


.summary.mrmc_lme <- function(object, conf.level, ...) {

  coef <- object$lme_fit$coef
  cov <- object$lme_fit$cov$R
  cov0 <- object$lme_fit$cov$R0

  n <- nrow(object$data)
  ntests <- nlevels(object$data[[object$vars["test"]]])

  x <- diag(ntests)
  x[, 1] <- 1
  x <- cbind(x, matrix(0, ntests, length(coef) - ntests))

  combs <- combinations(ntests, 2)
  x_diff <- x[combs[, 1], , drop = FALSE] - x[combs[, 2], , drop = FALSE]

  est_diff <- as.numeric(x_diff %*% coef)
  var_diff <- diag(x_diff %*% cov %*% t(x_diff))
  var0_diff <- diag(x_diff %*% cov0 %*% t(x_diff))

  aov <- object$aov
  aov_assign <- aov$assign[aov$qr$pivot[1:aov$rank]]
  df0 <- nrow(aov$model) - sum(aov_assign %in% c(0, 1, 3))
  df_satt <- var_diff^2 / (var0_diff^2 / df0)

  test_diffs <- data.frame(
    Comparison = paste(combs[, 1], "-", combs[, 2]),
    Estimate = est_diff,
    StdErr = sqrt(var_diff),
    df = df_satt
  )
  test_diffs$CI <- with(test_diffs, {
    Estimate + qt((1 + conf.level) / 2, df) * StdErr %o% c(Lower = -1, Upper = 1)
  })
  test_diffs$t <- with(test_diffs, Estimate / StdErr)
  test_diffs$`p-value` <- with(test_diffs, 2 * (1 - pt(abs(t), df)))

  res <- new_summary_mrmc(object,
                          conf.level = conf.level,
                          vcov_comps = summary(vcov_comps(object)),
                          test_diffs = test_diffs)
  structure(res, class = c("summary.mrmc_lme", class(res)))
}
