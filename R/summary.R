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
  .summary(object, conf.level = conf.level, ...)
}


.summary <- function(object, ...) {
  UseMethod(".summary")
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

  structure(
    list(design = object$design,
         vars = object$vars,
         conf.level = conf.level,
         vcov_comps = summary(comps),
         test_equality = test_equality,
         test_diffs = test_diffs,
         test_means = test_means),
    class = c("summary.mrmc_rrrc", "summary.mrmc")
  )
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

  structure(
    list(design = object$design,
         vars = object$vars,
         conf.level = conf.level,
         vcov_comps = summary(comps),
         test_equality = test_equality,
         test_diffs = test_diffs,
         test_means = test_means,
         reader_test_diffs = reader_test_diffs(object, conf.level)),
    class = c("summary.mrmc_frrc", "summary.mrmc")
  )
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

  structure(
    list(design = object$design,
         vars = object$vars,
         conf.level = conf.level,
         vcov_comps = summary(comps),
         test_equality = test_equality,
         test_diffs = test_diffs,
         test_means = test_means),
    class = c("summary.mrmc_rrfc", "summary.mrmc")
  )
}
