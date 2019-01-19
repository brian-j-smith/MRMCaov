mrmc_tests <- function(data, cov, design) {
  tests <- data[[2]]
  data <- data[-2]
  formula <- formula(data)
  
  lapply(levels(tests), function(test) {
    is_test <- tests == test
    structure(
      list(design = design,
           aov = aov(formula, data = subset(data, is_test)),
           cov = cov[is_test, is_test]),
      class = "mrmc_tests"
    )
  })
}


summary.mrmc_tests_rrrc <- function(object, conf.level = 0.95, ...) {
  comps <- vcov_comps(object)
  n <- comps$n
  MS <- comps$MS
  cov <- comps$cov

  df <- data.frame(
    Estimate = mean(object),
    `MS(R)` = MS[["R"]],
    Cov2 = cov[2],
    StdErr = sqrt(MS[["R"]] / n[["readers"]] + max(cov[2], 0)),
    df = (MS[["R"]] + n[["readers"]] * max(cov[2], 0))^2 /
      (MS[["R"]]^2 / (n[["readers"]] - 1)),
    check.names = FALSE
  )
  df$CI <- with(df, {
    Estimate + qt((1 + conf.level) / 2, df) * StdErr %o% c(Lower = -1, Upper = 1)
  }) %>% pmin(1) %>% pmax(0)

  df
}


summary.mrmc_tests_frrc <- function(object, conf.level = 0.95, ...) {
  comps <- vcov_comps(object)
  n <- comps$n
  MS <- comps$MS
  cov <- comps$cov

  df <- data.frame(
    Estimate = mean(object),
    `Var(Error)` = comps$var,
    Cov2 = cov[2],
    StdErr = sqrt((comps$var + (n[["readers"]] - 1) * max(cov[2], 0)) /
                    n[["readers"]]),
    check.names = FALSE
  )
  df$CI <- with(df, {
    Estimate + qnorm((1 + conf.level) / 2) * StdErr %o% c(Lower = -1, Upper = 1)
  }) %>% pmin(1) %>% pmax(0)

  df
}


summary.mrmc_tests_rrfc <- function(object, conf.level = 0.95, ...) {
  comps <- vcov_comps(object)
  n <- comps$n
  MS <- comps$MS
  cov <- comps$cov

  df <- data.frame(
    Estimate = mean(object),
    `MS(R)` = MS[["R"]],
    StdErr = sqrt(MS[["R"]] / n[["readers"]]),
    df = n[["readers"]] - 1,
    check.names = FALSE
  )
  df$CI <- with(df, {
    Estimate + qt((1 + conf.level) / 2, df) * StdErr %o% c(Lower = -1, Upper = 1)
  }) %>% pmin(1) %>% pmax(0)

  df
}
