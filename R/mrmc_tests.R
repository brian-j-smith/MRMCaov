mrmc_tests <- function(formula, data, method) {
  vars <- extract_vars(formula)
  by(data, data[vars["tests"]], function(split_test) {
    df <- by(split_test, split_test[vars["readers"]], function(split_reader) {
      eval(formula[[2]], split_reader)
    }) %>% as.data.frame.table(responseName = vars["metric"])
    fo <- reformulate(vars["readers"], names(df)[ncol(df)])
    structure(
      list(aov = aov(fo, data = df),
           vars = vars,
           cov = method(formula, droplevels(split_test))),
      class = "mrmc_tests"
    )
  })
}


summary.mrmc_tests <- function(object, design = 1, conf.level = 0.95, ...) {
  MS <- summary(object$aov)[[1]][["Mean Sq"]]
  
  model <- object$aov$model
  vars <- object$vars
  n_readers <- nlevels(model[[vars["readers"]]])
  
  cov2 <- ifelse(design == 2, 0,
                 get_cov_comps(object$cov, rep(1, n_readers),
                               model[[vars["readers"]]])[2])

  df <- data.frame(
    "Estimate" = mean(model[[1]]),
    "MS(R)" = MS,
    "Cov2" = cov2,
    "StdErr" = sqrt(MS / n_readers + max(cov2, 0)),
    "df" = (MS + n_readers * max(cov2, 0))^2 / (MS^2 / (n_readers - 1)),
    check.names = FALSE
  )
  df$CI <- (df$Estimate + qt((1 + conf.level) / 2, df$df) *
    df$StdErr %o% c(Lower = -1, Upper = 1)) %>% pmin(1) %>% pmax(0)

  df
}
