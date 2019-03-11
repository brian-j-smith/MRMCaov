#' Multi-Reader Multi-Case ROC Analysis
#' 
#' Estimation and comparison of ROC performance metrics for multi-reader
#' multi-case studies.
#' 
#' @param formula formula defining the response metric and variables on the left
#' hand side and the test and reader variables on the right hand side.
#' The response should be expressed in terms of a package-supplied performance
#' \code{\link[=metrics]{metric}}.  Test and reader variables should be specified in that
#' order and separated by the \code{+} operator.
#' @param cases variable containing the case identifiers.
#' @param data data frame containing the \code{formula} and \code{cases}
#' variables.
#' @param method function, function call, or character string naming the
#' \code{\link[=cov_methods]{method}} to use in calculating performance
#' metric covariances.
#' @param fixed formula whose right hand side specifies one of the reader or
#' case variables as being a fixed factor in the analysis, or \code{NULL} if
#' both are random.
#' @param design one of the following study designs: 1 = factorial, 2 = cases
#' nested within readers, or 3 = cases nested within tests.
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
mrmc <- function(formula, cases, data, method = jackknife, fixed = NULL,
                 design = 1) {
  formula <- update(formula, ~ .^2)
  vars <- extract_vars(formula)
  vars["cases"] <- as.character(substitute(cases))

  factor_vars <- c("observed", "tests", "readers")
  data[vars[factor_vars]] <- lapply(data[vars[factor_vars]], factor)
  data[["(cases)"]] <- factor(data[[vars["cases"]]])

  cov <- get_method(method)(formula, data)

  fo <- formula
  df_by <- by(data, data[vars[3:4]], function(split) {
    structure(
      c(nrow(split), eval(fo[[2]], split)),
      names = c("N", vars["metric"])
    )
  })
  df <- cbind(expand.grid(dimnames(df_by)), do.call(rbind, df_by))
  fo[[2]] <- as.name(names(df)[ncol(df)])
  aovfit <- aov(fo, data = df)
  
  mrmc_class <- if (!is(fixed, "formula")) "mrmc_rrrc" else {
    fixed <- labels(terms(fixed))
    if (all(vars["readers"] == fixed)) "mrmc_frrc" else
      if (all(vars["cases"] == fixed)) "mrmc_rrfc" else
        stop("fixed formula must identify one of the readers or cases variable")
  }
  
  structure(
    list(call = sys.call(),
         design = design,
         vars = vars,
         roc = do.call(roc, c(list(data[[vars["observed"]]],
                                   data[[vars["predicted"]]]),
                              data[vars[c("tests", "readers")]])),
         aov = aovfit,
         aov_data = df,
         cov = cov,
         mrmc_tests = mrmc_tests(aovfit$model, cov, design)
         ),
    class = c(mrmc_class, "mrmc")
  )
}
