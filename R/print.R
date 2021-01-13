#' Print ROC Objects
#'
#' Print ROC objects from the \pkg{MRMCaov} package.
#'
#' @rdname print
#'
#' @param x object to print.
#' @param n number of ROC curve points to print.
#' @param n_curves number of ROC curves to print.
#' @param ... arguments passed to other methods.
#'
#' @examples
#' curves <- with(VanDyke,
#'   roc_curves(truth, rating, groups = list(Test = treatment, Reader = reader))
#' )
#' print(curves)
#'
print.roc_curve <- function(x, n = 11, ...) {
  cat("Points:\n")
  print(points(x, values = seq(0, 1, length = n), ...))
  invisible(x)
}


print.binormal_curve <- function(x, ...) {
  params <- parameters(x)
  cat("Binormal Curve\n",
      "Parameters:",
      paste(names(params), format(c(params)), sep = " = ", collapse = ", "),
      "\n", sep = "")
  NextMethod()
}


print.binormalLR_curve <- function(x, ...) {
  params <- parameters(x)
  format_params <- function(x) {
    paste(names(x), format(c(x)), sep = " = ", collapse = ", ")
  }
  cat("Binormal Likelihood Ratio Curve\n",
      "Parameters\n",
      "  Metz and Pan: ", format_params(params$Metz), "\n",
      "  Bi-Chi-Square: ", format_params(params$bichisquare), "\n",
      "  Binormal: ", format_params(params$binormal), "\n", sep = "")
  NextMethod()
}


print.empirical_curve <- function(x, ...) {
  cat("Empirical Curve\n")
  NextMethod()
}


#' @rdname print
#'
print.roc_curves <- function(x, n_curves = 5, n = 11, ...) {
  cat("ROC Curves\n\n")
  n_more <- nrow(x) - n_curves
  n_curves <- min(n_curves, nrow(x))
  vsep <- strrep("-", 0.75 * getOption("width"))
  for (i in seq_len(n_curves)) {
    if (i != 1) cat(vsep, "\n")
    cat(paste0(names(x$Group), ": ",
               as.character(x$Group[i, ]),
               collapse = "\n"), "\n")
    print(x$Curve[[i]], n = n, ...)
  }
  if (n_more) cat("... with", n_more, "more curves\n")
  invisible(x)
}


print.cov_matrix <- function(x, ...) {
  print(as(x, "matrix"))
}


print.mrmc <- function(x, ...) {
  cat("Call:\n")
  print(x$call)

  cat("\nPositive truth status:", x$levels[2], "\n")

  cat("\nResponse metric data:\n\n")
  print(tibble(N = x$num_obs, data = x$data))
  cat("\nANOVA Table:\n\n")
  print(summary(x$aov))

  cat("\n\nObuchowski-Rockette error variance and covariance estimates:\n\n")
  if (is.null(x$cov)) {
    cat("Not applicable because cases are fixed\n")
  } else {
    comps <- vcov_comps(x, design = 1)
    vcov_comps <- data.frame(
      Estimate = c(comps$var, comps$cov),
      row.names = c("Error", "Cov1", "Cov2", "Cov3")
    )
    vcov_comps$Correlation = vcov_comps$Estimate / vcov_comps$Estimate[1]
    vcov_comps$Correlation[1] <- NA
    print(vcov_comps)
  }

  invisible(x)
}


print.summary.mrmc <- function(x, ...) {
  .print(x, ...)
}


.print <- function(x, ...) {
  UseMethod(".print")
}


.print.summary.mrmc_frrc <- function(x, ...) {
  cat(if (is_one_reader(x)) "Single" else "Multi",
      "-Reader Multi-Case Analysis of Variance\n",
      "Data: ", x$data_name, "\n",
      "Factor types: Fixed Readers and Random Cases\n",
      "Covariance method: ", x$cov_method, "\n\n",
      sep = "")

  .print.summary.mrmc(x)

  if (!is.null(x$reader_test_diffs)) {
    cat("\n\nReader-specific ", 100 * x$conf.level, "% CIs and tests for ",
        x$vars["metric"], " pairwise differences (each analysis based only on",
        " data for the specified reader):\n\n", sep = "")
    print(x$reader_test_diffs)
  }

  if (!is.null(x$reader_means)) {
    cat("\n\nSingle reader ", 100 * x$conf.level, "% CIs:\n\n", sep = "")
    print(x$reader_means)
  }

  invisible(x)
}


.print.summary.mrmc_rrfc <- function(x, ...) {
  cat("Multi-Reader Multi-Case Analysis of Variance\n",
      "Data: ", x$data_name, "\n",
      "Factor types: Random Readers and Fixed Cases\n",
      sep = "")

  .print.summary.mrmc(x)
}


.print.summary.mrmc_rrrc <- function(x, ...) {
  cat("Multi-Reader Multi-Case Analysis of Variance\n",
      "Data: ", x$data_name, "\n",
      "Factor types: Random Readers and Random Cases\n",
      "Covariance method: ", x$cov_method, "\n\n",
      sep = "")

  .print.summary.mrmc(x)
}


.print.summary.mrmc_lme <- function(x, ...) {
  cat("Multi-Reader Multi-Case Linear Mixed Effects Analysis\n",
      "Data: ", x$data_name, "\n",
      "Factor types: Random Readers and Random Cases\n",
      "Covariance method: ", x$cov_method, "\n\n",
      sep = "")

  .print.summary.mrmc(x)
}


.print.summary.mrmc <- function(x, ...) {
  cat("Experimental design:",
      switch(x$design,
             "factorial",
             paste("cases nested within", x$vars["reader"]),
             paste("cases nested within", x$vars["test"]),
             paste(x$vars["reader"], "nested within", x$vars["test"])),
      "\n")

  cat("\nObuchowski-Rockette variance component and covariance estimates:\n\n")
  if (is.null(x$vcov_comps)) {
    cat("Not applicable because cases are fixed\n")
  } else {
    print(if (is_one_reader(x)) x$vcov_comps[-(1:2), ] else x$vcov_comps)
  }

  test_metric <- paste(x$vars["test"], x$vars["metric"])

  if (!is.null(x$test_equality)) {
    cat("\n\nANOVA global test of equal ", test_metric, ":\n\n", sep = "")
    print(x$test_equality)
  }

  if (!is.null(x$test_diffs)) {
    cat("\n\n", 100 * x$conf.level, "% CIs and tests for ", test_metric,
        " pairwise differences:\n\n", sep = "")
    print(x$test_diffs)
  }

  if (!is.null(x$test_means)) {
    cat("\n\n", 100 * x$conf.level, "% ", test_metric, " CIs (each analysis",
        " based only on data for the specified treatment):\n\n", sep = "")
    print(x$test_means)
  }

  invisible(x)
}
