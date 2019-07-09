print.mrmc <- function(x, n = 20, ...) {
  model <- x$aov$model
  vars <- x$vars

  cat("Call:\n")
  print(x$call)
  
  cat("\nPositive truth status:", x$levels[2], "\n")
  
  cat("\nResponse metric estimates:\n\n")
  print(x$aov_data)
  cat("\nANOVA Table:\n\n")
  print(summary(x$aov))
  
  cat("\n\nObuchowski-Rockette error variance and covariance estimates:\n\n")
  comps <- vcov_comps(x, design = 1)
  vcov_comps <- data.frame(
    Estimate = c(comps$var, comps$cov),
    row.names = c("Error", "Cov1", "Cov2", "Cov3")
  )
  vcov_comps$Correlation = vcov_comps$Estimate / vcov_comps$Estimate[1]
  vcov_comps$Correlation[1] <- NA
  print(vcov_comps)
  
  cat("\n\nFirst", n, "ROC performance metrics:\n\n")
  print(head(x$roc, n = n))
  
  invisible()
}


print.summary.mrmc <- function(x, ...) {
  .print(x, ...)
}


.print <- function(x, ...) {
  UseMethod(".print")
}


.print.summary.mrmc_rrrc <- function(x, ...) {
  cat("Factor types: Random Readers and Random Cases\n\n")
  
  .print.summary.mrmc(x)
}


.print.summary.mrmc_frrc <- function(x, ...) {
  cat("Factor types: Fixed Readers and Random Cases\n\n")
  
  .print.summary.mrmc(x)
  
  cat("\n\nReader-specific ", 100 * x$conf.level, "% CIs and tests for ",
      x$vars["metric"], " pairwise differences (each analysis based only on",
      " data for the specified reader):\n\n", sep = "")
  print(x$reader_test_diffs)
  
  invisible()
}


.print.summary.mrmc_rrfc <- function(x, ...) {
  cat("Factor types: Random Readers and Fixed Cases\n\n")
  
  .print.summary.mrmc(x)
}


.print.summary.mrmc <- function(x, ...) {
  cat("Experimental design:",
      switch(x$design,
             "factorial",
             paste("cases nested within", x$vars["reader"]),
             paste("cases nested within", x$vars["test"])),
      "\n\n")
  
  cat("Obuchowski-Rockette variance component and covariance estimates:\n\n")
  print(x$vcov_comps)
  
  test_metric <- paste(x$vars["test"], x$vars["metric"])
  
  cat("\n\nANOVA global test of equal ", test_metric, ":\n\n", sep = "")
  print(x$test_equality)
  
  cat("\n\n", 100 * x$conf.level, "% CIs and tests for ", test_metric,
      " pairwise differences:\n\n", sep = "")
  print(x$test_diffs)
  
  cat("\n\n", 100 * x$conf.level, "% ", test_metric, " CIs (each analysis",
      " based only on data for the specified treatment):\n\n", sep = "")
  print(x$test_means)
  
  invisible()
}
