utils::globalVariables(c("group", "x", "y"))


extract_vars <- function(formula) {
  vars <- all.vars(formula)
  c(observed = vars[1], predicted = vars[2], tests = vars[3], readers = vars[4],
    metric = all.names(formula)[2])
}


get_method <- function(x) {
  if (is(x, "character")) x <- get(x, mode = "function")
  if (is(x, "function")) x <- x()
  if (!is(x, "cov_method")) stop("invalid covariance method")
  x
}
