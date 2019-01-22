utils::globalVariables(c("group", "x", "y"))


dim.mrmc <- function(x) {
  sapply(levels(x), length)
}


dim.mrmc_tests <- function(x) {
  sapply(levels(x), length)
}


extract_vars <- function(formula) {
  vars <- all.vars(formula)
  c(observed = vars[1], predicted = vars[2], tests = vars[3], readers = vars[4],
    metric = all.names(formula)[2])
}


get_method <- function(x) {
  if (is(x, "character")) x <- get(x, mode = "function")
  if (!is(x, "cov_method")) x <- x()
  if (!is(x, "cov_method")) stop("invalid covariance method")
  x
}


levels.mrmc <- function(x) {
  structure(x$aov$xlevels, names = c("tests", "readers"))
}


levels.mrmc_tests <- function(x) {
  structure(x$aov$xlevels, names = c("readers"))
}


mean.mrmc_tests <- function(x, ...) {
  mean(x$aov$model[[1]])
}


meansq <- function(x, ...) {
  UseMethod("meansq")
}


meansq.mrmc <- function(x, ...) {
  structure(
    summary(x$aov)[[1]][["Mean Sq"]],
    names = c("T", "R", "T:R")
  )
}


meansq.mrmc_tests <- function(x, ...) {
  c("T" = 0, "R" = summary(x$aov)[[1]][["Mean Sq"]], "T:R" = 0)
}


vcov_comps <- function(object, ...) {
  UseMethod("vcov_comps")
}


vcov_comps.mrmc <- function(object, design = object$design, test = NULL,
                            reader = NULL, ...) {
  model <- object$aov$model

  tests <- model[[2]]
  readers <- model[[3]]
  
  same_test <- outer(tests, tests, "==")
  same_reader <- outer(readers, readers, "==")
  
  is_group <- rep(TRUE, nrow(model))
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
  }
  
  structure(
    list(vars = object$vars,
         n = dim(object),
         MS = meansq(object),
         var = mean(diag(object$cov[is_group, is_group])),
         cov = cov),
    class ="vcov_comps"
  )
}


vcov_comps.mrmc_tests <- function(object, design = object$design, ...) {
  model <- object$aov$model

  cov2 <- if (design == 2) 0 else {
    readers <- model[[2]]
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
  
  vcov_comps <- data.frame(
    Estimate = c(
      (MS[["R"]] - MS[["T:R"]]) / n[["tests"]] - cov[1] + cov[3],
      MS[["T:R"]] - var_error + cov[1] + (cov[2] - cov[3]),
      var_error,
      cov
    ),
    row.names = c(object$vars["readers"],
                  paste0(object$vars[c("tests", "readers")], collapse = ":"),
                  "Error", "Cov1", "Cov2", "Cov3")
  )
  vcov_comps$Correlation <- vcov_comps$Estimate /
    vcov_comps["Error", "Estimate"]
  vcov_comps$Correlation[1:3] <- NA
  
  vcov_comps
}
