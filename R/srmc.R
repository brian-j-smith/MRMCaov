#' Single-Reader Multi-Case ROC Analysis
#'
#' Estimation and comparison of ROC performance metrics for single-reader
#' multi-case studies.
#'
#' @param response response metric expressed in terms of a package-supplied
#'   performance \code{\link[=metrics]{metric}}.
#' @param test variable of test identifiers.
#' @param case variable of case identifiers.
#' @param data data frame containing the \code{response}, \code{test}, and
#'   \code{case} variables.
#' @param method function, function call, or character string naming the
#'   \code{\link[=cov_methods]{method}} to use in calculating performance
#'   metric covariances.
#'
#' @seealso \code{\link{metrics}}, \code{\link{cov_methods}},
#' \code{\link{parameters}}, \code{\link{plot}}, \code{\link{roc_curves}},
#' \code{\link{summary}}
#'
srmc <- function(response, test, case, data, method = jackknife) {

  args <- c(substitute(response), substitute(test), substitute(case))
  data <- data[unique(unlist(mapply(all.vars, args)))]
  data <- cbind(data, reader = 1)
  names(data) <- make.unique(names(data), sep = "_")
  reader <- as.name(tail(names(data), 1))

  object <- eval(substitute(
    mrmc(response = response, test = test, reader = fixed(reader), case = case,
         data = data, method = method)
  ))
  object$call <- match.call()
  object

}
