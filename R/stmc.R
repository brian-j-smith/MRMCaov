#' Single-Test (Single-Reader) Multi-Case ROC Analysis
#'
#' Estimation of ROC performance metrics for a single test of multiple cases.
#'
#' @param response response metric expressed in terms of a package-supplied
#'   performance \code{\link[=metrics]{metric}}.
#' @param case optional variable of case identifiers.
#' @param data data frame containing the \code{response}, \code{test},
#'   \code{reader}, and \code{case} variables.
#' @param method function, function call, or character string naming the
#'   \code{\link[=cov_methods]{method}} to use in calculating performance
#'   metric covariances.
#'
#' @seealso \code{\link{metrics}}, \code{\link{cov_methods}},
#' \code{\link{parameters}}, \code{\link{plot}}, \code{\link{roc_curves}},
#' \code{\link{summary}}
#'
#' @examples
#' VanDyke1 <- subset(VanDyke, reader == 1)
#' est <- stmc(empirical_auc(truth, rating), data = VanDyke1)
#' summary(est)
#'
stmc <- function(response, case, data, method = jackknife) {

  response_call <- substitute(response)
  metric <- as.character(response_call[[1]])
  response_call <- match.call(get(metric), response_call)

  case <- if (missing(case)) 1:nrow(data) else eval(substitute(case), data)

  stmc_data <- data.frame(
    truth = factor(eval(response_call$truth, data)),
    rating = as.numeric(eval(response_call$rating, data)),
    test = factor(1),
    reader = factor(1),
    case = factor(case)
  )
  response_call[c(2, 3)] <- c(quote(truth), quote(rating))
  attr(stmc_data, "metric_call") <- response_call

  structure(
    list(
      metric = metric,
      est = eval(response_call, stmc_data),
      se = sqrt(get_method(method)(stmc_data)[1]),
      stmc_data = stmc_data
    ),
    class = "stmc"
  )
}
