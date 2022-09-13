#' Multi-Reader Multi-Case ROC Analysis
#'
#' Estimation and comparison of ROC performance metrics for multi-reader
#' multi-case studies.
#'
#' @param response response metric expressed in terms of a package-supplied
#'   performance \code{\link[=metrics]{metric}}.
#' @param test variable of test identifiers.
#' @param reader variable of reader identifiers.
#' @param case variable of case identifiers.
#' @param data data frame containing the \code{response}, \code{test},
#'   \code{reader}, and \code{case} variables.
#' @param cov function, function call, or character string naming the
#'   \code{\link[=cov_methods]{method}} to use in calculating performance
#'   metric covariances.
#' @param method deprecated argument that will be removed in a future package
#'   version; use \code{cov} instead.
#' @param design one of the following study designs: 1 = factorial, 2 = cases
#'   nested within readers, 3 = cases nested within tests, or \code{NULL} to
#'   automatically set the design based on variable codings in data.
#'
#' @details
#' Readers and cases are treated as random factors by default.  Either one may
#' be designated as fixed in calls to \code{mrmc} with the syntax
#' \code{fixed(<variable name>)}, where \code{<variable name>} is the name of
#' the reader or case variable.
#'
#' @return
#' Returns an \code{mrmc} class object with the following elements.
#' \describe{
#'   \item{\code{design}}{experimental study design: 1 = factorial, 2 = cases
#'     nested within readers, 3 = cases nested within tests.}
#'   \item{\code{vars}}{character names of the analysis factors and reader
#'     performance metric.}
#'   \item{\code{fixed}}{logicals indicating whether the reader and case factors
#'     are treated as fixed in the analysis.}
#'   \item{\code{aov}}{results from an ordinary analysis of variance.}
#'   \item{\code{data}}{data frame of computed reader performance metrics for
#'     the analysis of variance.}
#'   \item{\code{num_obs}}{number of case observations for each of the computed
#'     metrics.}
#'   \item{\code{cov}}{reader performance covariance matrix.}
#'   \item{\code{mrmc_data}}{data frame of case-specific reader ratings.}
#'   \item{\code{levels}}{character levels of the true case statuses.}
#' }
#'
#' @seealso \code{\link{metrics}}, \code{\link{cov_methods}},
#' \code{\link{parameters}}, \code{\link{plot}}, \code{\link{roc_curves}},
#' \code{\link{summary}}
#'
#' @references
#' Dorfman DD, Berbaum KS, and Metz CE (1992). Receiver operating characteristic
#' rating analysis. Generalization to the population of readers and patients
#' with the jackknife method. Investigative Radiology, 27: 723–731.
#'
#' Obuchowski NA and Rockette HE (1995). Hypothesis testing of diagnostic
#' accuracy for multiple readers and multiple tests: an ANOVA approach with
#' dependent observations. Communications in Statistics–Simulation and
#' Computation 24: 285–308.
#'
#' Hillis SL, Obuchowski NA, Schartz KM, and Berbaum KS (2005). A comparison of
#' the Dorfman-Berbaum-Metz and Obuchowski-Rockette methods for receiver
#' operating characteristic (ROC) data. Statisticsin Medicine, 24: 1579–1607.
#'
#' Hillis SL (2007). A comparison of denominator degrees of freedom methods for
#' multiple observer ROC analysis. Statistics in Medicine, 26: 596–619.
#'
#' Hillis SL, Berbaum KS, and Metz CE (2008). Recent developments in the
#' Dorfman-Berbaum-Metz procedure for multireader ROC study analysis. Academic
#' Radiology, 15: 647–661.
#'
#' @examples
#' \donttest{
#' ## Random readers and cases
#' (est <- mrmc(empirical_auc(truth, rating), treatment, reader, case,
#'              data = VanDyke))
#' plot(est)
#' summary(est)
#'
#' ## Fixed readers and random cases
#' est <- mrmc(empirical_auc(truth, rating), treatment, fixed(reader), case,
#'             data = VanDyke)
#' summary(est)
#' }
#'
mrmc <- function(
  response, test, reader, case, data, cov = method, method = jackknife,
  design = NULL
) {

  dep_methodarg(missing(method))

  object <- eval(substitute(
    new_mrmc(response, test, reader, case, data, cov = cov, design = design)
  ))
  object$mrmc_tests <- mrmc_tests(object$data, object$cov, object$design)
  object$call <- match.call()

  mrmc_class <- if (all(object$fixed)) {
    stop("only one of reader or case may be fixed")
  } else if (is_one_reader(object) && !object$fixed["reader"]) {
    stop("reader must be fixed if there is only one")
  } else if (object$fixed["reader"]) {
    "mrmc_frrc"
  } else if (object$fixed["case"]) {
    "mrmc_rrfc"
  } else {
    "mrmc_rrrc"
  }

  structure(object, class = c(mrmc_class, class(object)))

}


mrmc_lme <- function(
  formula, test, reader, case, data, cov = method, method = jackknife,
  design = NULL
) {

  stopifnot(is(formula, "formula"))
  if (length(formula) != 3) stop("formula requires left and right terms")

  dep_methodarg(missing(method))

  response <- formula[[2]]
  object <- eval(substitute(
    new_mrmc(
      response, test, reader, case, data, cov = cov, design = design,
      types = "random"
    )
  ))

  args <- get_lme_args(formula, object, data)
  args$f <- negloglik_lme_reml
  args$grad <- grad_lme_reml

  args0 <- args
  args0$R <- (args$var - sum(args$cov)) * diag(length(args$y))

  params <- do.call(get_lme_params, args)
  params0 <- do.call(get_lme_params, args0)

  object$data <- args$data
  object$lme_fit <- list(
    coef = params$coef,
    cov = list(R = params$cov, R0 = params0$cov),
    optim = list(R = params$optim, R0 = params0$optim)
  )
  object$call <- match.call()

  structure(object, class = c("mrmc_lme", class(object)))

}


new_mrmc <- function(
  response, test, reader, case, data, cov, design, types = c("random", "fixed")
) {

  terms <- eval(substitute(
    mrmc_terms(response, test, reader, case, types = types)
  ))

  response_call <- match.call(get(terms$metric), terms$response)
  mrmc_data <- data.frame(
    truth = factor(eval(response_call$truth, data)),
    rating = as.numeric(eval(response_call$rating, data)),
    test = factor(data[[terms$labels["test"]]]),
    reader = factor(data[[terms$labels["reader"]]]),
    case = factor(data[[terms$labels["case"]]])
  )
  response_call[c(2, 3)] <- c(quote(truth), quote(rating))
  attr(mrmc_data, "metric_call") <- response_call

  mrmc_data <- preprocess(mrmc_data)

  design_data <- get_design(mrmc_data)
  if (is.null(design_data)) {
    stop("data factor codings are not a supported study design")
  } else if (is.null(design)) {
    design <- design_data
  } else if (design_data != design) {
    stop("data factor codings do not match study design ", design)
  }

  covmat <- if (!terms$fixed["case"]) get_cov_method(cov)(mrmc_data)

  var_names <- c("test", "reader")
  aov_data <- unique(mrmc_data[var_names])
  sort_order <- do.call(order, aov_data)
  aov_data <- aov_data[sort_order, , drop = FALSE]
  y <- num_obs <- numeric(nrow(aov_data))
  for (i in 1:nrow(aov_data)) {
    split <- merge(mrmc_data, aov_data[i, , drop = FALSE])
    y[i] <- eval(response_call, split)
    num_obs[i] <- nrow(split)
  }
  aov_data <- cbind(y, aov_data)
  names(aov_data) <- c(terms$metric, terms$labels[var_names])
  aovfit <- aov_mrmc(update(formula(aov_data), . ~ .^2), aov_data)

  structure(
    list(
      design = design,
      vars = c(terms$labels, metric = terms$metric),
      fixed = terms$fixed,
      aov = aovfit,
      data = aov_data,
      num_obs = num_obs,
      cov = covmat,
      mrmc_data = mrmc_data,
      levels = levels(mrmc_data$truth)
    ),
    class = "mrmc"
  )

}


mrmc_terms <- function(response, test, reader, case, types) {
  args <- eval(substitute(
    alist(response = response, test = test, reader = reader, case = case)
  ))

  test <- extract_term(args$test, types = "fixed")
  reader <- extract_term(args$reader, types = types)
  case <- extract_term(args$case, types = types)

  list(
    response = args$response,
    metric = as.character(args$response[[1]]),
    labels = c(test = test$label, reader = reader$label, case = case$label),
    fixed = c(reader = reader$type, case = case$type) == "fixed"
  )
}


extract_term <- function(x, types) {
  if (is.symbol(x)) {
    term_symbol <- x
    term_type <- types[1]
  } else if (is.call(x) && length(x) == 2) {
    term_symbol <- x[[2]]
    term_type <- as.character(x[[1]])
  } else {
    term_symbol <- NULL
  }
  if (!(is.symbol(term_symbol) && (term_type %in% types))) {
    stop("invalid mrmc term syntax: ", deparse(x), call. = FALSE)
  }
  list(label = as.character(term_symbol), type = term_type)
}
