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
#' @param method function, function call, or character string naming the
#'   \code{\link[=cov_methods]{method}} to use in calculating performance
#'   metric covariances.
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
#' @seealso \code{\link{metrics}}, \code{\link{cov_methods}},
#' \code{\link{summary}}, \code{\link{plot}}
#'
#' @examples
#' ## Random readers and cases
#' (est <- mrmc(empirical_auc(truth, rating), treatment, reader, case,
#'              data = VanDyke))
#' summary(est)
#' plot(est)
#'
#' ## Fixed readers and fixed cases
#' est <- mrmc(empirical_auc(truth, rating), treatment, fixed(reader), case,
#'             data = VanDyke)
#' summary(est)
#'
mrmc <- function(response, test, reader, case, data, method = jackknife,
                 design = NULL) {

  object <- eval(substitute(
    new_mrmc(response, test, reader, case, data, method = method,
             design = design)
  ))
  object$mrmc_tests <- mrmc_tests(object$aov$model, object$cov, object$design)
  object$call <- sys.call()

  mrmc_class <- if (all(object$fixed)) {
    stop("only one of reader or case may be fixed")
  } else if (object$fixed["reader"]) {
    "mrmc_frrc"
  } else if (object$fixed["case"]) {
    "mrmc_rrfc"
  } else {
    "mrmc_rrrc"
  }

  structure(object, class = c(mrmc_class, class(object)))

}


mrmc_lme <- function(formula, test, reader, case, data, method = jackknife,
                     design = NULL) {

  stopifnot(is(formula, "formula"))
  if (length(formula) != 3) stop("formula requires left and right terms")

  response <- formula[[2]]
  object <- eval(substitute(
    new_mrmc(response, test, reader, case, data, method = method,
             design = design, types = "random")
  ))

  args <- get_lme_args(formula, object, data)
  args$f <- negloglik_lme_reml
  args$grad <- grad_lme_reml

  args0 <- args
  args0$R <- (args$var - sum(args$cov)) * diag(length(args$y))

  params <- do.call(get_lme_params, args)
  params0 <- do.call(get_lme_params, args0)

  object$lme_fit <- list(
    coef = params$coef,
    cov = list(R = params$cov, R0 = params0$cov),
    optim = list(R = params$optim, R0 = params0$optim)
  )
  object$call <- sys.call()

  structure(object, class = c("mrmc_lme", class(object)))

}


new_mrmc <- function(response, test, reader, case, data, method, design,
                     types = c("random", "fixed")) {

  terms <- eval(substitute(
    mrmc_terms(response, test, reader, case, types = types)
  ))

  response_call <- match.call(get(terms$metric), terms$formula[[2]])
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

  cov <- get_method(method)(mrmc_data)

  fo <- terms$formula
  vars <- terms$labels[c("test", "reader")]
  mrmc_groups <- structure(mrmc_data[names(vars)], names = vars)
  df_by <- by(mrmc_data, mrmc_groups, function(split) {
    structure(
      c(nrow(split), eval(response_call, split)),
      names = c("N", terms$metric)
    )
  })
  df <- cbind(expand.grid(dimnames(df_by)), do.call(rbind, df_by))
  fo[[2]] <- as.name(names(df)[ncol(df)])
  aovfit <- aov(fo, data = df)

  structure(
    list(design = design,
         vars = c(terms$labels, metric = terms$metric),
         fixed = terms$fixed,
         aov = aovfit,
         aov_data = df,
         cov = cov,
         mrmc_data = mrmc_data,
         levels = levels(mrmc_data$truth)),
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

  fo <- reformulate(c(test$label, reader$label), args$response)

  list(
    formula = update(fo, . ~ .^2),
    metric = all.names(fo)[2],
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
