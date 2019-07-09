#' Multi-Reader Multi-Case ROC Analysis
#' 
#' Estimation and comparison of ROC performance metrics for multi-reader
#' multi-case studies.
#' 
#' @param response response metric expressed in terms of a package-supplied
#' performance \code{\link[=metrics]{metric}}.
#' @param test variable of test identifiers.
#' @param reader variable of reader identifiers.
#' @param case variable of case identifiers.
#' @param data data frame containing the \code{response}, \code{test},
#' \code{reader}, and \code{case} variables.
#' @param method function, function call, or character string naming the
#' \code{\link[=cov_methods]{method}} to use in calculating performance
#' metric covariances.
#' @param design one of the following study designs: 1 = factorial, 2 = cases
#' nested within readers, 3 = cases nested within tests, or \code{NULL} to
#' automatically set the design based on variable codings in data.
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
  
  terms <- eval(substitute(mrmc_terms(response, test, reader, case)))

  response_call <- match.call(get(terms$metric), terms$formula[[2]])
  mrmc_data <- structure(
    data.frame(truth = factor(eval(response_call$truth, data)),
               rating = eval(response_call$rating, data),
               test = factor(data[[terms$labels["test"]]]),
               reader = factor(data[[terms$labels["reader"]]]),
               case = factor(data[[terms$labels["case"]]])), 
    metric = terms$metric
  )

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
  df_by <- by(data, data[terms$labels[c("test", "reader")]], function(split) {
    structure(
      c(nrow(split), eval(fo[[2]], split)),
      names = c("N", terms$metric)
    )
  })
  df <- cbind(expand.grid(dimnames(df_by)), do.call(rbind, df_by))
  fo[[2]] <- as.name(names(df)[ncol(df)])
  aovfit <- aov(fo, data = df)
  
  mrmc_class <- if (all(terms$fixed)) {
    stop("only one of reader or case may be fixed")
  } else if (terms$fixed["reader"]) {
    "mrmc_frrc"
  } else if (terms$fixed["case"]) {
    "mrmc_rrfc"
  } else {
    "mrmc_rrrc"
  }

  structure(
    list(call = sys.call(),
         design = design,
         vars = c(terms$labels, metric = terms$metric),
         roc = do.call(roc, list(mrmc_data$truth, mrmc_data$rating,
                                 mrmc_data$test, mrmc_data$reader)),
         aov = aovfit,
         aov_data = df,
         cov = cov,
         mrmc_tests = mrmc_tests(aovfit$model, cov, design),
         levels = levels(mrmc_data$truth)),
    class = c(mrmc_class, "mrmc")
  )
}


mrmc_terms <- function(response, test, reader, case) {
  args <- eval(substitute(
    alist(response = response, test = test, reader = reader, case = case)
  ))
  
  test <- extract_term(args$test, type = c("", "fixed"))
  reader <- extract_term(args$reader)
  case <- extract_term(args$case)
  
  fo <- reformulate(c(test$label, reader$label), args$response)
  
  list(
    formula = update(fo, . ~ .^2),
    metric = all.names(fo)[2],
    labels = c(test = test$label, reader = reader$label, case = case$label),
    fixed = c(reader = reader$type, case = case$type) == "fixed"
  )
}


extract_term <- function(x, type = c("random", "fixed")) {
  if (is.symbol(x)) {
    term_symbol <- x
    term_type <- type[1]
  } else if (is.call(x) && length(x) == 2) {
    term_symbol <- x[[2]]
    term_type <- as.character(x[[1]])
  } else {
    term_symbol <- NULL
  }
  if (!(is.symbol(term_symbol) && (term_type %in% type))) {
    stop("invalid mrmc term syntax: ", deparse(x), call. = FALSE)
  } 
  list(label = as.character(term_symbol), type = term_type)
}
