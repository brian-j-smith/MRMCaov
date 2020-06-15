#' ROC Performance Curves
#'
#' Calculation of TPR and FPR pairs for values of a numeric rating of a
#' true binary response.
#'
#' @rdname roc_curves
#'
#' @param truth vector of true binary statuses.
#' @param rating vector of numeric ratings.
#' @param groups list or data frame of grouping variables of the same lengths as
#'   \code{truth} and \code{rating}.
#' @param method character string indicating the curve type as
#'   \code{"binormal"}, \code{"empirical"}, \code{"trapezoidal"}, or
#'   \code{"proproc"} or the averaging of binormal curves over \code{"points"}
#'   or \code{"parameters"}.
#' @param x object returned by \code{\link{mrmc}} or \code{roc_curves} for which
#'   to compute points on or to average over the curves.
#' @param values numeric vector of values at which to compute the points.  If
#'   \code{NULL} then the intersection of all empirical points is used.
#' @param metric reader performance metric to which the \code{values}
#'   correspond.
#' @param ties function determining empirical roc points returned in cases of
#'   ties.
#' @param ... arguments passed from the \code{mean()} method to \code{points()}.
#'
#' @seealso \code{\link{plot}}
#'
#' @references
#' Chen W and Samuelson FW (2014). The average receiver operating characteristic
#' curve in multireader multicase imaging studies. The British Journal of
#' Radiology, 87(1040): 20140016.
#'
#' @examples
#' curves <- with(VanDyke,
#'   roc_curves(truth, rating, groups = list(Test = treatment, Reader = reader))
#' )
#' points(curves)
#' mean(curves)
#'
roc_curves <- function(...) {
  UseMethod("roc_curves")
}


#' @rdname roc_curves
#'
roc_curves.default <- function(truth, rating, groups = list(),
                               method = "empirical", ...) {
  method <- match.arg(method,
                      c("binormal", "empirical", "trapezoidal", "proproc"))

  if (method %in% c("empirical", "trapezoidal")) {
    method <- "empirical"
    roc_curve <- empirical_curve
  } else if (method %in% c("binormal", "proproc")) {
    roc_curve <- param_curve
  }

  if (length(groups)) {
    data <- tibble(truth = truth, rating = rating)
    curves <- by(data, groups, function(split) {
      roc_curve(split$truth, split$rating, method = method)
    })
    groups <- expand.grid(dimnames(curves))
    keep <- !is.null(curves)
    curves <- tibble(Group = groups[keep, , drop = FALSE], Curve = curves[keep])
    structure(curves,
              class = c(paste0(method, "_curves"), "roc_curves", class(curves)))
  } else {
    roc_curve(truth, rating, method = method)
  }
}


#' @rdname roc_curves
#'
roc_curves.mrmc <- function(x, ...) {
  roc_method <- unlist(strsplit(x$vars["metric"], "_"))[1]
  if (!(roc_method %in% c("binormal", "empirical", "proproc"))) {
    roc_method <- "empirical"
  }
  vars <- x$vars[c("reader", "test")]
  roc_curves(x$mrmc_data$truth, x$mrmc_data$rating,
             groups = structure(x$mrmc_data[names(vars)], names = vars),
             method = roc_method)
}


empirical_curve <- function(truth, rating, ...) {
  params <- empirical_params(truth, rating)
  roc <- pROC::roc(params$truth, params$rating, auc = FALSE, quiet = TRUE)
  curve <- tibble(FPR = 1 - roc$specificities, TPR = roc$sensitivities)
  structure(curve[nrow(curve):1, ], params = params,
            class = c("empirical_curve", "roc_curve", class(curve)))
}


param_curve <- function(truth, rating, method, ...) {
  roc_params <- switch(method,
                       binormal = binormal_params,
                       proproc = proproc_params)

  data <- tibble(truth = truth, rating = rating)
  roc <- pROC::roc(data$truth, data$rating, auc = FALSE, quiet = TRUE)
  spec <- rev(roc$specificities)
  params <- roc_params(data$truth, data$rating)
  curve <- tibble(FPR = 1 - spec, TPR = sensitivity(params, spec))

  structure(curve, params = params,
            class = c(paste0(method, "_curve"), "roc_curve", class(curve)))
}


binormal_params <- function(truth, rating) {
  truth <- as.factor(truth)
  is_pos <- truth == levels(truth)[2]
  pred_pos <- as.double(rating[is_pos])
  pred_neg <- as.double(rating[!is_pos])
  roc <- .Fortran("cvbmroc",
                  length(pred_neg), length(pred_pos), pred_neg, pred_pos,
                  a = double(1), b = double(1),
                  auc = double(1), auc_var = double(1))
  structure(list(a = roc$a, b = roc$b),
            class = c("binormal_params", "roc_params"))
}


empirical_params <- function(truth, rating) {
  structure(list(truth = truth, rating = rating),
            class = c("empirical_params", "roc_params"))
}


proproc_params <- function(truth, rating) {
  truth <- as.factor(truth)
  is_pos <- truth == levels(truth)[2]
  pred_pos <- as.double(rating[is_pos])
  pred_neg <- as.double(rating[!is_pos])
  roc <- .Fortran("pbmroc",
                  length(pred_neg), length(pred_pos), pred_neg, pred_pos,
                  d_a = double(1), c = double(1),
                  auc = double(1), auc_var = double(1))
  structure(list(d_a = roc$d_a, c = roc$c),
            class = c("proproc_params", "roc_params"))
}


#' @rdname roc_curves
#'
as.data.frame.roc_curves <- function(x, ...) {
  as.data.frame(as_tibble(x))
}


#' @rdname roc_curves
#'
as_tibble.roc_curves <- function(x, ...) {
  curves2tibble(x$Curve, x$Group)
}


#' @rdname roc_curves
#'
points.roc_curve <- function(x, metric = c("specificity", "sensitivity"),
                             values = seq(0, 1, length = 100), ...) {
  points(attr(x, "params"), metric = metric, values = values, ...)
}


#' @rdname roc_curves
#'
points.roc_curves <- function(x, metric = c("specificity", "sensitivity"),
                              values = seq(0, 1, length = 100), ...) {
  metric <- match.arg(metric)
  new_pts_list <- lapply(x$Curve, function(curve) {
    points(curve, metric = metric, values = values, ...)
  })
  new_pts <- curves2tibble(new_pts_list, x$Group)
  structure(new_pts, metric = metric,
            class = c("roc_points", class(new_pts)))
}


points.roc_params <- function(x, metric = c("specificity", "sensitivity"),
                              values = seq(0, 1, length = 100), ...) {
  metric <- match.arg(metric)
  switch(metric,
         specificity = {
           input_name <- "FPR"
           input <- sort(1 - values)
           output_name <- "TPR"
           output_fun <- function(x, input) sensitivity(x, 1 - input)
         },
         sensitivity = {
           input_name <- "TPR"
           input <- sort(values)
           output_name <- "FPR"
           output_fun <- function(x, input) 1 - specificity(x, input)
         })

  output <- output_fun(x, input)
  new_pts <- tibble(input, output)
  names(new_pts) = c(input_name, output_name)

  structure(new_pts[c("FPR", "TPR")], metric = metric,
            class = c("roc_points", class(new_pts)))
}


#' @rdname roc_curves
#'
points.empirical_curve <- function(x, metric = c("specificity", "sensitivity"),
                                   values = NULL, ties = max, ...) {
  metric <- match.arg(metric)

  new_pts <- if (!is.null(values)) {

    switch(metric,
           "specificity" = {
             input_name <- "FPR"
             output_name <- "TPR"
             output_fun <- function(x, input, ...) {
               sensitivity(x, 1 - input, ...)
             }
             sort_metric <- function(x) sort(1 - x)
           },
           sensitivity = {
             input_name <- "TPR"
             output_name <- "FPR"
             output_fun <- function(x, input, ...) {
               1 - specificity(x, input, ...)
             }
             sort_metric <- function(x) sort(x)
           })

    input <- sort_metric(values)
    output <- output_fun(x, input, ties = ties)
    pts <- tibble(input, output)
    names(pts) <- c(input_name, output_name)
    pts[c("FPR", "TPR")]

  } else as_tibble(x)

  structure(new_pts, metric = metric, class = c("roc_points", class(new_pts)))
}


#' @rdname roc_curves
#'
points.empirical_curves <- function(x, metric = c("specificity", "sensitivity"),
                                    values = NULL, ties = max, ...) {
  metric <- match.arg(metric)

  new_pts_list <- if (is.null(values)) {
    switch(metric,
           specificity = {
             input_name <- "FPR"
             pts_fun <- function(x, input, ...) {
               points(x, values = 1 - input, ...)
             }
           },
           sensitivity = {
             input_name <- "TPR"
             pts_fun <- function(x, input, ...) points(x, values = input, ...)
           })
    input_list <- mapply(getElement, x$Curve, input_name)
    input_dups_list <- lapply(input_list, function(input) {
      input[duplicated(input)]
    })
    input <- sort(unique(unlist(input_list)))
    input_dups <- sort(unique(unlist(input_dups_list)))
    lapply(x$Curve, function(curve) {
      pts <- pts_fun(curve, input, metric = metric, ties = max, ...)
      pts_dups <- pts_fun(curve, input_dups, metric = metric, ties = min, ...)
      pts <- rbind(pts, pts_dups)
      pts[order(pts$FPR, pts$TPR), c("FPR", "TPR")]
    })
  } else {
    lapply(x$Curve, function(curve) {
      points(curve, metric = metric, values = values, ties = ties, ...)
    })
  }
  new_pts <- curves2tibble(new_pts_list, x$Group)

  structure(new_pts, metric = metric, class = c("roc_points", class(new_pts)))
}


#' @rdname roc_curves
#'
mean.roc_curve <- function(x, ...) {
  points(x, ...)
}


#' @rdname roc_curves
#'
mean.roc_curves <- function(x, ...) {
  pts <- points(x, ...)
  metric <- attr(pts, "metric")
  switch(metric,
         specificity = {
           input_name <- "FPR"
           output_name <- "TPR"
         },
         sensitivity = {
           input_name <- "TPR"
           output_name <- "FPR"
         })
  output_list <- split(pts[[output_name]], pts[["Group"]])
  output <- rowMeans(do.call(cbind, output_list))
  input <- head(pts[[input_name]], n = length(output))
  new_pts <- tibble(input, output)
  names(new_pts) = c(input_name, output_name)
  structure(new_pts[c("FPR", "TPR")], metric = metric, class = class(pts))
}


mean.roc_params <- function(x, ...) {
  points(x, ...)
}


#' @rdname roc_curves
#'
mean.binormal_curves <- function(x, method = c("points", "parameters"), ...) {
  if (match.arg(method) == "parameters") {
    params_list <- mapply(attr, x$Curve, "params", SIMPLIFY = FALSE)
    auc_mean <- mean(sapply(params_list, auc))
    b_mean <- mean(mapply(getElement, params_list, "b"))
    a <- sqrt(1 + b_mean^2) * qnorm(auc_mean)
    params <- structure(list(a = a, b = b_mean),
                        class = class(params_list[[1]]))
    mean(params, ...)
  } else NextMethod()
}


auc <- function(x, ...) {
  UseMethod("auc")
}


auc.roc_curve <- function(x, ...) {
  auc(attr(x, "params"), ...)
}


auc.binormal_params <- function(x, partial = FALSE, min = 0, max = 1, ...) {
  if (isFALSE(partial)) {
    pnorm(x$a / sqrt(1 + x$b^2))
  } else {
    partial <- partial_auc_params(partial, min, max)
    .Fortran("cvbmrocpartial",
             x$a, x$b,
             partial$min, partial$max, partial$flag,
             est = double(1), err = integer(1))$est
  }
}


auc.empirical_params <- function(x, partial = FALSE, min = 0, max = 1, ...) {
  args <- list(pROC::roc(x$truth, x$rating, quiet = TRUE))
  if (!isFALSE(partial)) {
    partial <- match.arg(partial, c("sensitivity", "specificity"))
    args$partial.auc <- c(min, max)
    args$partial.auc.focus <- partial
  }
  as.numeric(do.call(pROC::auc, args))
}


auc.proproc_params <- function(x, partial = FALSE, min = 0, max = 1, ...) {
  if (isFALSE(partial)) {
    rho <- -1 * (1 - x$c^2) / (1 + x$c^2)
    rho <- rbind(c(1, rho), c(rho, 1))
    pnorm(x$d_a / sqrt(2)) +
      2 * as.numeric(pmvnorm(upper = c(-x$d_a / sqrt(2), 0), corr = rho))
  } else {
    partial <- partial_auc_params(partial, min, max)
    .Fortran("pbmrocpartial",
             x$d_a, x$c,
             partial$min, partial$max, partial$flag,
             est = double(1), err = integer(1))$est
  }
}


partial_auc_params <- function(x, min, max) {
  x <- match.arg(x, c("sensitivity", "specificity"))
  min <- as.double(min)
  max <- as.double(max)
  if (x == "specificity") {
    flag <- 1L
    min <- 1 - max
    max <- 1 - min
  } else {
    flag <- 2L
  }
  list(min = min, max = max, flag = flag)
}


sensitivity <- function(x, ...) {
  UseMethod("sensitivity")
}


sensitivity.roc_curve <- function(x, ...) {
  sensitivity(attr(x, "params"), ...)
}


sensitivity.binormal_params <- function(x, specificity, ...) {
  sapply(specificity, function(spec) {
    fpf <- 1.0 - spec
    .Fortran("cvbmrocfpf2tpf",
             x$a, x$b,
             fpf = fpf, tpf = double(1), double(1))$tpf
  })
}


sensitivity.empirical_curve <- function(x, specificity, ties = max, ...) {
  approx(x$FPR, x$TPR, 1 - specificity, ties = ties)$y
}


sensitivity.empirical_params <- function(x, specificity, ties = max, ...) {
  curve <- roc_curves(x$truth, x$rating, method = "empirical")
  sensitivity(curve, specificity, ties = ties, ...)
}


sensitivity.proproc_params <- function(x, specificity, ...) {
  sapply(specificity, function(spec) {
    fpf <- 1.0 - spec
    .Fortran("pbmrocfpf2tpf",
             x$d_a, x$c,
             fpf = fpf, tpf = double(1), double(1))$tpf
  })
}


specificity <- function(x, ...) {
  UseMethod("specificity")
}


specificity.roc_curve <- function(x, ...) {
  specificity(attr(x, "params"), ...)
}


specificity.binormal_params <- function(x, sensitivity, ...) {
  sapply(sensitivity, function(sens) {
    1 - .Fortran("cvbmroctpf2fpf",
                 x$a, x$b,
                 tpf = as.double(sens), fpf = as.double(1), double(1))$fpf
  })
}


specificity.empirical_curve <- function(x, sensitivity, ties = max, ...) {
  1 - approx(x$TPR, x$FPR, sensitivity, ties = ties)$y
}


specificity.empirical_params <- function(x, sensitivity, ties = max, ...) {
  curve <- roc_curves(x$truth, x$rating, method = "empirical")
  specificity(curve, sensitivity, ties = ties, ...)
}


specificity.proproc_params <- function(x, sensitivity, ...) {
  sapply(sensitivity, function(sens) {
    1 - .Fortran("pbmroctpf2fpf",
                 x$d_a, x$c,
                 tpf = as.double(sens), fpf = as.double(1), double(1))$fpf
  })
}
