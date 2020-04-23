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
#'   \code{"proproc"} or indicating the mean estimation method.
#' @param x object returned by \code{roc_curves} for which to compute points on
#'   or to average over the curves.
#' @param values numeric vector of values at which to compute the points.  If
#'   \code{NULL} then the intersection of all empirical points is used.
#' @param metric reader performance metric to which the \code{values}
#'   correspond.
#' @param ... arguments passed from the \code{mean()} method to \code{points()}.
#'
#' @seealso \code{\link{plot}}
#'
#' @references
#' Chen W and Samuelson FW (2014) The average receiver operating characteristics
#' curve in multireader multicase imaging studies. \emph{British Journal of
#' Radiology}, 87(1040) doi: 10.1259/bjr.20140016.
#'
#' @examples
#' curves <- with(VanDyke,
#'   roc_curves(truth, rating, groups = list(Test = treatment, Reader = reader))
#' )
#' points(curves)
#' mean(curves)
#'
roc_curves <- function(truth, rating, groups = list(), method = "empirical") {
  method <- match.arg(method,
                      c("binormal", "empirical", "trapezoidal", "proproc"))
  f <- switch(method,
              "binormal" = param_curves,
              "empirical" = empirical_curves,
              "trapezoidal" = empirical_curves,
              "proproc" = param_curves)
  f(truth, rating, groups, method = method)
}


empirical_curves <- function(truth, rating, groups, ...) {
  data <- tibble(truth = truth, rating = rating)

  get_pts <- function(data) {
    roc <- pROC::roc(data$truth, data$rating, auc = FALSE, quiet = TRUE)
    pts <- tibble(FPR = 1 - roc$specificities, TPR = roc$sensitivities)
    pts[nrow(pts):1, ]
  }

  curves <- if (length(groups)) {
    pts_list <- by(data, groups, get_pts)
    pts <- do.call(rbind, pts_list)

    labels <- expand.grid(dimnames(pts_list))
    labels_ind <- rep(1:nrow(labels), times = sapply(pts_list, nrow))

    tibble(Group = labels[labels_ind, , drop = FALSE],
           FPR = pts$FPR, TPR = pts$TPR)
  } else {
    get_pts(data)
  }

  structure(curves, class = c("empirical_curves", "roc_curves", class(curves)))
}


param_curves <- function(truth, rating, groups, method, ...) {
  roc_params <- switch(method,
                       "binormal" = binormal_params,
                       "proproc" = proproc_params)

  data <- tibble(truth = truth, rating = rating)

  get_curve <- function(data) {
    roc <- pROC::roc(data$truth, data$rating, auc = FALSE, quiet = TRUE)
    spec <- rev(roc$specificities)
    params <- roc_params(data$truth, data$rating)
    pts <- tibble(FPR = 1 - spec, TPR = sensitivity(params, spec))
    list(pts = pts, params = params)
  }

  curves <- if (length(groups)) {
    curve_list <- by(data, groups, get_curve)
    pts <- do.call(rbind, lapply(curve_list, getElement, name = "pts"))
    params <- lapply(curve_list, getElement, name = "params")

    labels <- expand.grid(dimnames(curve_list))
    labels_ind <- rep(1:nrow(labels),
                      times = sapply(curve_list, function(x) nrow(x$pts)))

    structure(tibble(Group = labels[labels_ind, , drop = FALSE],
                     FPR = pts$FPR, TPR = pts$TPR),
              models = tibble(Group = labels, Params = params))
  } else {
    curve <- get_curve(data)
    structure(curve$pts, models = tibble(Params = list(curve$params)))
  }

  structure(curves, class = c(paste0(method, "_curves"), "param_curves",
                              "roc_curves", class(curves)))
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
  structure(list(a = roc$a, b = roc$b), class = "binormal_params")
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
  structure(list(d_a = roc$d_a, c = roc$c), class = "proproc_params")
}


#' @rdname roc_curves
#'
points.empirical_curves <- function(x, values = NULL,
                                    metric = c("specificity", "sensitivity"),
                                    ...) {
  metric <- match.arg(metric)
  switch(metric,
         "specificity" = {
           input_name <- "FPR"
           output_name <- "TPR"
           sort_metric <- function(x) sort(1 - x)
         },
         sensitivity = {
           input_name <- "TPR"
           output_name <- "FPR"
           sort_metric <- function(x) sort(x)
         })

  new_pts <- if (!is.null(x[["Group"]])) {

    new_pts_list <- if (is.null(values)) {
      input <- sort(unique(x[[input_name]]))
      input_dups_list <- tapply(x[[input_name]], x[["Group"]], function(input) {
        input[duplicated(input)]
      })
      input_dups <- sort(unique(unlist(input_dups_list)))
      n <- length(input) + length(input_dups)
      by(x, x[["Group"]], function(curve) {
        output <- approx(curve[[input_name]], curve[[output_name]], input,
                         ties = "ordered")$y
        output_dups <- approx(curve[[input_name]], curve[[output_name]],
                              input_dups, ties = min)$y
        pts <- tibble(c(input, input_dups), c(output, output_dups))
        names(pts) <- c(input_name, output_name)
        pts[order(pts$FPR, pts$TPR), c("FPR", "TPR")]
      })

    } else {

      input <- sort_metric(values)
      n <- length(input)
      by(x, x[["Group"]], function(curve) {
        output <- approx(curve[[input_name]], curve[[output_name]], input,
                         ties = "ordered")$y
        pts <- tibble(input, output)
        names(pts) = c(input_name, output_name)
        pts[c("FPR", "TPR")]
      })
    }
    new_pts <- do.call(rbind, new_pts_list)

    labels <- expand.grid(dimnames(new_pts_list))
    labels_ind <- rep(seq(nrow(labels)), each = n)

    tibble(Group = labels[labels_ind, , drop = FALSE],
           FPR = new_pts$FPR, TPR = new_pts$TPR)

  } else if (!is.null(values)) {

    input <- sort_metric(values)
    output <- approx(x[[input_name]], x[[output_name]], input,
                     ties = "ordered")$y
    pts <- tibble(input, output)
    names(pts) <- c(input_name, output_name)
    pts[c("FPR", "TPR")]

  } else as_tibble(x)

  structure(new_pts, metric = metric, class = c("roc_points", class(new_pts)))
}


#' @rdname roc_curves
#'
points.param_curves <- function(x, values = seq(0, 1, length = 100),
                                metric = c("specificity", "sensitivity"), ...) {
  metric <- match.arg(metric)
  switch(metric,
         specificity = {
           input_name <- "FPR"
           input <- sort(1 - values)
           output_name <- "TPR"
           output_fun <- function(params, x) sensitivity(params, 1 - x)
         },
         sensitivity = {
           input_name <- "TPR"
           input <- sort(values)
           output_name <- "FPR"
           output_fun <- function(params, x) 1 - specificity(params, x)
         })

  n <- length(input)

  models <- attr(x, "models")
  new_pts_list <- lapply(models$Params, function(params) {
    output <- output_fun(params, input)
    pts <-tibble(input, output)
    names(pts) = c(input_name, output_name)
    pts[c("FPR", "TPR")]
  })
  new_pts <- do.call(rbind, new_pts_list)
  if (!is.null(models[["Group"]])) {
    inds <- rep(seq(new_pts_list), each = n)
    new_pts <- tibble(Group = models[["Group"]][inds, ],
                      FPR = new_pts$FPR, TPR = new_pts$TPR)
  }

  structure(new_pts, metric = metric, class = c("roc_points", class(new_pts)))
}


#' @rdname roc_curves
#'
mean.roc_curves <- function(x, ...) {
  pts <- points(x, ...)
  if (!is.null(pts[["Group"]])) {
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
  } else pts
}


#' @rdname roc_curves
#'
mean.binormal_curves <- function(x, method = c("parameters", "points"), ...) {
  if (match.arg(method) == "parameters") {
    params <- attr(x, "models")$Params
    auc_mean <- mean(sapply(params, auc))
    b_mean <- mean(mapply(getElement, params, "b"))
    a <- sqrt(1 + b_mean^2) * qnorm(auc_mean)
    params <- structure(list(a = a, b = b_mean), class = "binormal_params")
    x <- structure(tibble(), models = tibble(Params = list(params)),
                   class = class(x))
  }
  NextMethod()
}


auc <- function(x, ...) {
  UseMethod("auc")
}


auc.binormal_params <- function(x, partial = FALSE, min = 0, max = 1, ...) {
  if (isFALSE(partial)) {
    pnorm(x$a / sqrt(1 + x$b^2))
  } else {
    partial <- partial_params(partial, min, max)
    .Fortran("cvbmrocpartial",
             params$a, params$b,
             partial$min, partial$max, partial$flag,
             est = double(1), err = integer(1))$est
  }
}


auc.proproc_params <- function(x, partial = FALSE, min = 0, max = 1, ...) {
  fit <- if (isFALSE(partial)) {
    rho <- -1 * (1 - x$c^2) / (1 + x$c^2)
    rho <- rbind(c(1, rho), c(rho, 1))
    rho <- diag(2)
    rho[rbind(c(1, 2), c(2, 1))] <- -1 * (1 - x$c^2) / (1 + x$c^2)
    pnorm(x$d_a / sqrt(2)) +
      2 * pmvnorm(upper = c(-x$d_a / sqrt(2), 0), corr = rho)
  } else {
    partial <- partial_params(partial, min, max)
    .Fortran("pbmrocpartial",
             params$d_a, params$c,
             partial$min, partial$max, partial$flag,
             est = double(1), err = integer(1))$est
  }
}


partial_params <- function(x, min, max) {
  x <- match.arg(partial, c("sensitivity", "specificity"))
  min <- as.double(min)
  max <- as.double(max)
  if (x == "specificity") {
    flag <- 1L
    min <- 1 - max
    max <- 1 - min
  } else {
    flag <- 2L
  }
  list(min, max, flag)
}


sensitivity <- function(x, ...) {
  UseMethod("sensitivity")
}


sensitivity.binormal_params <- function(x, specificity, ...) {
  sapply(specificity, function(spec) {
    fpf <- 1.0 - spec
    .Fortran("cvbmrocfpf2tpf",
             x$a, x$b,
             fpf = fpf, tpf = double(1), double(1))$tpf
  })
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


specificity.binormal_params <- function(x, sensitivity, ...) {
  sapply(sensitivity, function(sens) {
    1 - .Fortran("cvbmroctpf2fpf",
                 x$a, x$b,
                 tpf = as.double(sens), fpf = as.double(1), double(1))$fpf
  })
}


specificity.proproc_params <- function(x, sensitivity, ...) {
  sapply(sensitivity, function(sens) {
    1 - .Fortran("pbmroctpf2fpf",
                 x$d_a, x$c,
                 tpf = as.double(sens), fpf = as.double(1), double(1))$fpf
  })
}
