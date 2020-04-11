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
#' @param method character string \code{"empirical"}, \code{"trapezoidal"}, or
#'   \code{"proproc"} indicating the curve type.
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
#' @examples
#' curves <- with(VanDyke,
#'   roc_curves(truth, rating, groups = list(Test = treatment, Reader = reader))
#' )
#' points(curves)
#' mean(curves)
#'
roc_curves <- function(truth, rating, groups = list(), method = "empirical") {
  method <- match.arg(method, c("empirical", "trapezoidal", "proproc"))
  switch(method,
         "empirical" = empirical_curves(truth, rating, groups),
         "trapezoidal" = empirical_curves(truth, rating, groups),
         "proproc" = proproc_curves(truth, rating, groups))
}


empirical_curves <- function(truth, rating, groups) {
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


proproc_curves <- function(truth, rating, groups) {
  data <- tibble(truth = truth, rating = rating)

  get_curve <- function(data) {
    roc <- pROC::roc(data$truth, data$rating, auc = FALSE, quiet = TRUE)
    spec <- rev(roc$specificities)
    params <- proproc_params(data$truth, data$rating)
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

  structure(curves, class = c("proproc_curves", "roc_curves", class(curves)))
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
points.empirical_curves <- function(x, values = NULL, metric = "specificity",
                                    ...) {
  metric <- match.arg(metric)

  new_pts <- if (!is.null(x[["Group"]])) {

    new_pts_list <- if (is.null(values)) {
      fpr <- sort(unique(x$FPR))
      fpr_dups_list <- tapply(x$FPR, x$Group, function(fpr) {
        fpr[duplicated(fpr)]
      })
      fpr_dups <- sort(unique(unlist(fpr_dups_list)))
      n <- length(fpr) + length(fpr_dups)
      by(x, x$Group, function(curve) {
        tpr <- approx(curve$FPR, curve$TPR, fpr, ties = "ordered")$y
        tpr_dups <- approx(curve$FPR, curve$TPR, fpr_dups, ties = min)$y
        pts <- tibble(FPR = c(fpr, fpr_dups), TPR = c(tpr, tpr_dups))
        pts[order(pts$FPR, pts$TPR), ]
      })
    } else {
      fpr <- sort(1 - values)
      n <- length(fpr)
      by(x, x$Group, function(curve) {
        tpr <- approx(curve$FPR, curve$TPR, fpr, ties = "ordered")$y
        tibble(FPR = fpr, TPR = tpr)
      })
    }
    new_pts <- do.call(rbind, new_pts_list)

    labels <- expand.grid(dimnames(new_pts_list))
    labels_ind <- rep(seq(nrow(labels)), each = n)

    tibble(Group = labels[labels_ind, , drop = FALSE],
           FPR = new_pts$FPR, TPR = new_pts$TPR)

  } else if (!is.null(values)) {

    fpr <- sort(1 - values)
    tpr <- approx(x$FPR, x$TPR, fpr, ties = "ordered")$y
    tibble(FPR = fpr, TPR = tpr)

  } else as_tibble(x)

  structure(new_pts, metric = metric, class = c("roc_points", class(new_pts)))
}


#' @rdname roc_curves
#'
points.proproc_curves <- function(x, values = seq(0, 1, length = 100),
                                  metric = "specificity", ...) {
  metric <- match.arg(metric)

  fpr <- sort(1 - values)
  n <- length(fpr)

  models <- attr(x, "models")
  new_pts_list <- lapply(models$Params, function(params) {
    tibble(FPR = fpr, TPR = sensitivity(params, 1 - fpr))
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
    tpr_list <- split(pts$TPR, pts[["Group"]])
    tpr <- rowMeans(do.call(cbind, tpr_list))
    structure(tibble(FPR = head(pts$FPR, n = length(tpr)), TPR = tpr),
              metric = attr(pts, "metric"), class = class(pts))
  } else pts
}


sensitivity <- function(x, ...) {
  UseMethod("sensitivity")
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


specificity.proproc_params <- function(x, sensitivity, ...) {
  sapply(sensitivity, function(sens) {
    1 - .Fortran("pbmroctpf2fpf",
                 x$d_a, x$c,
                 tpf = as.double(sens), fpf = as.double(1), double(1))$fpf
  })
}
