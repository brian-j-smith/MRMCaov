roc_counts <- function(truth, rating, collapse = FALSE, prior_count = 0) {
  events <- as.numeric(is_reference(truth))
  res <- rowsum(data.frame(events = events, nonevents = 1 - events), rating)
  res$cutoffs <- sort(unique(rating))

  if (prior_count) {
    prior <- 0.5 * prior_count / nrow(res)
    res$events <- res$events + prior
    res$nonevents <- res$nonevents + prior
  } else if (collapse) {
    res$events <- cumsum(res$events)
    res$nonevents <- cumsum(res$nonevents)
    drop <- (duplicated(res$events) & duplicated(res$events, fromLast = TRUE)) |
      (duplicated(res$nonevents) & duplicated(res$nonevents, fromLast = TRUE))
    res <- res[!drop, ]
    res$events <- c(res$events[1], diff(res$events))
    res$nonevents <- c(res$nonevents[1], diff(res$nonevents))
  }

  res
}


#' ROC Performance Curves
#'
#' Calculation of ROC curve true positive rate (TPR) and false positive rate
#' (FPR) pairs for values of a numeric rating of a true binary response.
#'
#' @rdname roc_curves
#'
#' @param truth vector of true binary statuses.
#' @param rating vector of numeric ratings.
#' @param groups list or data frame of grouping variables of the same lengths as
#'   \code{truth} and \code{rating}.
#' @param method character string indicating the curve type as
#'   \code{"binormal"}, \code{"binormalLR"}, \code{"empirical"}, or
#'   \code{"trapezoidal"} or the averaging of binormal curves over
#'   \code{"points"} or \code{"parameters"}.
#' @param x object returned by \code{\link{mrmc}} or \code{roc_curves} for which
#'   to compute points on or to average over the curves.
#' @param values numeric vector of values at which to compute ROC curve points,
#'   or \code{NULL} for default empirical values as determined by \code{which}.
#' @param which character string indicating whether to use curve-specific
#'  observed values and 0 and 1 (\code{"curve"}), the combination of these
#'  values over all curves (\code{"curves"}), or only the observed
#'  curve-specific values (\code{"observed"}).
#' @param metric reader performance metric to which the \code{values}
#'   correspond.
#' @param ties function determining empirical roc points returned in cases of
#'   ties.
#' @param ... arguments passed from the \code{mean()} method to \code{points()}.
#'
#' @return
#' Function \code{roc_curves} returns an \code{roc_curve} class object of a
#' single estimated ROC curve or an \code{roc_curves} class object of grouped
#' ROC curves.  Function \code{parameters} extracts the parameters that define
#' the curves, \code{points} returns a data frame of points on individual
#' curves, and \code{means} returns points on averaged curves (Chen and
#' Samuelson, 2014).
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
roc_curves.default <- function(
  truth, rating, groups = list(), method = "empirical", ...
) {
  choices <- c("binormal", "binormalLR", "empirical", "proproc", "trapezoidal")
  method <- match.arg(method, choices)

  method <- switch(method,
    "proproc" = "binormalLR",
    "trapezoidal" = "empirical",
    method
  )
  roc_curve <- get(paste0(method, "_curve"))

  data <- tibble(truth = as.factor(truth), rating = rating)
  if (length(groups)) {
    curves <- by(data, groups, function(split) {
      roc_curve(split, ...)
    }, simplify = FALSE)
    groups <- expand.grid(dimnames(curves))
    keep <- !is.null(curves)
    curves <- tibble(Group = groups[keep, , drop = FALSE], Curve = curves[keep])
    new_roc_curves(curves, method = method)
  } else {
    roc_curve(data, ...)
  }
}


roc_curves.roc_curve <- function(x, ...) {
  roc_curves(x$data$truth, x$data$rating, ...)
}


roc_curves.roc_curves <- function(x, method = "empirical", ...) {
  x$Curve <- lapply(x$Curve, function(curve) {
    roc_curves(curve, method = method, ...)
  })
  new_roc_curves(x, method = method)
}


new_roc_curves <- function(x, method) {
  x <- as_tibble(x)
  structure(x, class = c(paste0(method, "_curves"), "roc_curves", class(x)))
}


#' @rdname roc_curves
#'
roc_curves.mrmc <- function(x, ...) {
  roc_method <- unlist(strsplit(x$vars["metric"], "_"))[1]
  if (!(roc_method %in% c("binormal", "binormalLR", "proproc"))) {
    roc_method <- "empirical"
  }
  vars <- x$vars[c("reader", "test")]
  roc_curves(
    x$mrmc_data$truth, x$mrmc_data$rating,
    groups = structure(x$mrmc_data[names(vars)], names = vars),
    method = roc_method
  )
}


#' @rdname roc_curves
#'
roc_curves.stmc <- function(x, ...) {
  roc_method <- unlist(strsplit(x$metric, "_"))[1]
  if (!(roc_method %in% c("binormal", "binormalLR", "proproc"))) {
    roc_method <- "empirical"
  }
  roc_curves(x$stmc_data$truth, x$stmc_data$rating, method = roc_method)
}


binormal_curve <- function(data, ...) {
  params <- binormal_params(data$truth, data$rating, ...)
  structure(
    list(
      params = tibble(a = params$a, b = params$b),
      data = data
    ),
    class = c("binormal_curve", "roc_curve")
  )
}


binormalLR_curve <- function(data, ...) {
  params <- binormalLR_params(data$truth, data$rating, ...)
  structure(
    list(
      params = tibble(
        Metz = tibble(
          d_a = params$d_a,
          c = params$c
        ),
        bichisquared = tibble(
          lambda = ((1 - params$c) / (1 + params$c))^2,
          theta = (params$d_a * (1 + params$c) / (4 * params$c))^2 *
            (1 + params$c^2)
        ),
        binormal = tibble(
          a = params$d_a * sqrt(1 + params$c^2) / (1 - params$c),
          b = (1 + params$c) / (1 - params$c)
        )
      ),
      data = data
    ),
    class = c("binormalLR_curve", "roc_curve")
  )
}


empirical_curve <- function(data, prior_count = 0, ...) {
  counts <- roc_counts(data$truth, data$rating, prior_count = prior_count)
  f <- function(x) c(rev(cumsum(rev(x))) / sum(x), 0)
  structure(
    list(
      params = tibble(
        FPR = f(counts$nonevents),
        TPR = f(counts$events)
      ),
      data = data
    ),
    class = c("empirical_curve", "roc_curve")
  )
}


#' @rdname roc_curves
#'
parameters <- function(x, ...) {
  UseMethod("parameters")
}


#' @rdname roc_curves
#'
parameters.roc_curve <- function(x, ...) {
  x$params
}


#' @rdname roc_curves
#'
parameters.roc_curves <- function(x, ...) {
  params <- do.call(rbind, lapply(x$Curve, parameters))
  tibble(Group = x$Group, as_tibble(params))
}


parameters.empirical_curves <- function(x, ...) {
  params_list <- t(mapply(parameters, x$Curve))
  tibble(Group = x$Group, as_tibble(params_list))
}


#' @rdname roc_curves
#'
parameters.mrmc <- function(x, ...) {
  parameters(roc_curves(x))
}


#' @rdname roc_curves
#'
parameters.stmc <- function(x, ...) {
  parameters(roc_curves(x))
}


#' @rdname roc_curves
#' @aliases points
#'
points.roc_curve <- function(x, metric = c("specificity", "sensitivity"),
                             values = seq(0, 1, length = 101), ...) {
  metric <- match.arg(metric)
  switch(metric,
    "specificity" = {
      input_name <- "FPR"
      input <- sort(1 - values)
      output_name <- "TPR"
      output_fun <- function(x, input) sensitivity(x, 1 - input)
    },
    "sensitivity" = {
      input_name <- "TPR"
      input <- sort(values)
      output_name <- "FPR"
      output_fun <- function(x, input) 1 - specificity(x, input)
    }
  )

  output <- output_fun(x, input)
  new_pts <- tibble(input, output)
  names(new_pts) = c(input_name, output_name)

  structure(
    new_pts[c("FPR", "TPR")], metric = metric,
    class = c("roc_points", class(new_pts))
  )
}


#' @rdname roc_curves
#'
points.roc_curves <- function(
  x, metric = c("specificity", "sensitivity"), values = seq(0, 1, length = 101),
  ...
) {
  metric <- match.arg(metric)
  new_pts_list <- lapply(x$Curve, function(curve) {
    points(curve, metric = metric, values = values, ...)
  })
  new_pts <- curves2tibble(new_pts_list, x$Group)
  structure(
    new_pts, metric = metric,
    class = c("roc_points", class(new_pts))
  )
}


#' @rdname roc_curves
#'
points.empirical_curve <- function(
  x, metric = c("specificity", "sensitivity"), values = NULL,
  which = c("curve", "curves", "observed"), ties = max, ...
) {

  metric <- match.arg(metric)
  which <- match.arg(which)

  new_pts <- if (is.numeric(values)) {

    switch(metric,
      "specificity" = {
        input_name <- "FPR"
        output_name <- "TPR"
        output_fun <- function(x, input, ...) {
          sensitivity(x, 1 - input, ...)
        }
        sort_metric <- function(x) sort(1 - x)
      },
      "sensitivity" = {
        input_name <- "TPR"
        output_name <- "FPR"
        output_fun <- function(x, input, ...) {
          1 - specificity(x, input, ...)
        }
        sort_metric <- function(x) sort(x)
      }
    )

    input <- sort_metric(values)
    output <- output_fun(x, input, ties = ties)
    pts <- tibble(input, output)
    names(pts) <- c(input_name, output_name)
    pts[c("FPR", "TPR")]

  } else {

    params <- parameters(x)
    if (which == "observed") params <- params[-c(1, nrow(params)), ]
    params

  }

  structure(new_pts, metric = metric, class = c("roc_points", class(new_pts)))

}


#' @rdname roc_curves
#'
points.empirical_curves <- function(
  x, metric = c("specificity", "sensitivity"), values = NULL,
  which = c("curve", "curves", "observed"), ties = max, ...
) {

  metric <- match.arg(metric)
  which <- match.arg(which)

  new_pts_list <- if (is.null(values) && which == "curves") {

    switch(metric,
      "specificity" = {
        input_name <- "FPR"
        pts_fun <- function(x, input, ...) {
          points(x, values = 1 - input, ...)
        }
      },
      "sensitivity" = {
        input_name <- "TPR"
        pts_fun <- function(x, input, ...) points(x, values = input, ...)
      }
    )
    input_list <- parameters(x)[[input_name]]
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
      points(
        curve, metric = metric, values = values, which = which, ties = ties,
        ...
      )
    })

  }
  new_pts <- curves2tibble(new_pts_list, x$Group)

  structure(new_pts, metric = metric, class = c("roc_points", class(new_pts)))

}


#' @rdname roc_curves
#' @aliases mean
#'
mean.roc_curve <- function(x, ...) {
  points(x, ...)
}


#' @rdname roc_curves
#'
mean.roc_curves <- function(x, ...) {
  pts <- points(x, which = "curves", ...)
  metric <- attr(pts, "metric")
  switch(metric,
    "specificity" = {
      input_name <- "FPR"
      output_name <- "TPR"
    },
    "sensitivity" = {
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
    auc_mean <- mean(sapply(x$Curve, auc))
    b_mean <- mean(parameters(x$Curve)$b)
    a <- sqrt(1 + b_mean^2) * qnorm(auc_mean)
    curve <- structure(
      params = list(a = a, b = b_mean), class = class(x$Curve[[1]])
    )
    mean(curve, ...)
  } else NextMethod()
}


auc <- function(x, ...) {
  UseMethod("auc")
}


auc.binormal_curve <- function(
  x, partial = FALSE, min = 0, max = 1, normalize = FALSE, ...
) {
  params <- parameters(x)
  if (isFALSE(partial)) {
    .auc.binormal(params$a, params$b)
  } else {
    partial <- match.arg(partial, c("sensitivity", "specificity"))
    normalize <- if (normalize) max - min else 1
    if (partial == "sensitivity") {
      params <- list(a = params$a / params$b, b = 1 / params$b)
    }
    x <- qnorm(1 - c(max, min))
    rho <- -params$b / sqrt(1 + params$b^2)
    R <- rho + diag(1 - rho, 2)
    y <- -rho * params$a / params$b
    auc1 <- c(pmvnorm(upper = c(x[1], y), corr = R))
    auc2 <- c(pmvnorm(upper = c(x[2], y), corr = R))
    (auc2 - auc1) / normalize
  }
}


.auc.binormal <- function(a, b) {
  pnorm(a / sqrt(1 + b^2))
}


auc.binormalLR_curve <- function(
  x, partial = FALSE, min = 0, max = 1, normalize = FALSE, ...
) {
  params <- parameters(x)$Metz
  d_a <- params$d_a
  c <- params$c

  if (has_binormal_equiv(x)) {
    call <- match.call()
    call[[1]] <- quote(auc)
    call$x <- as.binormal_curve(x)
    eval.parent(call)
  } else if (isFALSE(partial)) {
    .auc.binormalLR(d_a, c)
  } else {
    partial <- match.arg(partial, c("sensitivity", "specificity"))
    normalize <- if (normalize) max - min else 1
    if (partial == "sensitivity") {
      c <- -c
    }

    y <- numeric(2)
    y[1] <- d_a / sqrt(2)
    y[2] <- y[1] * (1 + c^2) / (2 * c)
    rho <- -(1 + c) / sqrt(2 * (1 + c^2))
    R <- rho + (1 - rho) * diag(2)

    cutoffs <- numeric(2)
    cutoffs[1] <- binormalLR_cutoff(1 - max, d_a, c)
    cutoffs[2] <- binormalLR_cutoff(1 - min, d_a, c, init = cutoffs[1])

    x <- binormalLR_z1_fpr(cutoffs, d_a, c)
    xp <- binormalLR_z2_fpr(cutoffs, d_a, c)

    cdf <- function(x, y) pmvnorm(upper = c(x, y), corr = R)
    F1 <- cdf(x[2], y[1])
    F2 <- cdf(x[1], y[1])
    F3 <- cdf(x[2], y[2])
    F4 <- cdf(x[1], y[2])
    F5 <- cdf(xp[2], -y[2])
    F6 <- cdf(xp[1], -y[2])
    F7 <- cdf(xp[2], -y[1])
    F8 <- cdf(xp[1], -y[1])

    res <- c((F1 + F3 + F5 + F7) - (F2 + F4 + F6 + F8))
    if (c > 0) res <- res - max + min
    res / normalize
  }
}


.auc.binormalLR <- function(d_a, c) {
  z <- d_a / sqrt(2)
  rho <- -(1 - c^2) / (1 + c^2)
  R <- rho + (1 - rho) * diag(2)
  c(pnorm(z) + 2 * pmvnorm(upper = c(-z, 0), corr = R))
}


auc.empirical_curve <- function(
  x, partial = FALSE, min = 0, max = 1, normalize = FALSE, ...
) {
  x_coords <- rev(x$params$FPR)
  y_coords <- rev(x$params$TPR)
  if (!isFALSE(partial)) {
    partial <- match.arg(partial, c("sensitivity", "specificity"))
    switch(partial,
      "sensitivity" = {
        spec <- 1 - x_coords
        x_coords <- y_coords
        y_coords <- spec
        bounds <- data.frame(x = c(min, max))
        bounds$y <- specificity(x, bounds$x)
      },
      "specificity" = {
        bounds <- c(max, min)
        bounds <- data.frame(x = 1 - bounds, y = sensitivity(x, bounds))
      }
    )
    keep <- x_coords > bounds$x[1] & x_coords < bounds$x[2]
    x_coords <- c(bounds$x[1], x_coords[keep], bounds$x[2])
    y_coords <- c(bounds$y[1], y_coords[keep], bounds$y[2])
  }
  res <- sum(diff(x_coords) * (y_coords[-length(y_coords)] + y_coords[-1]) / 2)
  if (normalize) res <- res / diff(range(x_coords))
  res
}


roc_eu <- function(x, ...) {
  UseMethod("roc_eu")
}


roc_eu.roc_curve <- function(x, slope = 1, ...) {
  f <- function(sens) {
    spec <- specificity(x, sens)
    sens - slope * (1 - spec)
  }
  optimize(f, c(0, 1), maximum = TRUE)$objective
}


roc_eu.empirical_curve <- function(x, slope = 1, ...) {
  pts <- points(x)
  max(pts$TPR - slope * pts$FPR)
}


sensitivity <- function(x, ...) {
  UseMethod("sensitivity")
}


sensitivity.binormal_curve <- function(x, specificity, ...) {
  params <- parameters(x)
  ifelse(specificity %in% c(0, 1),
    1 - specificity,
    pnorm(params$a + params$b * qnorm(1.0 - specificity))
  )
}


sensitivity.binormalLR_curve <- function(x, specificity, ...) {
  if (has_binormal_equiv(x)) {
    call <- match.call()
    call[[1]] <- quote(sensitivity)
    call$x <- as.binormal_curve(x)
    eval.parent(call)
  } else {
    params <- parameters(x)$bichisquared
    ncp <- c(params$theta, params$lambda * params$theta)
    if (all(ncp <= 1e4)) {
      f <- function(p) {
        pchisq(qchisq(p, 1, ncp[1]) / params$lambda, 1, ncp[2])
      }
      if (params$lambda < 1) f(1 - specificity) else 1 - f(specificity)
    } else {
      params <- parameters(x)$Metz
      cutoff <- binormalLR_cutoff(1 - specificity, params$d_a, params$c)
      binormalLR_tpr(cutoff, params$d_a, params$c)
    }
  }
}


sensitivity.empirical_curve <- function(x, specificity, ties = max, ...) {
  params <- parameters(x)
  approx(params$FPR, params$TPR, 1 - specificity, ties = ties)$y
}


specificity <- function(x, ...) {
  UseMethod("specificity")
}


specificity.binormal_curve <- function(x, sensitivity, ...) {
  params <- parameters(x)
  ifelse(sensitivity %in% c(0, 1),
    1 - sensitivity,
    1 - pnorm((qnorm(sensitivity) - params$a) / params$b)
  )
}


specificity.binormalLR_curve <- function(x, sensitivity, ...) {
  if (has_binormal_equiv(x)) {
    call <- match.call()
    call[[1]] <- quote(specificity)
    call$x <- as.binormal_curve(x)
    eval.parent(call)
  } else {
    params <- parameters(x)$bichisquared
    ncp <- c(params$lambda * params$theta, params$theta)
    if (all(ncp <= 1e4)) {
      f <- function(p) {
        pchisq(qchisq(p, 1, ncp[1]) * params$lambda, 1, ncp[2])
      }
      if (params$lambda < 1) 1  - f(sensitivity) else f(1 - sensitivity)
    } else {
      params <- parameters(x)$Metz
      cutoff <- binormalLR_cutoff(1 - sensitivity, params$d_a, -params$c)
      1 - binormalLR_fpr(cutoff, params$d_a, params$c)
    }
  }
}


specificity.empirical_curve <- function(x, sensitivity, ties = max, ...) {
  params <- parameters(x)
  1 - approx(params$TPR, params$FPR, sensitivity, ties = ties)$y
}
