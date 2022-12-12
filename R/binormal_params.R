binormal_params <- function(
  truth, rating, prior_count = 0, control = list(), ...
) {
  counts <- roc_counts(
    truth, rating, collapse = TRUE, prior_count = prior_count
  )
  inits <- binormal_inits(
    truth, rating, counts$cutoffs, prior_count = prior_count
  )
  res <- inits[c("a", "b", "cutoffs")]

  if (!inits$optim) {
    res$cutoffs <- binormal_cutoffs(
      res, counts$events, counts$nonevents, method = control$init_cutoffs
    )
    res <- binormal_mle(res, counts$events, counts$nonevents, control = control)
  }

  res
}


binormal_inits <- function(
  truth, rating, cutoffs, intercept = TRUE, shift = 0, scale = 1,
  prior_count = 0
) {
  cutoffs <- 3 * (cutoffs - min(cutoffs)) / diff(range(cutoffs)) - 1
  res <- list(
    a = NA,
    b = 1,
    cutoffs = cutoffs[-1] - diff(cutoffs) / 2,
    optim = TRUE
  )

  transform <- function(x) lapply(x, function(col) (col + shift) / scale)
  make_data <- function(x) data.frame(x = qnorm(x$FPR), y = qnorm(x$TPR))

  emp_curve <- roc_curves(truth, rating, prior_count = prior_count)
  emp_pts <- transform(points(emp_curve))
  emp_data <- make_data(emp_pts)

  if (nrow(emp_data) <= 2) {
    res$a <- 0
  } else if (all(emp_data$x == -Inf | emp_data$y == Inf)) {
    res$a <- Inf
  } else if (all(emp_data$x == Inf | emp_data$y == -Inf)) {
    res$a <- -Inf
  } else {
    res$optim <- FALSE
  }

  if (!res$optim) {
    props <- (1:2) / 3
    emp_pts_fpr <- get_emp_points(emp_pts$FPR, emp_pts$TPR, props)
    emp_pts_tpr <- get_emp_points(emp_pts$TPR, emp_pts$FPR, props)
    emp_pts <- transform(unique(tibble(
      FPR = c(emp_pts_fpr$x, emp_pts_tpr$y),
      TPR = c(emp_pts_fpr$y, emp_pts_tpr$x)
    )))
    emp_data <- rbind(emp_data, make_data(emp_pts))
    keep <- is.finite(emp_data$x) & is.finite(emp_data$y)

    if (diff(range(emp_data$x[keep])) == 0) {
      ind <- match(Inf, emp_data$y)
      res$b <- 1e4
      res$a <- -emp_data$x[ind] * res$b
    } else {
      fo <- reformulate("x", "y", intercept = intercept)
      emp_coef <- coef(lm(fo, data = emp_data[keep, ]))
      res$a <- emp_coef["(Intercept)"]
      res$b <- emp_coef["x"]
      if (is.na(res$a)) res$a <- 0
      if (is.na(res$b)) res$b <- 1
    }
  }

  res
}


binormal_mle <- function(inits, events, nonevents, control = list()) {
  control <- do.call(
    function(
      method = c("trust", "CG", "BFGS", "L-BFGS-B", "Nelder-Mead"),
      gradient = TRUE, ...
    ) {
      list(method = match.arg(method), gradient = gradient)
    },
    control
  )

  inits <- c(inits$a, inits$b, inits$cutoffs)
  switch(control$method,
    "trust" = {
      control$gradient <- TRUE
      res <- trust::trust(
        function(x) {
          res <- binormal_objfun(
            x, events, nonevents, all = TRUE, loglik_min = -Inf
          )
          list(
            value = -res$loglik,
            gradient = -res$gradient,
            hessian = -res$hessian
          )
        },
        inits,
        1,
        10
      )
      estimates <- res$argument
    },
    "L-BFGS-B" = {
      num_cutoffs <- length(inits) - 2
      res <- optim(
        inits,
        fn = function(x) -binormal_loglik(x, events, nonevents),
        gr = if (control$gradient) {
          function(x) -binormal_gradient(x, events, nonevents)
        },
        method = control$method,
        lower = c(-Inf, 0, rep(-Inf, num_cutoffs)),
        upper = Inf
      )
      estimates <- res$par
    },
    {
      res <- optim(
        inits,
        fn = function(x) -binormal_loglik(x, events, nonevents),
        gr = if (control$gradient) {
          function(x) -binormal_gradient(x, events, nonevents)
        },
        method = control$method
      )
      estimates <- res$par
    }
  )
  res$control <- control

  list(
    a = estimates[1],
    b = estimates[2],
    cutoffs = estimates[-(1:2)],
    inits = inits,
    optim = res
  )
}


binormal_cutoffs <- function(
  inits, events, nonevents, method = c("uni", "multi")
) {
  switch(match.arg(method),
    "multi" = binormal_cutoffs_multi(inits, events, nonevents),
    "uni" = binormal_cutoffs_uni(inits, events, nonevents)
  )
}


binormal_cutoffs_multi <- function(inits, events, nonevents) {
  fixed <- c(inits$a, inits$b)
  optim(
    inits$cutoffs,
    fn = function(x) -binormal_loglik(x, events, nonevents, fixed = fixed),
    gr = function(x) -binormal_gradient(x, events, nonevents, fixed = fixed),
    method = "BFGS"
  )$par
}


binormal_cutoffs_uni <- function(inits, events, nonevents) {
  cum_events <- cumsum(events)
  tot_events <- tail(cum_events, 1)
  cum_nonevents <- cumsum(nonevents)
  tot_nonevents <- tail(cum_nonevents, 1)

  get_counts2 <- function(cum_events, cum_nonevents) {
    tibble(
      events = c(cum_events, tot_events - cum_events),
      nonevents = c(cum_nonevents, tot_nonevents - cum_nonevents)
    )
  }

  loglik2 <- function(cutoff) {
    binormal_loglik(
      cutoff, counts2$events, counts2$nonevents, fixed = c(inits$a, inits$b)
    )
  }

  res <- numeric(length(events) - 1)

  lower <- max(-5, (inits$a - 5) / inits$b)
  scale <- if (inits$b > 1) 1 / inits$b else 1
  finc <- function(x, scale, iter) x + scale * (1 + log(iter))

  for (ind in seq_along(res)) {
    counts2 <- get_counts2(cum_events[ind], cum_nonevents[ind])
    upper <- get_next_bound(lower, scale, loglik2, finc)
    lower <- get_next_cutoff(c(lower, upper), loglik2)
    res[ind] <- lower
  }

  res
}


binormal_objfun <- function(
  x, events, nonevents, fixed = NULL, all = FALSE, loglik = all, gradient = all,
  hessian = all, loglik_min = -1e100
) {
  res <- list()

  x <- c(fixed, x)
  a <- x[1]
  b <- x[2]
  cutoffs <- x[-(1:2)]

  in_support <- b > 0 && all(diff(cutoffs) > 0)

  if (!in_support) {
    n <- length(x) - length(fixed)
    if (loglik) res$loglik <- loglik_min
    if (gradient) res$gradient <- numeric(n)
    if (hessian) res$hessian <- matrix(0, n, n)
    return(res)
  }

  p <- binormal_probs(binormal_fpr, cutoffs, a, b)
  q <- binormal_probs(binormal_tpr, cutoffs, a, b)

  if (loglik) {
    res$loglik <- max(c(nonevents %*% log(p) + events %*% log(q)), loglik_min)
  }

  if (gradient || hessian) {
    keep <- seq(length(fixed) + 1L, length(x))

    dfpr_dcut <- -dnorm(-cutoffs)

    dtpr_da <- dnorm(a - b * cutoffs)
    dtpr_db <- dtpr_da * -cutoffs
    dtpr_dcut <- dtpr_da * -b
  }

  if (gradient) { 
    f <- function(counts, probs, dxr_da, dxr_db, dxr_dcut) {
      n <- length(counts)
      a <- -counts[-n] / probs[-n] + counts[-1] / probs[-1]
      c(
        if (length(dxr_da)) a %*% dxr_da else 0,
        if (length(dxr_db)) a %*% dxr_db else 0,
        a * dxr_dcut
      )[keep]
    }
    res$gradient <- f(nonevents, p, NULL, NULL, dfpr_dcut) +
      f(events, q, dtpr_da, dtpr_db, dtpr_dcut)
  }

  if (hessian) {
    f <- function(counts, probs, dxr_da, dxr_db, dxr_dcut) {
      num_cutoffs <- length(dxr_dcut)
      A <- matrix(0, num_cutoffs + 1, num_cutoffs + 2)
      if (length(dxr_da)) A[, 1] <- c(0, dxr_da) - c(dxr_da, 0)
      if (length(dxr_db)) A[, 2] <- c(0, dxr_db) - c(dxr_db, 0)
      inds1 <- seq_len(num_cutoffs)
      inds2 <- inds1 + 2
      A[cbind(inds1, inds2)] <- -dxr_dcut
      A[cbind(inds1 + 1, inds2)] <- dxr_dcut
      crossprod((-sum(counts) / probs) * A, A)[keep, keep]
    }
    res$hessian <- f(nonevents, p, NULL, NULL, dfpr_dcut) +
      f(events, q, dtpr_da, dtpr_db, dtpr_dcut)
  }

  res
}


binormal_loglik <- function(x, events, nonevents, fixed = NULL) {
  binormal_objfun(x, events, nonevents, fixed = fixed, loglik = TRUE)$loglik
}


binormal_gradient <- function(x, events, nonevents, fixed = NULL) {
  binormal_objfun(x, events, nonevents, fixed = fixed, gradient = TRUE)$gradient
}


binormal_probs <- function(f, cutoffs, a, b, min = 1e-10) {
  pmax(-diff(c(1, f(cutoffs, a, b), 0)), min)
}


binormal_fpr <- function(cutoff, a, b) {
  pnorm(-cutoff)
}


binormal_tpr <- function(cutoff, a, b) {
  pnorm(a - b * cutoff)
}


get_emp_points <- function(x, y, props) {
  data <- tibble(x = x, y = y)
  data <- data[order(data$x, data$y), ]
  props <- 1 - props
  tibble(
    x = c(data$x[-1] - diff(data$x) %o% props),
    y = c(data$y[-1] - diff(data$y) %o% props)
  )
}


get_next_bound <- function(
  x, scale, loglik2, finc1, finc2 = finc1, maxiter = 500
) {
  ll <- loglik2(x)
  old_bound <- x + scale
  away <- TRUE
  for (iter in 1:maxiter) {
    old_ll <- loglik2(old_bound)
    if (old_ll > ll) {
      break
    } else if (is_approx(old_ll, ll) && away) {
      old_bound <- finc1(old_bound, scale, iter)
    } else {
      old_bound <- (old_bound + x) / 2
      if (is_approx(old_bound, x)) {
        break
      }
      away <- FALSE
    }
  }

  for (iter in 1:maxiter) {
    bound <- finc2(old_bound, scale, iter)
    ll <- loglik2(bound)
    if (ll < old_ll) {
      break
    } else {
      old_bound <- bound
      old_ll <- ll
    }
  }

  bound
}


get_next_cutoff <- function(interval, loglik2) {
  optimize(
    function(x) -loglik2(x),
    interval = interval
  )$minimum
}
