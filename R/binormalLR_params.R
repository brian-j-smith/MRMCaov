binormalLR_params <- function(
  truth, rating, prior_count = 0, control = list(), ...
) {
  convert_inits <- function(x) {
    list(
      d_a = max(x$a * sqrt(2 / (1 + x$b^2)), 0),
      c = (x$b - 1) / (x$b + 1),
      cutoffs = x$cutoffs
    )
  }

  counts <- roc_counts(
    truth, rating, collapse = TRUE, prior_count = prior_count
  )
  inits <- binormal_inits(
    truth, rating, counts$cutoffs, prior_count = prior_count
  )
  res <- convert_inits(inits)

  if (!inits$optim) {

    get_mle <- function(x) {
      x$cutoffs <- binormalLR_cutoffs(
        x, counts$events, counts$nonevents, method = control$init_cutoffs
      )
      binormalLR_mle(x, counts$events, counts$nonevents, control = control)
    }

    inits_list <- list(
      res,
      convert_inits(binormal_inits(
        truth, rating, counts$cutoffs, intercept = FALSE, shift = 0, scale = 2,
        prior_count = prior_count
      )),
      convert_inits(binormal_inits(
        truth, rating, counts$cutoffs, intercept = FALSE, shift = 1, scale = 2,
        prior_count = prior_count
      ))
    )

    res <- lapply(inits_list, get_mle)

    d_a_range <- range(mapply(getElement, res, "d_a"))
    c_range <- range(mapply(getElement, res, "c"))

    if (diff(d_a_range) > 0.1 || diff(c_range) > 0.1) {

      if (prod(c_range) < 0 && all(d_a_range < 0.5 * abs(c_range))) {
        auc_range <- mapply(.auc.binormalLR, d_a_range, c_range)
        z <- qnorm(sum(auc_range) / 2)
        d_a1 <- d_a_range[1]
        d_a2 <- sqrt(2) * z - d_a_range[1]
        c1 <- c_range[2]
        c2 <- c_range[1] - c_range[1]
        d_a_phase <- 2
      } else if (d_a_range[2] > d_a_range[1]) {
        d_a1 <- d_a_range[1]
        d_a2 <- d_a_range[2] - d_a_range[1]
        c1 <- c_range[2]
        c2 <- c_range[1] - c_range[2]
        d_a_phase <- 1
      } else {
        d_a1 <- d_a_range[2]
        d_a2 <- d_a_range[1] - d_a_range[2]
        c1 <- c_range[1]
        c2 <- c_range[2] - c_range[1]
        d_a_phase <- 1
      }

      n <- 5
      theta <- 1:n * pi / (2 * (1 + n))
      inits2 <- Map(list,
        d_a = d_a1 + d_a2 * cos(d_a_phase * theta - (d_a_phase - 1) * pi / 2),
        c = c1 + c2 * sin(theta),
        cutoffs = list(inits_list[[1]]$cutoffs)
      )

      res <- c(res, lapply(inits2, get_mle))

    }

    values <- sapply(res, function(x) x$optim$value)
    res <- res[[which.min(values)]]
  } else {
    warning(
      "\nLikelihood ratio binormal curve fit has AUC = ", 
      0.5 * (inits$a %in% c(-Inf, 0)) + 1 * (inits$a == Inf),
      " due to lack of interior points.",
      if (inits$a != 0) "\nConsider fitting an empirical curve instead."
    )
  }

  res
}


binormalLR_mle <- function(inits, events, nonevents, control = list()) {
  control <- do.call(
    function(
      method = c("trust", "CG", "BFGS", "L-BFGS-B", "Nelder-Mead"),
      gradient = TRUE, ...
    ) {
      list(method = match.arg(method), gradient = gradient)
    },
    control
  )

  inits <- c(inits$d_a, inits$c, inits$cutoffs)
  switch(control$method,
    "trust" = {
      control$gradient <- TRUE
      res <- trust::trust(
        function(x) {
          res <- binormalLR_objfun(
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
        fn = function(x) -binormalLR_loglik(x, events, nonevents),
        gr = if (control$gradient) {
          function(x) -binormalLR_gradient(x, events, nonevents)
        },
        method = control$method,
        lower = c(0, -1, rep(-Inf, num_cutoffs)),
        upper = c(Inf, 1, rep(Inf, num_cutoffs))
      )
      estimates <- res$par
    },
    {
      res <- optim(
        inits,
        fn = function(x) -binormalLR_loglik(x, events, nonevents),
        gr = if (control$gradient) {
          function(x) -binormalLR_gradient(x, events, nonevents)
        },
        method = control$method
      )
      estimates <- res$par
    }
  )
  res$control <- control

  list(
    d_a = estimates[1],
    c = estimates[2],
    cutoffs = estimates[-(1:2)],
    inits = inits,
    optim = res
  )
}


binormalLR_cutoff <- function(fpr, d_a, c, init = 0) {
  bound <- d_a * sqrt(1 + c^2) / (4 * c)
  interval <- c(-Inf, Inf)
  interval[(c > 0) + 1] <- bound
  f <- function(x, target) abs(binormalLR_fpr(x, d_a, c) - target)
  sapply(fpr, function(target) {
    if (target <= 0) {
      interval[2]
    } else if (target >= 1) {
      interval[1]
    } else {
      optim(
        init,
        f,
        target = target,
        method = "L-BFGS-B",
        lower = interval[1],
        upper = interval[2]
      )$par
    }
  })
}


binormalLR_cutoffs <- function(
  inits, events, nonevents, method = c("uni", "multi")
) {
  switch(match.arg(method),
    "multi" = binormalLR_cutoffs_multi(inits, events, nonevents),
    "uni" = binormalLR_cutoffs_uni(inits, events, nonevents)
  )
}


binormalLR_cutoffs_multi <- function(inits, events, nonevents) {
  cutoffs <- inits$cutoffs
  bound <- inits$d_a * sqrt(1 + inits$c^2) / (4 * inits$c)
  cutoffs_range <- range(cutoffs)
  if (inits$c < 0 && cutoffs_range[1] < bound) {
    cutoffs <- cutoffs + (bound - cutoffs_range[1] + 1e-10)
  } else if (inits$c > 0 && cutoffs_range[2] > bound) {
    cutoffs <- cutoffs + (bound - cutoffs_range[2] - 1e-10)
  }
  fixed <- c(inits$d_a, inits$c)
  optim(
    cutoffs,
    fn = function(x) -binormalLR_loglik(x, events, nonevents, fixed = fixed),
    gr = function(x) -binormalLR_gradient(x, events, nonevents, fixed = fixed),
    method = "BFGS"
  )$par
}


binormalLR_cutoffs_uni <- function(inits, events, nonevents) {
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
    binormalLR_loglik(
      cutoff, counts2$events, counts2$nonevents, fixed = c(inits$d_a, inits$c)
    )
  }

  inds <- seq_len(length(events) - 1)
  if (inits$c < 0) {
    get_interval <- function(bound1, bound2) {
      c(lower = bound1 + 1e-3, upper = bound2)
    }
  } else {
    inds <- rev(inds)
    get_interval <- function(bound1, bound2) {
      c(lower = bound2, upper = bound1 - 1e-3)
    }
  }

  res <- numeric(length(inds))

  counts2 <- get_counts2(cum_events[inds[1]], cum_nonevents[inds[1]])
  old_bound <- 0
  old_ll <- loglik2(old_bound)
  scale <- sign(inits$c) / (1 + abs(inits$c))
  maxiter <- 1000
  for (iter in 1:maxiter) {
    bound1 <- old_bound + scale
    ll1 <- loglik2(bound1)
    if (ll1 < old_ll) {
      break
    }
    old_bound <- bound1
    old_ll <- ll1
  }
  if (iter == maxiter) bound1 <- inits$d_a * sqrt(1 + inits$c^2) / (4 * inits$c)

  finc1 <- function(x, scale, iter) x + scale * iter
  finc2 <- function(x, scale, iter) x + scale * iter^(1 + iter >= 10)
  for (ind in inds) {
    counts2 <- get_counts2(cum_events[ind], cum_nonevents[ind])
    bound2 <- get_next_bound(bound1, -scale, loglik2, finc1, finc2)
    bound1 <- get_next_cutoff(get_interval(bound1, bound2), loglik2)
    res[ind] <- bound1
  }

  res
}


binormalLR_objfun <- function(
  x, events, nonevents, fixed = NULL, all = FALSE, loglik = all, gradient = all,
  hessian = all, loglik_min = -1e100
) {
  res <- list()

  x <- c(fixed, x)
  d_a <- x[1]
  c <- x[2]
  cutoffs <- x[-(1:2)]

  bound <- d_a * sqrt(1 + c^2) / (4 * c)
  in_support <- d_a >= 0 && abs(c) < 1 && all(diff(cutoffs) > 0) &&
    ifelse(c < 0, min(cutoffs) >= bound, max(cutoffs) <= bound)

  if (!in_support) {
    n <- length(x) - length(fixed)
    if (loglik) res$loglik <- loglik_min
    if (gradient) res$gradient <- numeric(n)
    if (hessian) res$hessian <- matrix(0, n, n)
    return(res)
  }

  p <- binormal_probs(binormalLR_fpr, cutoffs, d_a, c)
  q <- binormal_probs(binormalLR_tpr, cutoffs, d_a, c)

  if (loglik) {
    res$loglik <- max(c(nonevents %*% log(p) + events %*% log(q)), loglik_min)
  }

  if (gradient || hessian) {
    keep <- seq(length(fixed) + 1L, length(x))

    d1_d_a <- sqrt(1 + c^2) / 2
    d1_c <- (d_a * c) / (2 * sqrt(1 + c^2))
    df1 <- dnorm(binormalLR_z1_fpr(cutoffs, d_a, c))
    dt1 <- dnorm(binormalLR_z1_tpr(cutoffs, d_a, c))

    if (is_approx(c, 0)) {
      d2_d_a <- 0
      d2_c <- 0
      df2 <- 0
      dt2 <- 0
    } else {
      d2_d_a <- d1_d_a / c
      d2_c <- -(d_a * sqrt(1 + c^2)) / (2 * c^2) + d_a / (2 * sqrt(1 + c^2))
      df2 <- dnorm(binormalLR_z2_fpr(cutoffs, d_a, c))
      dt2 <- dnorm(binormalLR_z2_tpr(cutoffs, d_a, c))
    }

    dfpr_dd_a <- df1 * -d1_d_a + df2 * d2_d_a
    dfpr_dc <- df1 * (cutoffs - d1_c) + df2 * (cutoffs + d2_c)
    dfpr_dcut <- -(1 - c) * (df1 + df2)

    dtpr_dd_a <- dt1 * d1_d_a + dt2 * d2_d_a
    dtpr_dc <- dt1 * (-cutoffs + d1_c) + dt2 * (-cutoffs + d2_c)
    dtpr_dcut <- -(1 + c) * (dt1 + dt2)
  }

  if (gradient) {
    f <- function(counts, probs, dxr_dd_a, dxr_dc, dxr_dcut) {
      n <- length(counts)
      a <- -counts[-n] / probs[-n] + counts[-1] / probs[-1]
      c(a %*% dxr_dd_a, a %*% dxr_dc, a * dxr_dcut)[keep]
    }
    res$gradient <- f(nonevents, p, dfpr_dd_a, dfpr_dc, dfpr_dcut) +
      f(events, q, dtpr_dd_a, dtpr_dc, dtpr_dcut)
  }

  if (hessian) {
    f <- function(counts, probs, dxr_dd_a, dxr_dc, dxr_dcut) {
      num_cutoffs <- length(dxr_dcut)
      A <- matrix(0, num_cutoffs + 1, num_cutoffs + 2)
      A[, 1] <- c(0, dxr_dd_a) - c(dxr_dd_a, 0)
      A[, 2] <- c(0, dxr_dc) - c(dxr_dc, 0)
      inds1 <- seq_len(num_cutoffs)
      inds2 <- inds1 + 2
      A[cbind(inds1, inds2)] <- -dxr_dcut
      A[cbind(inds1 + 1, inds2)] <- dxr_dcut
      crossprod((-sum(counts) / probs) * A, A)[keep, keep]
    }
    res$hessian <- f(nonevents, p, dfpr_dd_a, dfpr_dc, dfpr_dcut) +
      f(events, q, dtpr_dd_a, dtpr_dc, dtpr_dcut)
  }

  res
}


binormalLR_loglik <- function(x, events, nonevents, fixed = NULL) {
  binormalLR_objfun(x, events, nonevents, fixed = fixed, loglik = TRUE)$loglik
}


binormalLR_gradient <- function(x, events, nonevents, fixed = NULL) {
  binormalLR_objfun(
    x, events, nonevents, fixed = fixed, gradient = TRUE
  )$gradient
}


binormalLR_fpr <- function(cutoff, d_a, c) {
  binormalLR_rate(
    binormalLR_z1_fpr(cutoff, d_a, c),
    binormalLR_z2_fpr(cutoff, d_a, c),
    c
  )
}


binormalLR_tpr <- function(cutoff, d_a, c) {
  binormalLR_rate(
    binormalLR_z1_tpr(cutoff, d_a, c),
    binormalLR_z2_tpr(cutoff, d_a, c),
    c
  )
}


binormalLR_rate <- function(z1, z2, c) {
  res <- pnorm(z1)
  if (is_approx(c, 0)) {
    res
  } else {
    res <- res + pnorm(z2)
    if (c < 0) res else res - 1
  }
}


binormalLR_z1_fpr <- function(cutoff, d_a, c) {
  -(1 - c) * cutoff - d_a * sqrt(1 + c^2) / 2
}


binormalLR_z2_fpr <- function(cutoff, d_a, c) {
  -(1 - c) * cutoff + d_a * sqrt(1 + c^2) / (2 * c)
}


binormalLR_z1_tpr <- function(cutoff, d_a, c) {
  -(1 + c) * cutoff + d_a * sqrt(1 + c^2) / 2
}


binormalLR_z2_tpr <- function(cutoff, d_a, c) {
  -(1 + c) * cutoff + d_a * sqrt(1 + c^2) / (2 * c)
}
