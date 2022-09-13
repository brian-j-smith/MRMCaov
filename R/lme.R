get_lme_args <- function(formula, object, data) {

  test_name <- object$vars["test"]
  reader_name <- object$vars["reader"]
  metric_name <- object$vars["metric"]

  lme_data <- object$data
  vars <- all.vars(formula[[3]])
  if (length(vars)) {
    vars_data_by <- by(data[vars], data[reader_name], function(split) {
      if (any(lengths(lapply(split, unique)) != 1)) {
        stop("righ-hand side formula values vary within readers")
      }
      split[1, , drop = FALSE]
    }, simplify = FALSE)
    vars_data <- cbind(dimnames(vars_data_by), do.call(rbind, vars_data_by))
    lme_data <- merge(lme_data, vars_data)
  }

  y <- lme_data[[metric_name]]
  fo <- update(formula[-2], paste("~", test_name, "+ ."))
  X <- model.matrix(fo, lme_data)
  fo <- as.formula(paste("~", reader_name, "- 1"))
  Z <- cbind(
    model.matrix(fo, lme_data),
    diag(length(y))
  )

  comps <- vcov_comps.mrmc(object)
  comps_cov <- comps$cov
  comps_var <- comps$var

  same_test <- outer(lme_data[[test_name]], lme_data[[test_name]], "==")
  same_reader <- outer(lme_data[[reader_name]], lme_data[[reader_name]], "==")

  R <- comps_cov[1] * (same_reader & !same_test) +
    comps_cov[2] * (same_test & !same_reader) +
    comps_cov[3] * (!same_reader & !same_test) +
    comps_var * (same_reader & same_test)

  list(
    data = lme_data,
    y = y,
    X = X,
    Z = Z,
    R = R,
    var = comps_var,
    cov = comps_cov,
    inits = comps$MS[2:3]
  )

}


get_lme_params <- function(inits, f, grad, y, X, Z, R, ...) {
  res <- optim(
    inits, fn = f, gr = grad,
    y = y, X = X, Z = Z, R = R,
    method = "L-BFGS-B",
    lower = 0
  )

  n <- c(ncol(Z) - nrow(Z), nrow(Z))
  V_hat <- get_lme_V(Z, get_lme_G(res$par, n), R)
  V_hat_inv <- chol2inv(chol(V_hat))

  beta_hat_cov <- chol2inv(chol(t(X) %*% V_hat_inv %*% X))
  beta_hat <- c(beta_hat_cov %*% t(X) %*% V_hat_inv %*% y)

  beta_names <- colnames(X)
  dimnames(beta_hat_cov) <- list(beta_names, beta_names)
  names(beta_hat) <- beta_names

  list(coef = beta_hat, cov = beta_hat_cov, optim = res)
}


## Negative log-likelihood + constant
negloglik_lme_reml <- function(params, y, X, Z, R, ...) {
  n <- c(ncol(Z) - nrow(Z), nrow(Z))

  V <- get_lme_V(Z, get_lme_G(params, n), R)
  V_chol <- chol(V)
  V_inv <- chol2inv(V_chol)

  W <- t(X) %*% V_inv %*% X
  W_chol <- chol(W)

  r <- y - X %*% chol2inv(W_chol) %*% t(X) %*% V_inv %*% y

  chol2det(V_chol, log = TRUE) +
    chol2det(W_chol, log = TRUE) +
    c(t(r) %*% V_inv %*% r)
}


## Gradient
grad_lme_reml <- function(params, y, X, Z, R, ...) {
  n <- c(ncol(Z) - nrow(Z), nrow(Z))

  V <- get_lme_V(Z, get_lme_G(params, n), R)
  V_inv <- chol2inv(chol(V))

  W_chol_inv <- backsolve(chol(t(X) %*% V_inv %*% X), diag(ncol(X)))
  r <- y - X %*% tcrossprod(W_chol_inv) %*% t(X) %*% V_inv %*% y
  X_star <- X %*% W_chol_inv

  Z_t_V_inv <- t(Z) %*% V_inv
  G1 <- get_lme_G(c(1, params[2]), n)
  G2 <- get_lme_G(c(params[1], 1), n)

  W <- Z_t_V_inv %*% Z
  g1 <- c(sum(diag(W %*% G1)), sum(diag(W %*% G2)))

  W <- Z_t_V_inv %*% r
  g2 <- -c(t(W) %*% G1 %*% W, t(W) %*% G2 %*% W)

  W <- Z_t_V_inv %*% X_star
  g3 <- -c(sum(diag(t(W) %*% G1 %*% W)), sum(diag(t(W) %*% G2 %*% W)))

  g1 + g2 + g3
}


## Random effects covariance
get_lme_G <- function(params, n) {
  diag(rep(params, times = n))
}


## Response covariance
get_lme_V <- function(Z, G, R) {
  Z %*% G %*% t(Z) + R
}
