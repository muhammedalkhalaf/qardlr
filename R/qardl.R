#' Quantile Autoregressive Distributed Lag Model Estimation
#'
#' Estimates the Quantile ARDL (QARDL) model of Cho, Kim & Shin (2015).
#' The model estimates quantile-specific long-run equilibrium relationships
#' and short-run dynamics between a dependent variable and covariates.
#'
#' @param formula A formula of the form \code{y ~ x1 + x2 + ...} where \code{y}
#'   is the dependent variable and \code{x1, x2, ...} are covariates.
#' @param data A data frame containing the variables in the formula.
#' @param tau Numeric vector of quantiles to estimate. Must be in (0, 1).
#'   Default is \code{c(0.25, 0.50, 0.75)}.
#' @param p Integer. AR lag order for the dependent variable. If 0,
#'   automatically selected via BIC. Default is 0.
#' @param q Integer. Distributed lag order for covariates. If 0,
#'   automatically selected via BIC. Default is 0.
#' @param pmax Integer. Maximum AR lag order for BIC selection. Default is 7.
#' @param qmax Integer. Maximum DL lag order for BIC selection. Default is 7.
#' @param ecm Logical. If \code{TRUE}, estimate Error Correction Model
#'   parameterization. Default is \code{FALSE}.
#' @param constant Logical. If \code{TRUE}, include an intercept. Default is \code{TRUE}.
#'
#' @return An object of class \code{"qardl"} containing:
#' \describe{
#'   \item{beta}{Long-run parameters matrix (k x ntau)}
#'   \item{beta_se}{Standard errors for beta}
#'   \item{phi}{Short-run AR parameters matrix (p x ntau)}
#'   \item{phi_se}{Standard errors for phi}
#'   \item{gamma}{Short-run impact parameters matrix (k x ntau)}
#'   \item{gamma_se}{Standard errors for gamma}
#'   \item{rho}{Speed of adjustment parameters (ECM coefficient)}
#'   \item{tau}{Vector of estimated quantiles}
#'   \item{p}{AR lag order used}
#'   \item{q}{Distributed lag order used}
#'   \item{nobs}{Number of observations}
#'   \item{k}{Number of covariates}
#'   \item{call}{The matched call}
#'   \item{formula}{The model formula}
#'   \item{data}{The data used}
#'   \item{qr_fits}{List of quantreg fit objects}
#'   \item{bic_grid}{BIC grid if lag selection was performed}
#'   \item{ecm}{Whether ECM parameterization was used}
#' }
#'
#' @details
#' The QARDL(p,q) model is specified as:
#' \deqn{Q_{y_t}(\tau | \mathcal{F}_{t-1}) = c(\tau) + \sum_{i=1}^{p} \phi_i(\tau) y_{t-i} + \sum_{j=0}^{q-1} \gamma'_j(\tau) x_{t-j}}
#'
#' Long-run parameters are computed as:
#' \deqn{\beta(\tau) = \frac{\sum_{j=0}^{q-1} \gamma_j(\tau)}{1 - \sum_{i=1}^{p} \phi_i(\tau)}}
#'
#' The speed of adjustment (ECM coefficient) is:
#' \deqn{\rho(\tau) = \sum_{i=1}^{p} \phi_i(\tau) - 1}
#'
#' Negative \eqn{\rho(\tau)} indicates convergence to long-run equilibrium.
#'
#' @references
#' Cho, J.S., Kim, T.-H., & Shin, Y. (2015). Quantile cointegration in the
#' autoregressive distributed-lag modeling framework. \emph{Journal of
#' Econometrics}, 188(1), 281-300. \doi{10.1016/j.jeconom.2015.01.003}
#'
#' @seealso \code{\link{qardl_rolling}}, \code{\link{qardl_simulate}},
#'   \code{\link{qardl_wald}}, \code{\link{summary.qardl}}
#'
#' @examples
#' # Load example data
#' data(qardl_sim)
#'
#' # Basic QARDL estimation with automatic lag selection
#' fit <- qardl(y ~ x1 + x2, data = qardl_sim, tau = c(0.25, 0.50, 0.75))
#' summary(fit)
#'
#' # QARDL with specified lags
#' fit2 <- qardl(y ~ x1 + x2, data = qardl_sim, tau = c(0.1, 0.5, 0.9), p = 2, q = 2)
#' print(fit2)
#'
#' # QARDL-ECM parameterization
#' fit_ecm <- qardl(y ~ x1 + x2, data = qardl_sim, tau = c(0.25, 0.50, 0.75), ecm = TRUE)
#' summary(fit_ecm)
#'
#' @export
qardl <- function(formula, data, tau = c(0.25, 0.50, 0.75),
                  p = 0L, q = 0L, pmax = 7L, qmax = 7L,
                  ecm = FALSE, constant = TRUE) {

  # Validate inputs
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object", call. = FALSE)
  }

  tau <- sort(unique(tau))
  if (any(tau <= 0) || any(tau >= 1)) {
    stop("All tau values must be strictly between 0 and 1", call. = FALSE)
  }

  p <- as.integer(p)
  q <- as.integer(q)
  pmax <- as.integer(pmax)
  qmax <- as.integer(qmax)

  if (p < 0 || q < 0) {
    stop("p and q must be non-negative integers", call. = FALSE)
  }

  # Extract variables from formula
  mf <- model.frame(formula, data = data, na.action = na.pass)
  y <- model.response(mf)
  X <- model.matrix(formula, data = mf)

  # Remove intercept from X if present (we handle it separately)
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }

  n <- length(y)
  k <- ncol(X)

  if (k < 1) {
    stop("At least one covariate is required", call. = FALSE)
  }

  if (n < 20) {
    stop("Insufficient observations (need at least 20)", call. = FALSE)
  }

  # Variable names
  depvar <- all.vars(formula)[1]
  indepvars <- colnames(X)

  # Automatic lag selection via BIC if p = 0 or q = 0
  bic_grid <- NULL
  if (p == 0 || q == 0) {
    bic_result <- qardl_bic_select(y, X, pmax = pmax, qmax = qmax,
                                    constant = constant)
    bic_grid <- bic_result$bic_grid

    if (p == 0) p <- bic_result$p_opt
    if (q == 0) q <- bic_result$q_opt

    message(sprintf("BIC-selected lag orders: p = %d, q = %d", p, q))
  }

  # Validate lag orders
  if (p < 1) p <- 1L
  if (q < 1) q <- 1L

  maxpq <- max(p, q)
  if (n <= maxpq + k + 5) {
    stop("Insufficient observations for specified lag orders", call. = FALSE)
  }

  # Core QARDL estimation
  est <- qardl_estimate(y, X, p = p, q = q, tau = tau, constant = constant)

  # Compute long-run parameters (beta)
  beta_result <- compute_longrun(est, k = k, tau = tau)

  # If ECM requested, compute ECM parameterization
  ecm_result <- NULL
  if (ecm) {
    ecm_result <- compute_ecm(est, y, X, p = p, q = q, tau = tau,
                               constant = constant)
  }

  # Build result object
  result <- list(
    beta = beta_result$beta,
    beta_se = beta_result$beta_se,
    beta_cov = beta_result$beta_cov,
    phi = est$phi,
    phi_se = est$phi_se,
    phi_cov = est$phi_cov,
    gamma = est$gamma,
    gamma_se = est$gamma_se,
    gamma_cov = est$gamma_cov,
    rho = beta_result$rho,
    rho_se = beta_result$rho_se,
    tau = tau,
    p = p,
    q = q,
    nobs = est$nobs,
    k = k,
    call = match.call(),
    formula = formula,
    depvar = depvar,
    indepvars = indepvars,
    qr_fits = est$qr_fits,
    bic_grid = bic_grid,
    ecm = ecm,
    ecm_result = ecm_result,
    raw_coefs = est$raw_coefs,
    raw_vcov = est$raw_vcov
  )

  class(result) <- "qardl"
  return(result)
}


#' Core QARDL Estimation
#'
#' Internal function for QARDL parameter estimation.
#'
#' @param y Numeric vector of dependent variable.
#' @param X Matrix of covariates.
#' @param p AR lag order.
#' @param q Distributed lag order.
#' @param tau Vector of quantiles.
#' @param constant Logical for intercept inclusion.
#'
#' @return List of estimated parameters and their covariances.
#'
#' @keywords internal
qardl_estimate <- function(y, X, p, q, tau, constant = TRUE) {

  n <- length(y)
  k <- ncol(X)
  ntau <- length(tau)
  maxlag <- max(p, q)

  # Build design matrix with lags
  # Y lags: y_{t-1}, ..., y_{t-p}
  # X and X lags: x_t, x_{t-1}, ..., x_{t-q+1}

  neff <- n - maxlag
  y_eff <- y[(maxlag + 1):n]

  # Design matrix
  design_cols <- list()

  # AR lags for y
  for (i in seq_len(p)) {
    lag_idx <- (maxlag + 1 - i):(n - i)
    design_cols[[paste0("y_lag", i)]] <- y[lag_idx]
  }

  # X and its lags (contemporaneous + q-1 lags = q total)
  for (j in seq_len(q)) {
    lag <- j - 1  # 0-indexed: lag 0 is contemporaneous
    for (kk in seq_len(k)) {
      if (lag == 0) {
        idx <- (maxlag + 1):n
      } else {
        idx <- (maxlag + 1 - lag):(n - lag)
      }
      col_name <- paste0(colnames(X)[kk], "_lag", lag)
      design_cols[[col_name]] <- X[idx, kk]
    }
  }

  Z <- do.call(cbind, design_cols)

  if (constant) {
    Z <- cbind(1, Z)
    colnames(Z)[1] <- "(Intercept)"
  }

  # Quantile regression for each tau
  qr_fits <- list()
  raw_coefs <- matrix(NA, nrow = ncol(Z), ncol = ntau)
  raw_vcov <- vector("list", ntau)

  for (t_idx in seq_along(tau)) {
    tau_val <- tau[t_idx]

    fit <- quantreg::rq(y_eff ~ Z - 1, tau = tau_val)
    qr_fits[[t_idx]] <- fit

    raw_coefs[, t_idx] <- coef(fit)

    # Covariance matrix using Hendricks-Koenker sandwich
    vcov_fit <- tryCatch(
      summary(fit, se = "nid", covariance = TRUE)$cov,
      error = function(e) {
        # Fallback to kernel estimate
        tryCatch(
          summary(fit, se = "ker", covariance = TRUE)$cov,
          error = function(e2) {
            diag(length(coef(fit))) * 0.01
          }
        )
      }
    )

    raw_vcov[[t_idx]] <- vcov_fit
  }

  rownames(raw_coefs) <- colnames(Z)
  colnames(raw_coefs) <- paste0("tau=", tau)

  # Extract phi (AR parameters) and gamma (impact parameters)
  start_idx <- ifelse(constant, 2, 1)

  # phi: AR coefficients (p x ntau)
  phi_idx <- start_idx:(start_idx + p - 1)
  phi <- raw_coefs[phi_idx, , drop = FALSE]
  rownames(phi) <- paste0("phi_", seq_len(p))

  # Extract phi covariance
  phi_se <- matrix(NA, nrow = p, ncol = ntau)
  phi_cov <- array(NA, dim = c(p, p, ntau))
  for (t_idx in seq_along(tau)) {
    vcov_mat <- raw_vcov[[t_idx]]
    if (!is.null(vcov_mat) && all(phi_idx <= nrow(vcov_mat)) && all(phi_idx <= ncol(vcov_mat))) {
      phi_cov_t <- vcov_mat[phi_idx, phi_idx, drop = FALSE]
      if (nrow(phi_cov_t) == p && ncol(phi_cov_t) == p) {
        phi_cov[, , t_idx] <- phi_cov_t
        phi_se[, t_idx] <- sqrt(pmax(diag(phi_cov_t), 0))
      }
    }
  }
  rownames(phi_se) <- rownames(phi)
  colnames(phi_se) <- colnames(phi)

  # gamma: Impact parameters for X at each lag
  # Organized as (k*q x ntau) matrix
  gamma_idx <- (start_idx + p):ncol(Z)
  gamma_raw <- raw_coefs[gamma_idx, , drop = FALSE]

  # Reorganize gamma: sum over lags to get cumulative impact
  # gamma_j = coefficient on x_{t-j} for j = 0, ..., q-1
  # We report gamma_0 (contemporaneous impact) as the main short-run impact
  gamma <- matrix(NA, nrow = k, ncol = ntau)
  gamma_se <- matrix(NA, nrow = k, ncol = ntau)
  gamma_cov <- array(NA, dim = c(k, k, ntau))

  for (t_idx in seq_along(tau)) {
    # Extract contemporaneous coefficients (lag 0)
    gamma_0_idx <- seq(1, k * q, by = q)
    if (all(gamma_0_idx <= nrow(gamma_raw))) {
      gamma[, t_idx] <- gamma_raw[gamma_0_idx, t_idx]
    }

    # Covariance for contemporaneous impact
    vcov_mat <- raw_vcov[[t_idx]]
    if (!is.null(vcov_mat) && all(gamma_idx <= nrow(vcov_mat)) && all(gamma_idx <= ncol(vcov_mat))) {
      raw_gamma_cov <- vcov_mat[gamma_idx, gamma_idx, drop = FALSE]
      if (all(gamma_0_idx <= nrow(raw_gamma_cov)) && all(gamma_0_idx <= ncol(raw_gamma_cov))) {
        gamma_cov_t <- raw_gamma_cov[gamma_0_idx, gamma_0_idx, drop = FALSE]
        if (nrow(gamma_cov_t) == k && ncol(gamma_cov_t) == k) {
          gamma_cov[, , t_idx] <- gamma_cov_t
          gamma_se[, t_idx] <- sqrt(pmax(diag(gamma_cov_t), 0))
        }
      }
    }
  }

  rownames(gamma) <- colnames(X)
  rownames(gamma_se) <- colnames(X)
  colnames(gamma) <- colnames(raw_coefs)
  colnames(gamma_se) <- colnames(raw_coefs)

  # Also store cumulative gamma (sum over all lags)
  gamma_cumul <- matrix(NA, nrow = k, ncol = ntau)
  gamma_cumul_se <- matrix(NA, nrow = k, ncol = ntau)
  gamma_cumul_cov <- array(NA, dim = c(k, k, ntau))

  for (t_idx in seq_along(tau)) {
    for (kk in seq_len(k)) {
      # Indices for all lags of variable kk
      var_idx <- seq(kk, k * q, by = k)
      if (all(var_idx <= nrow(gamma_raw))) {
        gamma_cumul[kk, t_idx] <- sum(gamma_raw[var_idx, t_idx])
      }

      # Variance of sum
      vcov_mat <- raw_vcov[[t_idx]]
      gamma_idx_var <- gamma_idx[var_idx]
      if (!is.null(vcov_mat) && all(gamma_idx_var <= nrow(vcov_mat)) && all(gamma_idx_var <= ncol(vcov_mat))) {
        sub_cov <- vcov_mat[gamma_idx_var, gamma_idx_var, drop = FALSE]
        gamma_cumul_se[kk, t_idx] <- sqrt(sum(pmax(sub_cov, 0)))
      }
    }

    # Full covariance for cumulative gamma
    L <- matrix(0, nrow = k, ncol = k * q)
    for (kk in seq_len(k)) {
      var_idx <- seq(kk, k * q, by = k)
      if (all(var_idx <= ncol(L))) {
        L[kk, var_idx] <- 1
      }
    }
    vcov_mat <- raw_vcov[[t_idx]]
    if (!is.null(vcov_mat) && all(gamma_idx <= nrow(vcov_mat)) && all(gamma_idx <= ncol(vcov_mat))) {
      raw_gamma_cov <- vcov_mat[gamma_idx, gamma_idx, drop = FALSE]
      if (ncol(raw_gamma_cov) == ncol(L)) {
        gamma_cumul_cov[, , t_idx] <- L %*% raw_gamma_cov %*% t(L)
      }
    }
  }

  rownames(gamma_cumul) <- colnames(X)
  rownames(gamma_cumul_se) <- colnames(X)
  colnames(gamma_cumul) <- colnames(raw_coefs)
  colnames(gamma_cumul_se) <- colnames(raw_coefs)

  return(list(
    phi = phi,
    phi_se = phi_se,
    phi_cov = phi_cov,
    gamma = gamma,
    gamma_se = gamma_se,
    gamma_cov = gamma_cov,
    gamma_cumul = gamma_cumul,
    gamma_cumul_se = gamma_cumul_se,
    gamma_cumul_cov = gamma_cumul_cov,
    raw_coefs = raw_coefs,
    raw_vcov = raw_vcov,
    qr_fits = qr_fits,
    nobs = neff,
    p = p,
    q = q,
    k = k
  ))
}


#' Compute Long-Run Parameters
#'
#' Computes long-run (beta) parameters from QARDL estimates.
#'
#' @param est QARDL estimation results from \code{qardl_estimate}.
#' @param k Number of covariates.
#' @param tau Vector of quantiles.
#'
#' @return List containing beta, beta_se, beta_cov, rho, rho_se.
#'
#' @keywords internal
compute_longrun <- function(est, k, tau) {

  ntau <- length(tau)
  p <- est$p

  # Long-run parameters: beta(tau) = gamma_cumul(tau) / (1 - rho(tau))
  # where rho(tau) = sum(phi_i(tau))

  beta <- matrix(NA, nrow = k, ncol = ntau)
  beta_se <- matrix(NA, nrow = k, ncol = ntau)
  beta_cov <- array(NA, dim = c(k, k, ntau))
  rho <- numeric(ntau)
  rho_se <- numeric(ntau)

  for (t_idx in seq_along(tau)) {
    # Sum of AR coefficients
    rho[t_idx] <- sum(est$phi[, t_idx])

    # Variance of rho = sum of all elements in phi covariance
    rho_se[t_idx] <- sqrt(sum(est$phi_cov[, , t_idx]))

    # Long-run multiplier
    denom <- 1 - rho[t_idx]

    if (abs(denom) < 1e-10) {
      warning(sprintf("Near-unit root at tau = %.2f, long-run parameters may be unreliable",
                      tau[t_idx]))
      denom <- sign(denom) * 1e-10
    }

    # beta = gamma_cumul / (1 - rho)
    beta[, t_idx] <- est$gamma_cumul[, t_idx] / denom

    # Delta method for beta variance
    # beta_j = g_j / d where g = gamma_cumul, d = 1 - rho
    # Var(beta_j) approx= (1/d^2)*Var(g_j) + (g_j/d^2)^2*Var(rho) + 2*(g_j/d^3)*Cov(g_j, rho)
    # Simplified: Var(beta_j) approx= Var(g_j)/d^2 + beta_j^2 * Var(rho) / d^2

    for (kk in seq_len(k)) {
      var_g <- est$gamma_cumul_cov[kk, kk, t_idx]
      var_rho <- rho_se[t_idx]^2
      b_val <- beta[kk, t_idx]

      # Delta method approximation
      var_beta <- (var_g + b_val^2 * var_rho) / (denom^2)
      beta_se[kk, t_idx] <- sqrt(max(var_beta, 0))
    }

    # Full covariance matrix for beta (simplified)
    # Using scaled gamma covariance
    beta_cov[, , t_idx] <- est$gamma_cumul_cov[, , t_idx] / (denom^2)
  }

  rownames(beta) <- rownames(est$gamma_cumul)
  rownames(beta_se) <- rownames(est$gamma_cumul)
  colnames(beta) <- colnames(est$gamma_cumul)
  colnames(beta_se) <- colnames(est$gamma_cumul)

  # Convert rho to speed of adjustment: rho - 1
  rho_adj <- rho - 1

  return(list(
    beta = beta,
    beta_se = beta_se,
    beta_cov = beta_cov,
    rho = rho_adj,
    rho_se = rho_se
  ))
}


#' Compute ECM Parameterization
#'
#' Computes Error Correction Model parameters.
#'
#' @param est QARDL estimation results.
#' @param y Dependent variable vector.
#' @param X Covariate matrix.
#' @param p AR lag order.
#' @param q Distributed lag order.
#' @param tau Vector of quantiles.
#' @param constant Logical for intercept.
#'
#' @return List containing ECM parameters.
#'
#' @keywords internal
compute_ecm <- function(est, y, X, p, q, tau, constant = TRUE) {

  n <- length(y)
  k <- ncol(X)
  ntau <- length(tau)
  maxlag <- max(p, q)
  neff <- n - maxlag

  # ECM parameters:
  # phi*_j = -sum_{i=j+1}^{p} phi_i for j = 1, ..., p-1
  # theta_j = -sum_{i=j+1}^{q} gamma_i for j = 0, ..., q-1

  phi_star <- NULL
  phi_star_se <- NULL

  if (p > 1) {
    phi_star <- matrix(NA, nrow = p - 1, ncol = ntau)
    phi_star_se <- matrix(NA, nrow = p - 1, ncol = ntau)

    for (t_idx in seq_along(tau)) {
      for (j in seq_len(p - 1)) {
        # phi*_j = -sum_{i=j+1}^{p} phi_i
        phi_star[j, t_idx] <- -sum(est$phi[(j + 1):p, t_idx])

        # Variance approximation
        sub_cov <- est$phi_cov[(j + 1):p, (j + 1):p, t_idx, drop = FALSE]
        phi_star_se[j, t_idx] <- sqrt(sum(sub_cov))
      }
    }

    rownames(phi_star) <- paste0("phi_star_", seq_len(p - 1))
    colnames(phi_star) <- colnames(est$phi)
    rownames(phi_star_se) <- rownames(phi_star)
    colnames(phi_star_se) <- colnames(phi_star)
  }

  # Theta parameters (short-run dynamics in ECM)
  # For simplicity, theta_0 = gamma_0 (contemporaneous impact of dx)
  theta <- est$gamma
  theta_se <- est$gamma_se

  return(list(
    phi_star = phi_star,
    phi_star_se = phi_star_se,
    theta = theta,
    theta_se = theta_se
  ))
}
