#' Monte Carlo Simulation for QARDL
#'
#' Performs Monte Carlo simulation to assess the finite-sample properties
#' of QARDL estimators under specified data generating processes.
#'
#' @param nobs Integer. Sample size for each simulation. Default is 200.
#' @param reps Integer. Number of Monte Carlo replications. Default is 1000.
#' @param tau Numeric vector of quantiles. Default is \code{c(0.25, 0.50, 0.75)}.
#' @param p Integer. AR lag order. Default is 1.
#' @param q Integer. Distributed lag order. Default is 1.
#' @param k Integer. Number of covariates. Default is 1.
#' @param beta_true Numeric vector. True long-run parameters (length k).
#'   Default is \code{rep(1, k)}.
#' @param phi_true Numeric vector. True AR parameters (length p).
#'   Default is \code{rep(0.5, p)}.
#' @param gamma_true Numeric vector. True impact parameters (length k).
#'   Default is \code{rep(0.3, k)}.
#' @param sigma_u Numeric. Standard deviation of the error term. Default is 1.
#' @param sigma_x Numeric. Standard deviation of covariate innovations. Default is 1.
#' @param seed Integer. Random seed for reproducibility. Default is \code{NULL}.
#' @param parallel Logical. Use parallel processing. Default is \code{FALSE}.
#' @param ncores Integer. Number of cores for parallel processing.
#'   Default is \code{parallel::detectCores() - 1}.
#'
#' @return An object of class \code{"qardl_mc"} containing:
#' \describe{
#'   \item{beta_sim}{Array of simulated beta estimates (k x ntau x reps)}
#'   \item{phi_sim}{Array of simulated phi estimates (p x ntau x reps)}
#'   \item{gamma_sim}{Array of simulated gamma estimates (k x ntau x reps)}
#'   \item{beta_true}{True beta values}
#'   \item{phi_true}{True phi values}
#'   \item{gamma_true}{True gamma values}
#'   \item{bias_beta}{Bias in beta estimates}
#'   \item{rmse_beta}{RMSE of beta estimates}
#'   \item{coverage_beta}{Empirical coverage of 95\% CI for beta}
#'   \item{reps}{Number of replications}
#'   \item{nobs}{Sample size}
#'   \item{tau}{Vector of quantiles}
#' }
#'
#' @details
#' The data generating process is:
#' \deqn{y_t = \sum_{i=1}^{p} \phi_i y_{t-i} + \sum_{j=1}^{k} \gamma_j x_{jt} + u_t}
#'
#' where \eqn{u_t \sim N(0, \sigma_u^2)} and \eqn{x_{jt}} follows a random walk
#' with innovations \eqn{\sim N(0, \sigma_x^2)}.
#'
#' @references
#' Cho, J.S., Kim, T.-H., & Shin, Y. (2015). Quantile cointegration in the
#' autoregressive distributed-lag modeling framework. \emph{Journal of
#' Econometrics}, 188(1), 281-300. \doi{10.1016/j.jeconom.2015.01.003}
#'
#' @seealso \code{\link{qardl}}, \code{\link{print.qardl_mc}}
#'
#' @examples
#' # Small simulation for illustration
#' mc <- qardl_simulate(nobs = 100, reps = 50, tau = c(0.25, 0.50, 0.75),
#'                      p = 1, q = 1, k = 1, seed = 123)
#' print(mc)
#'
#' @export
qardl_simulate <- function(nobs = 200L, reps = 1000L,
                           tau = c(0.25, 0.50, 0.75),
                           p = 1L, q = 1L, k = 1L,
                           beta_true = NULL, phi_true = NULL, gamma_true = NULL,
                           sigma_u = 1, sigma_x = 1,
                           seed = NULL, parallel = FALSE, ncores = NULL) {

  nobs <- as.integer(nobs)
  reps <- as.integer(reps)
  p <- as.integer(p)
  q <- as.integer(q)
  k <- as.integer(k)
  ntau <- length(tau)

  # Set defaults for true parameters
  if (is.null(beta_true)) beta_true <- rep(1, k)
  if (is.null(phi_true)) phi_true <- rep(0.5, p)
  if (is.null(gamma_true)) gamma_true <- rep(0.3, k)

  if (length(beta_true) != k) stop("beta_true must have length k")
  if (length(phi_true) != p) stop("phi_true must have length p")
  if (length(gamma_true) != k) stop("gamma_true must have length k")

  # Check stationarity
  rho_sum <- sum(phi_true)
  if (abs(rho_sum) >= 1) {
    warning("DGP may be non-stationary (sum of phi >= 1)")
  }

  # Set seed
  if (!is.null(seed)) set.seed(seed)

  # Initialize storage
  beta_sim <- array(NA, dim = c(k, ntau, reps))
  phi_sim <- array(NA, dim = c(p, ntau, reps))
  gamma_sim <- array(NA, dim = c(k, ntau, reps))
  se_beta <- array(NA, dim = c(k, ntau, reps))
  se_phi <- array(NA, dim = c(p, ntau, reps))
  se_gamma <- array(NA, dim = c(k, ntau, reps))

  message(sprintf("Running %d Monte Carlo replications (n = %d)...", reps, nobs))

  # Main simulation loop
  for (r in seq_len(reps)) {

    # Generate data
    sim_data <- dgp_qardl(nobs = nobs, p = p, q = q, k = k,
                          phi = phi_true, gamma = gamma_true,
                          sigma_u = sigma_u, sigma_x = sigma_x)

    y <- sim_data$y
    X <- sim_data$X

    # Estimate QARDL
    est <- tryCatch({
      qardl_estimate(y, X, p = p, q = q, tau = tau, constant = TRUE)
    }, error = function(e) NULL)

    if (is.null(est)) next

    # Compute long-run parameters
    lr <- tryCatch({
      compute_longrun(est, k = k, tau = tau)
    }, error = function(e) NULL)

    if (is.null(lr)) next

    # Store estimates
    beta_sim[, , r] <- lr$beta
    phi_sim[, , r] <- est$phi
    gamma_sim[, , r] <- est$gamma
    se_beta[, , r] <- lr$beta_se
    se_phi[, , r] <- est$phi_se
    se_gamma[, , r] <- est$gamma_se

    # Progress indicator
    if (r %% 100 == 0) {
      message(sprintf("  Completed %d/%d replications", r, reps))
    }
  }

  # Compute summary statistics
  # Bias
  bias_beta <- apply(beta_sim, c(1, 2), mean, na.rm = TRUE) -
    matrix(beta_true, nrow = k, ncol = ntau)
  bias_phi <- apply(phi_sim, c(1, 2), mean, na.rm = TRUE) -
    matrix(phi_true, nrow = p, ncol = ntau)
  bias_gamma <- apply(gamma_sim, c(1, 2), mean, na.rm = TRUE) -
    matrix(gamma_true, nrow = k, ncol = ntau)

  # RMSE
  rmse_beta <- matrix(NA, nrow = k, ncol = ntau)
  rmse_phi <- matrix(NA, nrow = p, ncol = ntau)
  rmse_gamma <- matrix(NA, nrow = k, ncol = ntau)

  for (kk in seq_len(k)) {
    for (t_idx in seq_len(ntau)) {
      vals <- beta_sim[kk, t_idx, ]
      rmse_beta[kk, t_idx] <- sqrt(mean((vals - beta_true[kk])^2, na.rm = TRUE))
    }
  }

  for (pp in seq_len(p)) {
    for (t_idx in seq_len(ntau)) {
      vals <- phi_sim[pp, t_idx, ]
      rmse_phi[pp, t_idx] <- sqrt(mean((vals - phi_true[pp])^2, na.rm = TRUE))
    }
  }

  for (kk in seq_len(k)) {
    for (t_idx in seq_len(ntau)) {
      vals <- gamma_sim[kk, t_idx, ]
      rmse_gamma[kk, t_idx] <- sqrt(mean((vals - gamma_true[kk])^2, na.rm = TRUE))
    }
  }

  # Coverage (95% CI)
  coverage_beta <- matrix(NA, nrow = k, ncol = ntau)
  for (kk in seq_len(k)) {
    for (t_idx in seq_len(ntau)) {
      est_vals <- beta_sim[kk, t_idx, ]
      se_vals <- se_beta[kk, t_idx, ]
      lower <- est_vals - 1.96 * se_vals
      upper <- est_vals + 1.96 * se_vals
      coverage_beta[kk, t_idx] <- mean(lower <= beta_true[kk] & beta_true[kk] <= upper,
                                        na.rm = TRUE)
    }
  }

  # Set names
  rownames(bias_beta) <- paste0("beta_", seq_len(k))
  rownames(rmse_beta) <- paste0("beta_", seq_len(k))
  rownames(coverage_beta) <- paste0("beta_", seq_len(k))
  colnames(bias_beta) <- paste0("tau=", tau)
  colnames(rmse_beta) <- paste0("tau=", tau)
  colnames(coverage_beta) <- paste0("tau=", tau)

  result <- list(
    beta_sim = beta_sim,
    phi_sim = phi_sim,
    gamma_sim = gamma_sim,
    beta_true = beta_true,
    phi_true = phi_true,
    gamma_true = gamma_true,
    bias_beta = bias_beta,
    bias_phi = bias_phi,
    bias_gamma = bias_gamma,
    rmse_beta = rmse_beta,
    rmse_phi = rmse_phi,
    rmse_gamma = rmse_gamma,
    coverage_beta = coverage_beta,
    reps = reps,
    nobs = nobs,
    tau = tau,
    p = p,
    q = q,
    k = k,
    call = match.call()
  )

  class(result) <- "qardl_mc"
  return(result)
}


#' Data Generating Process for QARDL
#'
#' Internal function to generate data from a QARDL DGP.
#'
#' @param nobs Sample size.
#' @param p AR lag order.
#' @param q Distributed lag order.
#' @param k Number of covariates.
#' @param phi AR parameters.
#' @param gamma Impact parameters.
#' @param sigma_u Error standard deviation.
#' @param sigma_x Covariate innovation standard deviation.
#'
#' @return List with y and X.
#'
#' @keywords internal
dgp_qardl <- function(nobs, p, q, k, phi, gamma, sigma_u, sigma_x) {

  # Burn-in period
  burnin <- 100
  n_total <- nobs + burnin

  # Generate covariates (random walks)
  X <- matrix(0, nrow = n_total, ncol = k)
  for (kk in seq_len(k)) {
    innov <- rnorm(n_total, sd = sigma_x)
    X[, kk] <- cumsum(innov)
  }

  # Generate dependent variable
  y <- numeric(n_total)
  u <- rnorm(n_total, sd = sigma_u)

  # Initialize
  for (t in 1:max(p, q)) {
    y[t] <- u[t]
  }

  # Generate y
  for (t in (max(p, q) + 1):n_total) {
    # AR component
    ar_component <- sum(phi * y[(t - 1):(t - p)])

    # X impact (contemporaneous)
    x_component <- sum(gamma * X[t, ])

    y[t] <- ar_component + x_component + u[t]
  }

  # Remove burn-in
  y <- y[(burnin + 1):n_total]
  X <- X[(burnin + 1):n_total, , drop = FALSE]
  colnames(X) <- paste0("x", seq_len(k))

  return(list(y = y, X = X))
}


#' Print Monte Carlo Results
#'
#' @param x Object of class \code{"qardl_mc"}.
#' @param digits Number of decimal places. Default is 4.
#' @param ... Additional arguments (unused).
#'
#' @return Invisible \code{x}.
#'
#' @export
print.qardl_mc <- function(x, digits = 4, ...) {

  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("  QARDL Monte Carlo Simulation Results\n")
  cat(strrep("=", 70), "\n\n")

  cat(sprintf("  Replications:  %d\n", x$reps))
  cat(sprintf("  Sample size:   %d\n", x$nobs))
  cat(sprintf("  QARDL(%d, %d)\n", x$p, x$q))
  cat(sprintf("  Covariates:    %d\n", x$k))
  cat(sprintf("  Quantiles:     %s\n", paste(x$tau, collapse = ", ")))
  cat("\n")

  # True parameters
  cat(strrep("-", 70), "\n")
  cat("  True Parameters:\n")
  cat(sprintf("    beta:  %s\n", paste(round(x$beta_true, digits), collapse = ", ")))
  cat(sprintf("    phi:   %s\n", paste(round(x$phi_true, digits), collapse = ", ")))
  cat(sprintf("    gamma: %s\n", paste(round(x$gamma_true, digits), collapse = ", ")))
  cat(strrep("-", 70), "\n\n")

  # Beta results
  cat("  Long-Run Parameter (beta) Results:\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("  %12s %12s %12s %12s %12s\n",
              "Quantile", "True", "Bias", "RMSE", "Coverage"))
  cat(strrep("-", 70), "\n")

  ntau <- length(x$tau)
  k <- x$k

  for (kk in seq_len(k)) {
    if (k > 1) cat(sprintf("  Variable %d:\n", kk))
    for (t_idx in seq_len(ntau)) {
      cat(sprintf("  %12.2f %12.4f %12.4f %12.4f %12.2f%%\n",
                  x$tau[t_idx], x$beta_true[kk],
                  x$bias_beta[kk, t_idx], x$rmse_beta[kk, t_idx],
                  x$coverage_beta[kk, t_idx] * 100))
    }
  }

  cat(strrep("-", 70), "\n\n")
  cat(strrep("=", 70), "\n")

  invisible(x)
}
