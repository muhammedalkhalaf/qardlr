#' Rolling Window QARDL Estimation
#'
#' Performs rolling or recursive window QARDL estimation to assess
#' parameter stability over time.
#'
#' @param formula A formula of the form \code{y ~ x1 + x2 + ...}.
#' @param data A data frame containing the variables.
#' @param tau Numeric vector of quantiles. Default is \code{c(0.25, 0.50, 0.75)}.
#' @param p Integer. AR lag order. Default is 1.
#' @param q Integer. Distributed lag order. Default is 1.
#' @param window Integer. Rolling window size. If 0, uses 10\% of sample size.
#' @param method Character. Either \code{"rolling"} (fixed window) or
#'   \code{"recursive"} (expanding window). Default is \code{"rolling"}.
#' @param constant Logical. Include intercept. Default is \code{TRUE}.
#'
#' @return An object of class \code{"qardl_rolling"} containing:
#' \describe{
#'   \item{beta}{3D array of long-run parameters (k x ntau x nwindows)}
#'   \item{phi}{3D array of AR parameters (p x ntau x nwindows)}
#'   \item{gamma}{3D array of impact parameters (k x ntau x nwindows)}
#'   \item{rho}{Matrix of ECM coefficients (nwindows x ntau)}
#'   \item{wald_beta}{Matrix of beta constancy Wald statistics}
#'   \item{wald_phi}{Matrix of phi constancy Wald statistics}
#'   \item{wald_gamma}{Matrix of gamma constancy Wald statistics}
#'   \item{dates}{Vector of end dates for each window}
#'   \item{window}{Window size used}
#'   \item{method}{Method used ("rolling" or "recursive")}
#'   \item{tau}{Vector of quantiles}
#' }
#'
#' @details
#' Rolling window estimation helps detect structural breaks and assess
#' parameter stability. The function estimates QARDL models on successive
#' windows of data and tracks how parameters evolve over time.
#'
#' For \code{method = "rolling"}, a fixed window of size \code{window} is used.
#' For \code{method = "recursive"}, the window expands from \code{window}
#' to the full sample.
#'
#' @references
#' Cho, J.S., Kim, T.-H., & Shin, Y. (2015). Quantile cointegration in the
#' autoregressive distributed-lag modeling framework. \emph{Journal of
#' Econometrics}, 188(1), 281-300. \doi{10.1016/j.jeconom.2015.01.003}
#'
#' @seealso \code{\link{qardl}}, \code{\link{plot.qardl_rolling}}
#'
#' @examples
#' data(qardl_sim)
#'
#' # Rolling estimation with 50-observation window
#' roll <- qardl_rolling(y ~ x1 + x2, data = qardl_sim,
#'                       tau = c(0.25, 0.50, 0.75), p = 2, q = 2, window = 50)
#' print(roll)
#'
#' # Recursive estimation
#' recur <- qardl_rolling(y ~ x1 + x2, data = qardl_sim,
#'                        tau = c(0.50), p = 2, q = 2,
#'                        window = 50, method = "recursive")
#' print(recur)
#'
#' @export
qardl_rolling <- function(formula, data, tau = c(0.25, 0.50, 0.75),
                          p = 1L, q = 1L, window = 0L,
                          method = c("rolling", "recursive"),
                          constant = TRUE) {

  method <- match.arg(method)

  # Extract variables
  mf <- model.frame(formula, data = data, na.action = na.pass)
  y <- model.response(mf)
  X <- model.matrix(formula, data = mf)

  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }

  n <- length(y)
  k <- ncol(X)
  ntau <- length(tau)
  maxlag <- max(p, q)

  # Determine window size
  if (window == 0) {
    window <- max(as.integer(n * 0.1), maxlag + k + 10)
  }
  window <- as.integer(window)

  if (window < maxlag + k + 5) {
    window <- maxlag + k + 10
    message(sprintf("Window size increased to %d for stability", window))
  }

  if (window >= n) {
    stop("Window size must be smaller than sample size", call. = FALSE)
  }

  # Number of windows
  if (method == "rolling") {
    n_windows <- n - window + 1
  } else {
    n_windows <- n - window + 1
  }

  # Initialize storage
  beta_array <- array(NA, dim = c(k, ntau, n_windows))
  phi_array <- array(NA, dim = c(p, ntau, n_windows))
  gamma_array <- array(NA, dim = c(k, ntau, n_windows))
  rho_mat <- matrix(NA, nrow = n_windows, ncol = ntau)
  wald_beta <- matrix(NA, nrow = n_windows, ncol = 1)
  wald_phi <- matrix(NA, nrow = n_windows, ncol = 1)
  wald_gamma <- matrix(NA, nrow = n_windows, ncol = 1)
  end_dates <- integer(n_windows)

  message(sprintf("Running %s QARDL with window = %d, %d iterations",
                  method, window, n_windows))

  # Rolling estimation
  for (w in seq_len(n_windows)) {
    if (method == "rolling") {
      start_idx <- w
      end_idx <- w + window - 1
    } else {
      start_idx <- 1
      end_idx <- window + w - 1
    }

    y_sub <- y[start_idx:end_idx]
    X_sub <- X[start_idx:end_idx, , drop = FALSE]

    end_dates[w] <- end_idx

    # Estimate QARDL on subsample
    est <- tryCatch({
      qardl_estimate(y_sub, X_sub, p = p, q = q, tau = tau, constant = constant)
    }, error = function(e) NULL)

    if (is.null(est)) {
      next
    }

    # Compute long-run parameters
    lr <- tryCatch({
      compute_longrun(est, k = k, tau = tau)
    }, error = function(e) NULL)

    if (is.null(lr)) {
      next
    }

    # Store results
    beta_array[, , w] <- lr$beta
    phi_array[, , w] <- est$phi
    gamma_array[, , w] <- est$gamma
    rho_mat[w, ] <- lr$rho

    # Compute Wald statistics if ntau >= 2
    if (ntau >= 2) {
      wald_b <- tryCatch({
        wald_constancy_test(lr$beta, lr$beta_cov, est$nobs, "beta")$statistic
      }, error = function(e) NA)

      wald_p <- tryCatch({
        wald_constancy_test(est$phi, est$phi_cov, est$nobs, "phi")$statistic
      }, error = function(e) NA)

      wald_g <- tryCatch({
        wald_constancy_test(est$gamma, est$gamma_cov, est$nobs, "gamma")$statistic
      }, error = function(e) NA)

      wald_beta[w, 1] <- wald_b
      wald_phi[w, 1] <- wald_p
      wald_gamma[w, 1] <- wald_g
    }

    # Progress indicator
    if (w %% 50 == 0) {
      message(sprintf("  Completed %d/%d windows", w, n_windows))
    }
  }

  # Set dimension names
  dimnames(beta_array) <- list(colnames(X), paste0("tau=", tau), NULL)
  dimnames(phi_array) <- list(paste0("phi_", seq_len(p)), paste0("tau=", tau), NULL)
  dimnames(gamma_array) <- list(colnames(X), paste0("tau=", tau), NULL)
  colnames(rho_mat) <- paste0("tau=", tau)

  result <- list(
    beta = beta_array,
    phi = phi_array,
    gamma = gamma_array,
    rho = rho_mat,
    wald_beta = wald_beta,
    wald_phi = wald_phi,
    wald_gamma = wald_gamma,
    end_dates = end_dates,
    window = window,
    method = method,
    tau = tau,
    p = p,
    q = q,
    k = k,
    n_windows = n_windows,
    call = match.call()
  )

  class(result) <- "qardl_rolling"
  return(result)
}


#' Print Rolling QARDL Results
#'
#' @param x Object of class \code{"qardl_rolling"}.
#' @param ... Additional arguments (unused).
#'
#' @return Invisible \code{x}.
#'
#' @export
print.qardl_rolling <- function(x, ...) {

  cat("\n")
  cat(strrep("=", 70), "\n")
  cat(sprintf("  %s QARDL Estimation Results\n",
              ifelse(x$method == "rolling", "Rolling", "Recursive")))
  cat(strrep("=", 70), "\n\n")

  cat(sprintf("  Window size:  %d\n", x$window))
  cat(sprintf("  N windows:    %d\n", x$n_windows))
  cat(sprintf("  QARDL(%d, %d)\n", x$p, x$q))
  cat(sprintf("  Quantiles:    %s\n", paste(x$tau, collapse = ", ")))
  cat("\n")

  # Summary statistics for beta
  cat(strrep("-", 70), "\n")
  cat("  Long-Run Parameter (beta) Summary:\n")
  cat(strrep("-", 70), "\n")

  k <- x$k
  ntau <- length(x$tau)

  for (kk in seq_len(k)) {
    var_name <- dimnames(x$beta)[[1]][kk]
    cat(sprintf("\n  Variable: %s\n", var_name))
    cat(sprintf("  %10s %12s %12s %12s %12s\n",
                "Quantile", "Mean", "SD", "Min", "Max"))
    cat(strrep("-", 60), "\n")

    for (t_idx in seq_len(ntau)) {
      vals <- x$beta[kk, t_idx, ]
      vals <- vals[!is.na(vals)]
      if (length(vals) > 0) {
        cat(sprintf("  %10.2f %12.4f %12.4f %12.4f %12.4f\n",
                    x$tau[t_idx], mean(vals), sd(vals), min(vals), max(vals)))
      }
    }
  }

  cat("\n")
  cat(strrep("=", 70), "\n")

  invisible(x)
}


#' Plot Rolling QARDL Results
#'
#' Creates time series plots of rolling QARDL parameter estimates.
#'
#' @param x Object of class \code{"qardl_rolling"}.
#' @param which Character. Which parameter to plot: \code{"beta"}, \code{"phi"},
#'   \code{"gamma"}, or \code{"rho"}. Default is \code{"beta"}.
#' @param var Integer or character. Which variable to plot (for beta/gamma).
#'   Default is 1.
#' @param tau_idx Integer. Which quantile index to plot. Default is all.
#' @param ... Additional arguments passed to \code{plot}.
#'
#' @return Invisible \code{NULL}.
#'
#' @export
plot.qardl_rolling <- function(x, which = c("beta", "phi", "gamma", "rho"),
                               var = 1, tau_idx = NULL, ...) {

  which <- match.arg(which)
  n_windows <- x$n_windows
  ntau <- length(x$tau)

  if (is.null(tau_idx)) {
    tau_idx <- seq_len(ntau)
  }

  # Get data to plot
  if (which == "beta") {
    if (is.character(var)) {
      var_idx <- which(dimnames(x$beta)[[1]] == var)
    } else {
      var_idx <- var
    }
    plot_data <- x$beta[var_idx, tau_idx, , drop = FALSE]
    main_title <- sprintf("Rolling Beta: %s", dimnames(x$beta)[[1]][var_idx])
    ylab <- "beta"
  } else if (which == "phi") {
    var_idx <- var
    plot_data <- x$phi[var_idx, tau_idx, , drop = FALSE]
    main_title <- sprintf("Rolling Phi: lag %d", var_idx)
    ylab <- "phi"
  } else if (which == "gamma") {
    if (is.character(var)) {
      var_idx <- which(dimnames(x$gamma)[[1]] == var)
    } else {
      var_idx <- var
    }
    plot_data <- x$gamma[var_idx, tau_idx, , drop = FALSE]
    main_title <- sprintf("Rolling Gamma: %s", dimnames(x$gamma)[[1]][var_idx])
    ylab <- "gamma"
  } else {
    plot_data <- t(x$rho[, tau_idx, drop = FALSE])
    dim(plot_data) <- c(length(tau_idx), n_windows)
    main_title <- "Rolling Rho (ECM coefficient)"
    ylab <- "rho"
  }

  # Determine y-axis limits
  ylim <- range(plot_data, na.rm = TRUE)
  ylim <- ylim + c(-0.1, 0.1) * diff(ylim)

  # Colors for different quantiles
  colors <- grDevices::rainbow(length(tau_idx))

  # Create plot
  plot(x$end_dates, plot_data[1, ], type = "l", col = colors[1],
       ylim = ylim, xlab = "End Date Index", ylab = ylab,
       main = main_title, ...)

  if (length(tau_idx) > 1) {
    for (i in 2:length(tau_idx)) {
      lines(x$end_dates, plot_data[i, ], col = colors[i])
    }

    legend("topright", legend = paste0("tau=", x$tau[tau_idx]),
           col = colors, lty = 1, bty = "n", cex = 0.8)
  }

  invisible(NULL)
}
