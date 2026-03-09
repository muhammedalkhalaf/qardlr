#' Print QARDL Results
#'
#' Print method for QARDL estimation results.
#'
#' @param x An object of class \code{"qardl"}.
#' @param digits Number of decimal places. Default is 4.
#' @param ... Additional arguments (unused).
#'
#' @return Invisible \code{x}.
#'
#' @export
print.qardl <- function(x, digits = 4, ...) {

  cat("\n")
  cat(strrep("=", 70), "\n")
  if (x$ecm) {
    cat("  QARDL-ECM Estimation Results\n")
  } else {
    cat("  QARDL Estimation Results\n")
  }
  cat(strrep("=", 70), "\n")
  cat("  Cho, Kim & Shin (2015), Journal of Econometrics\n")
  cat(strrep("=", 70), "\n\n")

  cat(sprintf("  Dep. variable  : %s\n", x$depvar))
  cat(sprintf("  Indep. vars    : %s\n", paste(x$indepvars, collapse = ", ")))
  cat(sprintf("  Observations   : %d\n", x$nobs))
  cat(sprintf("  QARDL(%d, %d)\n", x$p, x$q))
  cat(sprintf("  Quantiles      : %s\n", paste(x$tau, collapse = " ")))
  cat("\n")

  # Long-run parameters (beta)
  cat(strrep("-", 70), "\n")
  cat("  Long-Run Parameters: beta(tau)\n")
  cat("  beta_j(tau) = gamma_j(tau) / (1 - sum(phi_i(tau)))\n")
  cat(strrep("-", 70), "\n")

  print_param_table(x$beta, x$beta_se, x$tau, x$nobs, "Variable", digits)

  # Short-run AR parameters (phi)
  cat("\n")
  cat(strrep("-", 70), "\n")
  cat("  Short-Run AR Parameters: phi(tau)\n")
  cat(strrep("-", 70), "\n")

  print_param_table(x$phi, x$phi_se, x$tau, x$nobs, "Lag", digits)

  # Short-run impact parameters (gamma)
  cat("\n")
  cat(strrep("-", 70), "\n")
  cat("  Short-Run Impact Parameters: gamma(tau)\n")
  cat(strrep("-", 70), "\n")

  print_param_table(x$gamma, x$gamma_se, x$tau, x$nobs, "Variable", digits)

  cat("\n")
  cat(strrep("=", 70), "\n")

  invisible(x)
}


#' Summary of QARDL Results
#'
#' Provides a detailed summary of QARDL estimation results including
#' parameter estimates, standard errors, t-statistics, p-values, and
#' diagnostic tests.
#'
#' @param object An object of class \code{"qardl"}.
#' @param wald Logical. Include Wald tests for parameter constancy.
#'   Default is \code{TRUE}.
#' @param digits Number of decimal places. Default is 4.
#' @param ... Additional arguments (unused).
#'
#' @return An object of class \code{"summary.qardl"} (invisibly).
#'
#' @export
summary.qardl <- function(object, wald = TRUE, digits = 4, ...) {

  x <- object

  cat("\n")
  cat(strrep("=", 70), "\n")
  if (x$ecm) {
    cat("  QARDL-ECM Estimation Summary\n")
  } else {
    cat("  QARDL Estimation Summary\n")
  }
  cat(strrep("=", 70), "\n")
  cat("  Reference: Cho, Kim & Shin (2015), Journal of Econometrics\n")
  cat("  DOI: 10.1016/j.jeconom.2015.01.003\n")
  cat(strrep("=", 70), "\n\n")

  # Model info
  cat(sprintf("  Dependent variable  : %s\n", x$depvar))
  cat(sprintf("  Covariates          : %s\n", paste(x$indepvars, collapse = ", ")))
  cat(sprintf("  Observations        : %d\n", x$nobs))
  cat(sprintf("  QARDL(%d, %d)\n", x$p, x$q))
  cat(sprintf("  Quantiles           : %s\n", paste(x$tau, collapse = ", ")))
  cat("\n")

  # BIC grid if available
  if (!is.null(x$bic_grid)) {
    cat(strrep("-", 70), "\n")
    cat("  Lag selection performed via BIC (see print_bic_grid())\n")
    cat(strrep("-", 70), "\n\n")
  }

  ntau <- length(x$tau)
  k <- x$k
  p <- x$p

  # ============================================================
  # Long-run parameters (beta)
  # ============================================================
  cat(strrep("=", 70), "\n")
  cat("  LONG-RUN PARAMETERS: beta(tau)\n")
  cat("  beta_j(tau) = cumulative_gamma_j(tau) / (1 - rho(tau))\n")
  cat(strrep("=", 70), "\n\n")

  print_detailed_table(x$beta, x$beta_se, x$tau, x$nobs,
                       rownames(x$beta), "Variable", digits)

  # ============================================================
  # Short-run AR parameters (phi)
  # ============================================================
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("  SHORT-RUN AR PARAMETERS: phi(tau)\n")
  cat(strrep("=", 70), "\n\n")

  print_detailed_table(x$phi, x$phi_se, x$tau, x$nobs,
                       rownames(x$phi), "Lag", digits)

  # ============================================================
  # Short-run impact parameters (gamma)
  # ============================================================
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("  SHORT-RUN IMPACT PARAMETERS: gamma(tau)\n")
  cat(strrep("=", 70), "\n\n")

  print_detailed_table(x$gamma, x$gamma_se, x$tau, x$nobs,
                       rownames(x$gamma), "Variable", digits)

  # ============================================================
  # ECM Speed of Adjustment (rho)
  # ============================================================
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("  ECM SPEED OF ADJUSTMENT: rho(tau) = sum(phi_i(tau)) - 1\n")
  cat("  (rho < 0 implies convergence to long-run equilibrium)\n")
  cat(strrep("=", 70), "\n\n")

  cat(sprintf("  %10s %12s %12s %12s %12s %10s\n",
              "Quantile", "rho(tau)", "Std.Err.", "t-stat", "p-value", "Signal"))
  cat(strrep("-", 70), "\n")

  for (t_idx in seq_len(ntau)) {
    rho_val <- x$rho[t_idx]
    se_val <- x$rho_se[t_idx]
    tstat <- rho_val / se_val
    pval <- 2 * (1 - stats::pnorm(abs(tstat)))
    signal <- ifelse(rho_val < 0, "Converge", "Diverge")

    stars <- get_stars(pval)
    cat(sprintf("  %10.2f %12.4f %12.4f %12.3f %11.4f%s %10s\n",
                x$tau[t_idx], rho_val, se_val, tstat, pval, stars, signal))
  }
  cat(strrep("-", 70), "\n")

  # ============================================================
  # ECM parameters if requested
  # ============================================================
  if (x$ecm && !is.null(x$ecm_result)) {
    cat("\n")
    cat(strrep("=", 70), "\n")
    cat("  ECM PARAMETERIZATION\n")
    cat(strrep("=", 70), "\n\n")

    if (!is.null(x$ecm_result$phi_star)) {
      cat("  ECM Short-Run AR Parameters: phi*(tau)\n")
      cat(strrep("-", 70), "\n")
      print_detailed_table(x$ecm_result$phi_star, x$ecm_result$phi_star_se,
                           x$tau, x$nobs, rownames(x$ecm_result$phi_star),
                           "Lag", digits)
    }

    cat("\n  ECM Short-Run Impact: theta(tau)\n")
    cat(strrep("-", 70), "\n")
    print_detailed_table(x$ecm_result$theta, x$ecm_result$theta_se,
                         x$tau, x$nobs, rownames(x$ecm_result$theta),
                         "Variable", digits)
  }

  # ============================================================
  # Wald tests
  # ============================================================
  if (wald && ntau >= 2) {
    cat("\n")
    wald_result <- qardl_wald(x, type = "all", pairwise = FALSE)
    print(wald_result)
  }

  cat("\n")
  cat("  Significance codes: *** p<0.01, ** p<0.05, * p<0.10\n")
  cat(strrep("=", 70), "\n\n")

  invisible(x)
}


#' Print Parameter Table
#'
#' Internal function to print formatted parameter table.
#'
#' @keywords internal
print_param_table <- function(est, se, tau, nobs, row_label, digits = 4) {

  nrows <- nrow(est)
  ntau <- length(tau)

  cat(sprintf("  %12s", row_label))
  for (t_idx in seq_len(ntau)) {
    cat(sprintf(" %12s", paste0("tau=", tau[t_idx])))
  }
  cat("\n")
  cat(strrep("-", 70), "\n")

  for (i in seq_len(nrows)) {
    cat(sprintf("  %12s", rownames(est)[i]))
    for (t_idx in seq_len(ntau)) {
      cat(sprintf(" %12.4f", est[i, t_idx]))
    }
    cat("\n")

    # Standard errors in parentheses
    cat(sprintf("  %12s", ""))
    for (t_idx in seq_len(ntau)) {
      cat(sprintf(" (%10.4f)", se[i, t_idx]))
    }
    cat("\n")
  }
}


#' Print Detailed Parameter Table
#'
#' Internal function to print detailed parameter table with t-stats and p-values.
#'
#' @keywords internal
print_detailed_table <- function(est, se, tau, nobs, row_names, row_label, digits = 4) {

  nrows <- nrow(est)
  ntau <- length(tau)

  for (t_idx in seq_len(ntau)) {
    cat(sprintf("  ---- tau = %.2f %s\n", tau[t_idx], strrep("-", 52)))
    cat(sprintf("  %12s %12s %10s %10s %10s\n",
                row_label, "Estimate", "Std.Err.", "t-stat", "p-value"))
    cat(strrep("-", 60), "\n")

    for (i in seq_len(nrows)) {
      est_val <- est[i, t_idx]
      se_val <- se[i, t_idx]
      tstat <- est_val / se_val
      pval <- 2 * (1 - stats::pnorm(abs(tstat)))
      stars <- get_stars(pval)

      cat(sprintf("  %12s %12.4f %10.4f %10.3f %9.4f%s\n",
                  row_names[i], est_val, se_val, tstat, pval, stars))
    }
    cat("\n")
  }
}


#' Get Significance Stars
#'
#' @keywords internal
get_stars <- function(pval) {
  if (is.na(pval)) return("")
  if (pval < 0.01) return("***")
  if (pval < 0.05) return("**")
  if (pval < 0.10) return("*")
  return("")
}


#' Coefficients Method for QARDL
#'
#' Extract coefficients from a QARDL model.
#'
#' @param object An object of class \code{"qardl"}.
#' @param type Character. Which coefficients to extract: \code{"beta"},
#'   \code{"phi"}, \code{"gamma"}, or \code{"all"}. Default is \code{"all"}.
#' @param ... Additional arguments (unused).
#'
#' @return Matrix or list of coefficient matrices.
#'
#' @export
coef.qardl <- function(object, type = c("all", "beta", "phi", "gamma"), ...) {

  type <- match.arg(type)

  if (type == "beta") {
    return(object$beta)
  } else if (type == "phi") {
    return(object$phi)
  } else if (type == "gamma") {
    return(object$gamma)
  } else {
    return(list(
      beta = object$beta,
      phi = object$phi,
      gamma = object$gamma
    ))
  }
}


#' Variance-Covariance Method for QARDL
#'
#' Extract variance-covariance matrices from a QARDL model.
#'
#' @param object An object of class \code{"qardl"}.
#' @param type Character. Which covariance to extract: \code{"beta"},
#'   \code{"phi"}, \code{"gamma"}, or \code{"all"}. Default is \code{"all"}.
#' @param ... Additional arguments (unused).
#'
#' @return Array or list of covariance arrays.
#'
#' @export
vcov.qardl <- function(object, type = c("all", "beta", "phi", "gamma"), ...) {

  type <- match.arg(type)

  if (type == "beta") {
    return(object$beta_cov)
  } else if (type == "phi") {
    return(object$phi_cov)
  } else if (type == "gamma") {
    return(object$gamma_cov)
  } else {
    return(list(
      beta = object$beta_cov,
      phi = object$phi_cov,
      gamma = object$gamma_cov
    ))
  }
}
