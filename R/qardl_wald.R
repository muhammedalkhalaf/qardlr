#' Wald Tests for QARDL Parameter Constancy
#'
#' Performs Wald tests for parameter constancy across quantiles in a QARDL model.
#' Tests whether parameters are equal across different quantile levels.
#'
#' @param object An object of class \code{"qardl"}.
#' @param type Character string specifying which parameters to test:
#'   \code{"all"} (default), \code{"beta"} (long-run), \code{"phi"} (AR),
#'   \code{"gamma"} (short-run impact), or \code{"rho"} (ECM speed of adjustment).
#' @param pairwise Logical. If \code{TRUE}, perform pairwise tests between
#'   adjacent quantiles. Default is \code{FALSE}.
#'
#' @return An object of class \code{"qardl_wald"} containing:
#' \describe{
#'   \item{tests}{Data frame of test results with columns: test, statistic, df, pvalue}
#'   \item{pairwise_tests}{Data frame of pairwise test results (if pairwise = TRUE)}
#'   \item{type}{Type of test performed}
#'   \item{tau}{Vector of quantiles}
#' }
#'
#' @details
#' The Wald test statistic is computed as:
#' \deqn{W = (R\hat{\theta} - r)' [R \hat{V} R']^{-1} (R\hat{\theta} - r) \sim \chi^2(q)}
#'
#' where \eqn{R} is a restriction matrix testing equality across quantiles,
#' \eqn{\hat{\theta}} is the vector of parameter estimates, and \eqn{\hat{V}}
#' is the estimated covariance matrix.
#'
#' @references
#' Cho, J.S., Kim, T.-H., & Shin, Y. (2015). Quantile cointegration in the
#' autoregressive distributed-lag modeling framework. \emph{Journal of
#' Econometrics}, 188(1), 281-300. \doi{10.1016/j.jeconom.2015.01.003}
#'
#' @seealso \code{\link{qardl}}, \code{\link{print.qardl_wald}}
#'
#' @examples
#' data(qardl_sim)
#' fit <- qardl(y ~ x1 + x2, data = qardl_sim, tau = c(0.25, 0.50, 0.75), p = 2, q = 2)
#' wald_results <- qardl_wald(fit)
#' print(wald_results)
#'
#' # Pairwise tests
#' wald_pairwise <- qardl_wald(fit, pairwise = TRUE)
#' print(wald_pairwise)
#'
#' @export
qardl_wald <- function(object, type = c("all", "beta", "phi", "gamma", "rho"),
                       pairwise = FALSE) {

  if (!inherits(object, "qardl")) {
    stop("'object' must be of class 'qardl'", call. = FALSE)
  }

  type <- match.arg(type)
  tau <- object$tau
  ntau <- length(tau)

  if (ntau < 2) {
    stop("At least two quantiles are required for constancy tests", call. = FALSE)
  }

  tests <- list()
  pairwise_tests <- NULL

  # Test beta constancy
  if (type %in% c("all", "beta")) {
    beta_test <- wald_constancy_test(object$beta, object$beta_cov,
                                      object$nobs, "beta")
    tests$beta <- beta_test
  }

  # Test phi constancy
  if (type %in% c("all", "phi")) {
    phi_test <- wald_constancy_test(object$phi, object$phi_cov,
                                     object$nobs, "phi")
    tests$phi <- phi_test
  }

  # Test gamma constancy
  if (type %in% c("all", "gamma")) {
    gamma_test <- wald_constancy_test(object$gamma, object$gamma_cov,
                                       object$nobs, "gamma")
    tests$gamma <- gamma_test
  }

  # Test rho constancy (scalar per quantile)
  if (type %in% c("all", "rho")) {
    rho_mat <- matrix(object$rho, nrow = 1)
    rho_cov <- array(object$rho_se^2, dim = c(1, 1, ntau))
    rho_test <- wald_constancy_test(rho_mat, rho_cov, object$nobs, "rho")
    tests$rho <- rho_test
  }

  # Combine into data frame
  test_df <- do.call(rbind, lapply(names(tests), function(nm) {
    t <- tests[[nm]]
    data.frame(
      parameter = nm,
      statistic = t$statistic,
      df = t$df,
      pvalue = t$pvalue,
      decision = ifelse(t$pvalue < 0.01, "Reject***",
                        ifelse(t$pvalue < 0.05, "Reject**",
                               ifelse(t$pvalue < 0.10, "Reject*", "Fail to reject"))),
      stringsAsFactors = FALSE
    )
  }))

  # Pairwise tests if requested
  if (pairwise) {
    pairwise_tests <- wald_pairwise_tests(object, type)
  }

  result <- list(
    tests = test_df,
    pairwise_tests = pairwise_tests,
    type = type,
    tau = tau
  )

  class(result) <- "qardl_wald"
  return(result)
}


#' Wald Constancy Test
#'
#' Internal function to perform Wald test for parameter constancy.
#'
#' @param params Parameter matrix (dim x ntau).
#' @param cov_array Covariance array (dim x dim x ntau).
#' @param nobs Number of observations.
#' @param param_name Name of parameter for labeling.
#'
#' @return List with statistic, df, pvalue.
#'
#' @keywords internal
wald_constancy_test <- function(params, cov_array, nobs, param_name) {

  dim_param <- nrow(params)
  ntau <- ncol(params)

  if (ntau < 2) {
    return(list(statistic = NA, df = NA, pvalue = NA))
  }

  # Number of restrictions: (ntau - 1) * dim_param
  # Testing: param(tau_i) = param(tau_{i+1}) for all i

  n_restr <- (ntau - 1) * dim_param

  # Build R matrix and compute Wald statistic
  # Stack parameters: vec(params) = [param_tau1; param_tau2; ...]

  theta <- as.vector(params)  # Column-major stacking

  # Build block-diagonal covariance (assuming independence across quantiles)
  V <- matrix(0, nrow = length(theta), ncol = length(theta))
  for (t_idx in seq_len(ntau)) {
    idx <- ((t_idx - 1) * dim_param + 1):(t_idx * dim_param)
    V[idx, idx] <- cov_array[, , t_idx]
  }

  # Build R matrix: differences between adjacent quantiles
  R <- matrix(0, nrow = n_restr, ncol = length(theta))
  row_idx <- 1
  for (t_idx in seq_len(ntau - 1)) {
    for (d in seq_len(dim_param)) {
      # theta_{t_idx, d} - theta_{t_idx+1, d} = 0
      col1 <- (t_idx - 1) * dim_param + d
      col2 <- t_idx * dim_param + d
      R[row_idx, col1] <- 1
      R[row_idx, col2] <- -1
      row_idx <- row_idx + 1
    }
  }

  # r vector (all zeros for equality test)
  r <- rep(0, n_restr)

  # Wald statistic: W = (R*theta - r)' * inv(R*V*R') * (R*theta - r)
  diff <- R %*% theta - r
  RVR <- R %*% V %*% t(R)

  # Regularize if needed
  RVR <- RVR + diag(1e-12, nrow(RVR))

  RVR_inv <- tryCatch(
    solve(RVR),
    error = function(e) MASS::ginv(RVR)
  )

  W <- as.numeric(t(diff) %*% RVR_inv %*% diff)

  if (W < 0) W <- abs(W)

  pval <- stats::pchisq(W, df = n_restr, lower.tail = FALSE)

  return(list(
    statistic = W,
    df = n_restr,
    pvalue = pval
  ))
}


#' Pairwise Wald Tests
#'
#' Internal function for pairwise parameter equality tests.
#'
#' @param object QARDL object.
#' @param type Parameter type.
#'
#' @return Data frame of pairwise test results.
#'
#' @keywords internal
wald_pairwise_tests <- function(object, type) {

  tau <- object$tau
  ntau <- length(tau)
  k <- object$k
  p <- object$p
  indepvars <- object$indepvars

  results <- list()

  # Beta pairwise tests
  if (type %in% c("all", "beta")) {
    for (v_idx in seq_len(k)) {
      var_name <- indepvars[v_idx]
      for (i in seq_len(ntau - 1)) {
        for (j in (i + 1):ntau) {
          b_i <- object$beta[v_idx, i]
          b_j <- object$beta[v_idx, j]
          diff <- b_i - b_j

          v_ii <- object$beta_cov[v_idx, v_idx, i]
          v_jj <- object$beta_cov[v_idx, v_idx, j]
          var_diff <- v_ii + v_jj

          if (var_diff > 1e-15) {
            W <- diff^2 / var_diff
            pval <- stats::pchisq(W, df = 1, lower.tail = FALSE)

            results[[length(results) + 1]] <- data.frame(
              parameter = "beta",
              variable = var_name,
              tau_i = tau[i],
              tau_j = tau[j],
              statistic = W,
              df = 1,
              pvalue = pval,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  # Phi pairwise tests
  if (type %in% c("all", "phi")) {
    for (lag in seq_len(p)) {
      for (i in seq_len(ntau - 1)) {
        for (j in (i + 1):ntau) {
          phi_i <- object$phi[lag, i]
          phi_j <- object$phi[lag, j]
          diff <- phi_i - phi_j

          v_ii <- object$phi_cov[lag, lag, i]
          v_jj <- object$phi_cov[lag, lag, j]
          var_diff <- v_ii + v_jj

          if (var_diff > 1e-15) {
            W <- diff^2 / var_diff
            pval <- stats::pchisq(W, df = 1, lower.tail = FALSE)

            results[[length(results) + 1]] <- data.frame(
              parameter = "phi",
              variable = paste0("lag_", lag),
              tau_i = tau[i],
              tau_j = tau[j],
              statistic = W,
              df = 1,
              pvalue = pval,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  # Gamma pairwise tests
  if (type %in% c("all", "gamma")) {
    for (v_idx in seq_len(k)) {
      var_name <- indepvars[v_idx]
      for (i in seq_len(ntau - 1)) {
        for (j in (i + 1):ntau) {
          g_i <- object$gamma[v_idx, i]
          g_j <- object$gamma[v_idx, j]
          diff <- g_i - g_j

          v_ii <- object$gamma_cov[v_idx, v_idx, i]
          v_jj <- object$gamma_cov[v_idx, v_idx, j]
          var_diff <- v_ii + v_jj

          if (var_diff > 1e-15) {
            W <- diff^2 / var_diff
            pval <- stats::pchisq(W, df = 1, lower.tail = FALSE)

            results[[length(results) + 1]] <- data.frame(
              parameter = "gamma",
              variable = var_name,
              tau_i = tau[i],
              tau_j = tau[j],
              statistic = W,
              df = 1,
              pvalue = pval,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  if (length(results) > 0) {
    return(do.call(rbind, results))
  } else {
    return(NULL)
  }
}


#' Print QARDL Wald Test Results
#'
#' @param x Object of class \code{"qardl_wald"}.
#' @param digits Number of decimal places. Default is 4.
#' @param ... Additional arguments (unused).
#'
#' @return Invisible \code{x}.
#'
#' @export
print.qardl_wald <- function(x, digits = 4, ...) {

  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("  Wald Tests for Parameter Constancy Across Quantiles\n")
  cat("  H0: parameter(tau_i) = parameter(tau_j) for all i, j\n")
  cat(strrep("=", 70), "\n\n")

  cat("  Quantiles: ", paste(round(x$tau, 2), collapse = ", "), "\n\n")

  # Main tests
  cat(strrep("-", 70), "\n")
  cat(sprintf("  %-15s %12s %8s %12s %15s\n",
              "Parameter", "Wald stat", "df", "p-value", "Decision"))
  cat(strrep("-", 70), "\n")

  for (i in seq_len(nrow(x$tests))) {
    row <- x$tests[i, ]
    pval_str <- format(round(row$pvalue, digits), nsmall = digits)
    cat(sprintf("  %-15s %12.3f %8d %12s %15s\n",
                row$parameter, row$statistic, row$df, pval_str, row$decision))
  }

  cat(strrep("-", 70), "\n")
  cat("  *** p<0.01, ** p<0.05, * p<0.10\n\n")

  # Pairwise tests if present
  if (!is.null(x$pairwise_tests) && nrow(x$pairwise_tests) > 0) {
    cat(strrep("=", 70), "\n")
    cat("  Pairwise Equality Tests\n")
    cat(strrep("=", 70), "\n\n")

    cat(sprintf("  %-10s %-12s %8s %8s %10s %6s %10s\n",
                "Parameter", "Variable", "tau_i", "tau_j", "Wald", "df", "p-value"))
    cat(strrep("-", 70), "\n")

    for (i in seq_len(nrow(x$pairwise_tests))) {
      row <- x$pairwise_tests[i, ]
      stars <- ifelse(row$pvalue < 0.01, "***",
                      ifelse(row$pvalue < 0.05, "**",
                             ifelse(row$pvalue < 0.10, "*", "")))
      cat(sprintf("  %-10s %-12s %8.2f %8.2f %10.3f %6d %9.4f%s\n",
                  row$parameter, row$variable, row$tau_i, row$tau_j,
                  row$statistic, row$df, row$pvalue, stars))
    }

    cat(strrep("-", 70), "\n")
  }

  invisible(x)
}
