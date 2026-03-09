#' BIC-Based Lag Order Selection for QARDL
#'
#' Automatically selects optimal lag orders (p, q) for the QARDL model
#' using the Bayesian Information Criterion (BIC) evaluated at the median
#' quantile (tau = 0.5).
#'
#' @param y Numeric vector of dependent variable.
#' @param X Matrix of covariates.
#' @param pmax Integer. Maximum AR lag order to consider. Default is 7.
#' @param qmax Integer. Maximum distributed lag order to consider. Default is 7.
#' @param constant Logical. Include intercept. Default is \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{p_opt}{Optimal AR lag order}
#'   \item{q_opt}{Optimal distributed lag order}
#'   \item{bic_grid}{Matrix of BIC values (pmax x qmax)}
#'   \item{bic_min}{Minimum BIC value}
#' }
#'
#' @details
#' The BIC is computed using the Schwarz criterion at the median quantile:
#' \deqn{BIC(p, q) = \log(\hat{\sigma}^2_{\tau=0.5}) + \frac{k_{pq} \log(n)}{n}}
#'
#' where \eqn{k_{pq}} is the number of parameters (p AR terms + q*k impact terms
#' + constant) and \eqn{\hat{\sigma}^2} is the estimated residual variance.
#'
#' @references
#' Cho, J.S., Kim, T.-H., & Shin, Y. (2015). Quantile cointegration in the
#' autoregressive distributed-lag modeling framework. \emph{Journal of
#' Econometrics}, 188(1), 281-300. \doi{10.1016/j.jeconom.2015.01.003}
#'
#' @seealso \code{\link{qardl}}
#'
#' @examples
#' data(qardl_sim)
#' y <- qardl_sim$y
#' X <- as.matrix(qardl_sim[, c("x1", "x2")])
#' bic_result <- qardl_bic_select(y, X, pmax = 5, qmax = 5)
#' print(bic_result$bic_grid)
#'
#' @export
qardl_bic_select <- function(y, X, pmax = 7L, qmax = 7L, constant = TRUE) {

  n <- length(y)
  k <- ncol(X)

  pmax <- as.integer(pmax)
  qmax <- as.integer(qmax)

  if (pmax < 1) pmax <- 1L
  if (qmax < 1) qmax <- 1L

  # BIC grid

  bic_grid <- matrix(NA, nrow = pmax, ncol = qmax)
  rownames(bic_grid) <- paste0("p=", seq_len(pmax))
  colnames(bic_grid) <- paste0("q=", seq_len(qmax))

  # Estimate at median (tau = 0.5) for each (p, q) combination
  for (p in seq_len(pmax)) {
    for (q in seq_len(qmax)) {
      maxlag <- max(p, q)

      # Check if we have enough observations
      neff <- n - maxlag
      num_params <- p + q * k + ifelse(constant, 1, 0)

      if (neff <= num_params + 5) {
        bic_grid[p, q] <- Inf
        next
      }

      # Build design matrix
      y_eff <- y[(maxlag + 1):n]

      design_cols <- list()

      # AR lags
      for (i in seq_len(p)) {
        lag_idx <- (maxlag + 1 - i):(n - i)
        design_cols[[paste0("y_lag", i)]] <- y[lag_idx]
      }

      # X and its lags
      for (j in seq_len(q)) {
        lag <- j - 1
        for (kk in seq_len(k)) {
          if (lag == 0) {
            idx <- (maxlag + 1):n
          } else {
            idx <- (maxlag + 1 - lag):(n - lag)
          }
          design_cols[[paste0("x", kk, "_lag", lag)]] <- X[idx, kk]
        }
      }

      Z <- do.call(cbind, design_cols)
      if (constant) {
        Z <- cbind(1, Z)
      }

      # Quantile regression at median
      fit <- tryCatch(
        quantreg::rq(y_eff ~ Z - 1, tau = 0.5),
        error = function(e) NULL
      )

      if (is.null(fit)) {
        bic_grid[p, q] <- Inf
        next
      }

      # Compute BIC
      resid <- residuals(fit)
      sigma2 <- mean(resid^2)

      if (sigma2 <= 0) sigma2 <- 1e-10

      bic_val <- log(sigma2) + (num_params * log(neff)) / neff
      bic_grid[p, q] <- bic_val
    }
  }

  # Find optimal (p, q)
  min_idx <- which(bic_grid == min(bic_grid, na.rm = TRUE), arr.ind = TRUE)
  p_opt <- min_idx[1, 1]
  q_opt <- min_idx[1, 2]
  bic_min <- bic_grid[p_opt, q_opt]

  return(list(
    p_opt = p_opt,
    q_opt = q_opt,
    bic_grid = bic_grid,
    bic_min = bic_min
  ))
}


#' Print BIC Grid
#'
#' Prints a formatted BIC grid for lag selection.
#'
#' @param bic_result Result from \code{qardl_bic_select}.
#' @param digits Number of decimal places. Default is 3.
#'
#' @return Invisible \code{NULL}. Called for side effect of printing.
#'
#' @examples
#' data(qardl_sim)
#' y <- qardl_sim$y
#' X <- as.matrix(qardl_sim[, c("x1", "x2")])
#' bic_result <- qardl_bic_select(y, X, pmax = 4, qmax = 4)
#' print_bic_grid(bic_result)
#'
#' @export
print_bic_grid <- function(bic_result, digits = 3) {

  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("  QARDL Lag Order Selection (BIC)\n")
  cat(strrep("=", 70), "\n\n")

  bic_grid <- bic_result$bic_grid
  p_opt <- bic_result$p_opt
  q_opt <- bic_result$q_opt

  pmax <- nrow(bic_grid)
  qmax <- ncol(bic_grid)

  # Header
  cat("  BIC Grid: rows = p (AR lags), columns = q (DL lags)\n")
  cat(strrep("-", 70), "\n")

  # Column headers
  cat(sprintf("%8s", "p \\ q"))
  for (j in seq_len(qmax)) {
    cat(sprintf("%12s", paste0("q=", j)))
  }
  cat("\n")
  cat(strrep("-", 70), "\n")

  # Grid values
  for (i in seq_len(pmax)) {
    cat(sprintf("%8s", paste0("p=", i)))
    for (j in seq_len(qmax)) {
      val <- bic_grid[i, j]
      if (is.infinite(val)) {
        cat(sprintf("%12s", "Inf"))
      } else if (i == p_opt && j == q_opt) {
        cat(sprintf("%11s*", format(round(val, digits), nsmall = digits)))
      } else {
        cat(sprintf("%12s", format(round(val, digits), nsmall = digits)))
      }
    }
    cat("\n")
  }

  cat(strrep("-", 70), "\n")
  cat(sprintf("  Optimal: p = %d, q = %d (min BIC = %.3f)\n",
              p_opt, q_opt, bic_result$bic_min))
  cat("  * denotes minimum BIC\n")
  cat(strrep("=", 70), "\n\n")

}
