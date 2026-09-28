#' Asymptotic Covariances of the QARDL Estimators
#'
#' Computes the covariance matrices of the QARDL estimators within and
#' across quantiles following Cho, Kim and Shin (2015): Theorems 1 and 3
#' (autoregressive parameters \eqn{\phi} and the level coefficient
#' \eqn{\gamma}), Theorems 2 and 4 (long-run parameter \eqn{\beta}), and
#' the estimators \eqn{\hat f_\tau} and \eqn{\hat\Pi_n(\tau)} of their
#' Section 5 and equation (13).
#'
#' @param y Dependent variable.
#' @param X Matrix of covariates.
#' @param est Output of \code{qardl_estimate}.
#' @param beta Long-run parameters (k x ntau).
#' @param tau Quantiles.
#' @param constant Logical for intercept inclusion.
#'
#' @return A list with joint covariance matrices \code{phi}, \code{gamma},
#'   \code{beta} and \code{rho}, ordered by quantile (all parameters at
#'   the first quantile, then the second, ...), and the density estimates
#'   \code{fhat}.
#'
#' @details
#' With \eqn{X_t} the covariate levels, \eqn{\tilde W_t = (1, \Delta
#' X_t', \dots, \Delta X_{t-q+2}')'} and \eqn{P} the projection off
#' \eqn{\tilde W}, the long-run covariance between quantiles
#' \eqn{\tau_i} and \eqn{\tau_j} is
#' \deqn{\frac{\min(\tau_i,\tau_j) - \tau_i\tau_j}{f_{\tau_i} f_{\tau_j}
#' (1 - \sum\phi(\tau_i))(1 - \sum\phi(\tau_j))} (X'PX)^{-1}.}
#' \eqn{\hat K_{t,i}(\tau)} is the quantile regression error of
#' \eqn{Y_{t-i} - X_t'\hat\beta(\tau)} on \eqn{\tilde W_t}, and
#' \eqn{L(\tau_i,\tau_j) = n^{-1}\hat K(\tau_i)' P \hat K(\tau_j)}; the
#' covariance of \eqn{\hat\phi(\tau_i)} and \eqn{\hat\phi(\tau_j)} is
#' \eqn{n^{-1}[\min(\tau_i,\tau_j) - \tau_i\tau_j] f_{\tau_i}^{-1}
#' f_{\tau_j}^{-1} L(\tau_i,\tau_i)^{-1} L(\tau_i,\tau_j)
#' L(\tau_j,\tau_j)^{-1}}, and \eqn{\hat\gamma(\tau) = \hat\beta(\tau)(1 -
#' \sum\hat\phi(\tau))} has covariance \eqn{\beta_i \iota'\,
#' \mathrm{Cov}(\hat\phi_i, \hat\phi_j)\,\iota\, \beta_j'} (Theorem 1(ii)).
#' The density is the Gaussian kernel estimate at zero of the quantile
#' regression residuals with the Bofinger (1975) bandwidth.
#'
#' @references
#' Cho, J.S., Kim, T.-H. and Shin, Y. (2015). Quantile cointegration in the
#' autoregressive distributed-lag modeling framework. \emph{Journal of
#' Econometrics}, 188(1), 281-300. \doi{10.1016/j.jeconom.2015.05.003}
#'
#' @keywords internal
cks_covariance <- function(y, X, est, beta, tau, constant = TRUE) {
  n_all <- length(y)
  k <- ncol(X)
  p <- est$p
  q <- est$q
  ntau <- length(tau)
  maxlag <- max(p, q)
  rows <- (maxlag + 1):n_all
  n <- length(rows)

  Xt <- X[rows, , drop = FALSE]
  dX <- rbind(NA_real_, diff(X))
  Wt <- if (constant) matrix(1, n, 1) else matrix(0, n, 0)
  if (q >= 2) {
    for (j in 0:(q - 2)) Wt <- cbind(Wt, dX[rows - j, , drop = FALSE])
  }
  proj <- function(A) {
    if (ncol(Wt) == 0) return(A)
    A - Wt %*% qr.coef(qr(Wt), A)
  }

  ## density at the tau-quantile (Bofinger bandwidth, Gaussian kernel)
  fhat <- vapply(seq_len(ntau), function(i) {
    u <- as.numeric(stats::residuals(est$qr_fits[[i]]))
    z <- stats::qnorm(tau[i])
    h <- n^(-1/5) * (4.5 * stats::dnorm(z)^4 / (2 * z^2 + 1)^2)^(1/5)
    mean(stats::dnorm(-u / h)) / h
  }, numeric(1))

  ## K-hat(tau): quantile regression errors of Y_{t-i} - X_t' beta(tau) on W
  Kh <- lapply(seq_len(ntau), function(i) {
    sapply(seq_len(p), function(j) {
      v <- y[rows - j] - as.numeric(Xt %*% beta[, i])
      if (ncol(Wt) == 0) return(v)
      as.numeric(stats::residuals(quantreg::rq.fit(Wt, v, tau = tau[i])))
    })
  })
  Kh <- lapply(Kh, function(m) matrix(m, n, p))
  PK <- lapply(Kh, proj)
  L <- function(i, j) crossprod(Kh[[i]], PK[[j]]) / n

  one_minus <- 1 - colSums(est$phi)
  XPX_inv <- solve(crossprod(Xt, proj(Xt)))
  iota <- rep(1, p)

  phiV <- matrix(0, p * ntau, p * ntau)
  betaV <- matrix(0, k * ntau, k * ntau)
  gamV <- matrix(0, k * ntau, k * ntau)
  rhoV <- matrix(0, ntau, ntau)
  Linv <- lapply(seq_len(ntau), function(i) solve(L(i, i)))
  for (i in seq_len(ntau)) for (j in seq_len(ntau)) {
    w <- min(tau[i], tau[j]) - tau[i] * tau[j]
    Xi <- w / (fhat[i] * fhat[j]) * Linv[[i]] %*% L(i, j) %*% Linv[[j]] / n
    pi_ <- (i - 1) * p + seq_len(p); pj <- (j - 1) * p + seq_len(p)
    ki <- (i - 1) * k + seq_len(k); kj <- (j - 1) * k + seq_len(k)
    phiV[pi_, pj] <- Xi
    rhoV[i, j] <- as.numeric(t(iota) %*% Xi %*% iota)
    gamV[ki, kj] <- (beta[, i] %o% beta[, j]) * rhoV[i, j]
    betaV[ki, kj] <- w / (fhat[i] * fhat[j] * one_minus[i] * one_minus[j]) * XPX_inv
  }

  list(phi = phiV, gamma = gamV, beta = betaV, rho = rhoV, fhat = fhat)
}


#' Extract a Per-Quantile Covariance Array from a Joint Covariance Matrix
#' @param V Joint covariance matrix.
#' @param d Dimension of the parameter vector at each quantile.
#' @param ntau Number of quantiles.
#' @return Array d x d x ntau.
#' @keywords internal
cks_blocks <- function(V, d, ntau) {
  out <- array(NA_real_, dim = c(d, d, ntau))
  for (i in seq_len(ntau)) {
    idx <- (i - 1) * d + seq_len(d)
    out[, , i] <- V[idx, idx, drop = FALSE]
  }
  out
}
