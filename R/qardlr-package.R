#' qardlr: Quantile Autoregressive Distributed Lag Model
#'
#' @description
#' The qardlr package implements the Quantile Autoregressive Distributed Lag
#' (QARDL) model of Cho, Kim and Shin (2015). It provides tools for estimating
#' quantile-specific long-run equilibrium relationships and short-run dynamics.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{qardl}}: Estimate QARDL model
#'   \item \code{\link{qardl_wald}}: Wald tests for parameter constancy
#'   \item \code{\link{qardl_rolling}}: Rolling window QARDL estimation
#'   \item \code{\link{qardl_simulate}}: Monte Carlo simulation
#'   \item \code{\link{qardl_table}}: Publication-ready tables
#'   \item \code{\link{qardl_bic_select}}: BIC-based lag selection
#' }
#'
#' @section Key Features:
#' \itemize{
#'   \item Quantile regression across multiple tau values
#'   \item BIC-based automatic lag selection (p, q)
#'   \item Error Correction Model (ECM) parameterization
#'   \item Long-run (beta), short-run AR (phi), and impact (gamma) parameters
#'   \item Wald tests for parameter constancy across quantiles
#'   \item Rolling/recursive QARDL estimation
#'   \item Monte Carlo simulation for finite-sample properties
#'   \item Publication-ready output tables (text, LaTeX, HTML)
#' }
#'
#' @references
#' Cho, J.S., Kim, T.-H., and Shin, Y. (2015). Quantile cointegration in the
#' autoregressive distributed-lag modeling framework. \emph{Journal of
#' Econometrics}, 188(1), 281-300. \doi{10.1016/j.jeconom.2015.01.003}
#'
#' @keywords internal
"_PACKAGE"

#' @name qardlr-package
#' @aliases qardlr-package NULL
#'
#' @importFrom quantreg rq
#' @importFrom stats coef fitted model.frame model.matrix model.response
#' @importFrom stats na.pass pchisq pnorm residuals rnorm sd
#' @importFrom MASS ginv
#' @importFrom grDevices rainbow
#' @importFrom graphics abline lines legend plot
#'
NULL
