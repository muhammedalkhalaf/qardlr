#' Simulated QARDL Dataset
#'
#' A simulated dataset for demonstrating QARDL estimation. The data
#' is generated from a QARDL(2,2) process with two covariates.
#'
#' @format A data frame with 200 observations and 3 variables:
#' \describe{
#'   \item{y}{Dependent variable generated from QARDL process}
#'   \item{x1}{First covariate (I(1) random walk)}
#'   \item{x2}{Second covariate (I(1) random walk)}
#' }
#'
#' @details
#' The data generating process follows:
#' \deqn{y_t = 0.4 y_{t-1} + 0.2 y_{t-2} + 0.5 x_{1t} + 0.3 x_{2t} + u_t}
#'
#' where \eqn{u_t \sim N(0, 1)} and \eqn{x_{it}} are independent random walks.
#'
#' True parameters:
#' \itemize{
#'   \item \eqn{\phi_1 = 0.4}, \eqn{\phi_2 = 0.2}
#'   \item \eqn{\gamma_1 = 0.5}, \eqn{\gamma_2 = 0.3}
#'   \item \eqn{\beta_1 = 0.5/(1-0.6) = 1.25}, \eqn{\beta_2 = 0.3/(1-0.6) = 0.75}
#' }
#'
#' @references
#' Cho, J.S., Kim, T.-H., & Shin, Y. (2015). Quantile cointegration in the
#' autoregressive distributed-lag modeling framework. \emph{Journal of
#' Econometrics}, 188(1), 281-300. \doi{10.1016/j.jeconom.2015.01.003}
#'
#' @examples
#' data(qardl_sim)
#' head(qardl_sim)
#' summary(qardl_sim)
#'
#' # Estimate QARDL model
#' fit <- qardl(y ~ x1 + x2, data = qardl_sim, tau = c(0.25, 0.50, 0.75), p = 2, q = 2)
#' summary(fit)
#'
"qardl_sim"
