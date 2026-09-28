test_that("qardl works with basic input", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 + 0.8 * x + rnorm(n)
  data <- data.frame(y = y, x = x)

  fit <- qardl(y ~ x, data = data, tau = c(0.25, 0.5, 0.75), p = 1, q = 1)

  expect_s3_class(fit, "qardl")
  expect_equal(length(fit$tau), 3)
  expect_equal(fit$p, 1)
  expect_equal(fit$q, 1)
})

test_that("qardl_bic_select works", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 + 0.8 * x + rnorm(n)
  X <- matrix(x, ncol = 1)
  colnames(X) <- "x"

  bic_sel <- qardl_bic_select(y, X, pmax = 3, qmax = 3)

  expect_true(bic_sel$p_opt >= 1 && bic_sel$p_opt <= 3)
  expect_true(bic_sel$q_opt >= 1 && bic_sel$q_opt <= 3)
  expect_true(!is.na(bic_sel$bic_min))
})

test_that("qardl_wald returns correct structure", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 + 0.8 * x + rnorm(n)
  data <- data.frame(y = y, x = x)

  fit <- qardl(y ~ x, data = data, tau = c(0.25, 0.5, 0.75), p = 1, q = 1)
  wald <- qardl_wald(fit)

  expect_s3_class(wald, "qardl_wald")
  expect_true(!is.null(wald$tests))
  expect_true(nrow(wald$tests) > 0)
})

test_that("qardl_rolling works", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 + 0.8 * x + rnorm(n)
  data <- data.frame(y = y, x = x)

  roll <- qardl_rolling(y ~ x, data = data, tau = c(0.5), p = 1, q = 1, window = 40)

  expect_s3_class(roll, "qardl_rolling")
  expect_true(roll$n_windows > 0)
})

test_that("qardl_simulate works", {
  set.seed(123)

  mc <- qardl_simulate(nobs = 50, reps = 10, tau = c(0.5), p = 1, q = 1, k = 1)

  expect_s3_class(mc, "qardl_mc")
  expect_equal(mc$reps, 10)
  expect_equal(mc$nobs, 50)
})

test_that("coef and vcov methods work", {
  set.seed(123)
  n <- 100
  x <- cumsum(rnorm(n))
  y <- 0.5 + 0.8 * x + rnorm(n)
  data <- data.frame(y = y, x = x)

  fit <- qardl(y ~ x, data = data, tau = c(0.25, 0.5, 0.75), p = 1, q = 1)

  beta <- coef(fit, type = "beta")
  expect_true(is.matrix(beta))

  vcov_beta <- vcov(fit, type = "beta")
  expect_true(is.array(vcov_beta))
})

test_that("gamma0 picks the coefficients on x_t and beta = gamma / (1 - sum(phi))", {
  data(qardl_sim)
  fit <- qardl(y ~ x1 + x2, data = qardl_sim, tau = c(0.25, 0.5, 0.75), p = 2, q = 2)
  expect_equal(unname(fit$gamma0), unname(fit$raw_coefs[c("x1_lag0", "x2_lag0"), ]))
  expect_equal(unname(fit$gamma),
               unname(fit$raw_coefs[c("x1_lag0", "x2_lag0"), ] +
                        fit$raw_coefs[c("x1_lag1", "x2_lag1"), ]))
  expect_equal(fit$beta, sweep(fit$gamma, 2, 1 - colSums(fit$phi), "/"))
  ## Cho, Kim and Shin (2015) covariance of beta across quantiles:
  ## (min(ti,tj) - ti tj) / (fi fj (1-sum phi_i)(1-sum phi_j)) (X'PX)^-1
  V <- fit$cov_joint$beta
  om <- unname(1 - colSums(fit$phi))
  r12 <- V[1, 3] / V[1, 1]
  expect_equal(r12, (0.25 - 0.25 * 0.5) / (0.25 * 0.75) *
                 fit$fhat[1] * om[1] / (fit$fhat[2] * om[2]))
})
