## Inference notice for version 1.0.1 and earlier

**The point estimates in qardlr 1.0.1 are usable. The standard errors, t ratios, confidence intervals and Wald p-values are not.** If you have published inference from this package, please read on.

**1. The long-run standard errors have no asymptotic justification.** In the Cho, Kim and Shin (2015) framework the long-run coefficient beta(tau) converges at rate n with a mixture-normal limit, standardised by the random matrix M = n^-2 X'[I - W(W'W)^-1 W']X. Version 1.0.1 instead applies a delta method to summary.rq with se = "nid", which is the stationary sandwich derived for rate-root-n asymptotics with stationary regressors. It does not converge to the right object at the right rate, and because the projection of x_t on the lagged differences is absent, the reported standard errors will usually be too small when x is endogenous. Long-run coefficients therefore look more significant than they are.

**2. The Wald tests for constancy across quantiles assume independence.** Quantile-regression estimates from one sample are positively correlated, with covariance factor min(tau_i, tau_j) - tau_i tau_j. The correct variance of a difference is V_ii + V_jj - 2 V_ij; the package uses V_ii + V_jj. For tau = 0.25 against tau = 0.50 the correct denominator is about 43 percent of the one used, so the statistic is roughly 43 percent of its correct value and the test under-rejects constancy. Since rejection is the evidence for quantile heterogeneity, the package systematically produces evidence against the phenomenon it exists to detect.

**3. beta_se and beta_cov are different variance concepts for the same parameter.** beta_se includes a term for the variance of rho; beta_cov does not. Printed t ratios and Wald tests therefore use inconsistent variances.

**4. The short-run gamma is read from the wrong design columns whenever the number of covariates differs from the lag order.** The design is built lag-major and the extraction uses a variable-major stride, so for k = 2, q = 3 the reported contemporaneous coefficients are x1_lag0 and x2_lag1. The long-run beta is not affected.

**5. The model is not the Cho, Kim and Shin QARDL.** Their specification includes the differences Delta x_t through Delta x_t-q+1 alongside the level x_t; qardlr fits a plain levels ARDL. The difference block is what controls the endogeneity of the I(1) regressors, and it is the reason the asymptotic theory holds.

**What is reliable in 1.0.1**: the point estimates of beta and rho. On a design with a true long-run coefficient of 0.80 and a true adjustment speed of -0.35, the package returned 0.792 / 0.786 / 0.806 and -0.373 / -0.346 / -0.351 at tau = 0.25 / 0.50 / 0.75.

**Recommended for now**: report the point estimates, do not report the standard errors or the Wald tests, and set p and q explicitly rather than relying on the BIC selector, which uses a different criterion, a different design and a different estimator from the Stata original. Version 2.0.0 will implement the Cho, Kim and Shin covariance properly.

---

# qardlr: Quantile Autoregressive Distributed Lag Model

R implementation of the Quantile Autoregressive Distributed Lag (QARDL) model by Cho, Kim & Shin (2015).

## Overview

The qardlr package provides tools for estimating quantile-specific long-run equilibrium relationships and short-run dynamics using the QARDL framework. This approach extends classical ARDL cointegration analysis to the quantile regression setting, allowing researchers to examine how relationships vary across different points of the conditional distribution.

## Installation

```r
# Install from CRAN (when available)
install.packages("qardlr")
```

## Key Features

- **Quantile regression** across multiple tau values
- **BIC-based automatic lag selection** (p, q)
- **Error Correction Model (ECM)** parameterization
- **Long-run (β), short-run AR (φ), and impact (γ)** parameters
- **Wald tests** for parameter constancy across quantiles
- **Rolling/recursive QARDL** estimation for stability analysis
- **Monte Carlo simulation** for finite-sample properties
- **Publication-ready output** tables (text, LaTeX, HTML)

## Quick Start

```r
library(qardlr)

# Load example data
data(qardl_sim)

# Basic QARDL estimation with automatic lag selection
fit <- qardl(y ~ x1 + x2, data = qardl_sim, 
             tau = c(0.25, 0.50, 0.75))

# View results
summary(fit)

# Wald tests for parameter constancy
wald_results <- qardl_wald(fit)
print(wald_results)

# Generate publication-ready table
cat(qardl_table(fit, type = "latex"))
```

## The QARDL Model

The QARDL(p,q) model is specified as:

$$Q_{y_t}(\tau | \mathcal{F}_{t-1}) = c(\tau) + \sum_{i=1}^{p} \phi_i(\tau) y_{t-i} + \sum_{j=0}^{q-1} \gamma'_j(\tau) x_{t-j}$$

**Key Parameters:**
- **β(τ)**: Long-run parameters = Σγ(τ) / (1 - Σφ(τ))
- **φ(τ)**: Short-run AR coefficients
- **γ(τ)**: Short-run impact parameters
- **ρ(τ)**: Speed of adjustment (ECM) = Σφ(τ) - 1

## Functions

| Function | Description |
|----------|-------------|
| `qardl()` | Main QARDL estimation |
| `qardl_wald()` | Wald tests for parameter constancy |
| `qardl_rolling()` | Rolling/recursive window estimation |
| `qardl_simulate()` | Monte Carlo simulation |
| `qardl_bic_select()` | BIC-based lag selection |
| `qardl_table()` | Publication-ready tables |

## ECM Parameterization

Use `ecm = TRUE` for Error Correction Model form:

```r
fit_ecm <- qardl(y ~ x1 + x2, data = qardl_sim, 
                 tau = c(0.25, 0.50, 0.75), 
                 ecm = TRUE)
summary(fit_ecm)
```

## Rolling Window Analysis

```r
# Rolling QARDL with 50-observation window
roll <- qardl_rolling(y ~ x1 + x2, data = qardl_sim,
                      tau = c(0.25, 0.50, 0.75), 
                      p = 2, q = 2, window = 50)
plot(roll, which = "beta", var = 1)
```

## Monte Carlo Simulation

```r
# Assess finite-sample properties
mc <- qardl_simulate(nobs = 200, reps = 1000, 
                     tau = c(0.25, 0.50, 0.75),
                     p = 1, q = 1, k = 1)
print(mc)
```

## Citation

If you use this package, please cite:

> Cho, J.S., Kim, T.-H., & Shin, Y. (2015). Quantile cointegration in the autoregressive distributed-lag modeling framework. *Journal of Econometrics*, 188(1), 281-300. https://doi.org/10.1016/j.jeconom.2015.01.003

## Author

## License

GPL-3
