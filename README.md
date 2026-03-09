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

**Dr. Merwan Roudane**  
Email: merwanroudane920@gmail.com

## License

GPL-3
