# qardlr 1.1.1

* Removed an empty item from the author list of `qardlr-package.Rd`, which caused an HTML validation NOTE ("trimming empty <li>"). No changes to code.

# qardlr 1.1.0

* Standard errors and Wald tests now follow Cho, Kim and Shin (2015): the covariance of the long-run parameters is Theorem 2 (with (X'PX)^-1 and the error density estimated by a Gaussian kernel with the Bofinger bandwidth), the covariance of phi and gamma is Theorem 1 with the estimator of equation (13), and the Wald tests across quantiles use the joint covariances of Theorems 3 and 4. The previous version assumed independence across quantiles, computed the long-run standard errors by a delta method without covariance terms, and fell back to a fixed covariance of 0.01 when the quantile regression covariance failed.
* Fixed the contemporaneous coefficients with several covariates and q > 1: the coefficient of x1 at lag 1 was reported as the coefficient of x2 at lag 0.
* `gamma` is now the level coefficient gamma(tau) of Cho, Kim and Shin (2015), the sum of the coefficients on x_t, ..., x_{t-q+1}, so that beta = gamma / (1 - sum(phi)) as printed; the coefficients on x_t are in `gamma0`. Constancy tests for gamma are carried out covariate by covariate because its covariance has rank one at each quantile.
* `qardl_simulate()`: the true long-run parameter is now gamma_true / (1 - sum(phi_true)), the value implied by the simulated process (the default beta_true of 1 was not the parameter of the DGP).
* New outputs `cov_joint` and `fhat`.

