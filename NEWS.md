# nimbleNeonorm (development version)

## Documentation

* Added real-data case studies in `analysis/stocks-case-study/`: the Bayesian
  workflow applied to the daily returns of ADRO and PTBA (2010-2019). These
  scripts are not part of the installed package.
* `vignette("bayesian-workflow")` is now described as a simulated example, in
  which the true parameter values are known.

# nimbleNeonorm 0.2.0

## Breaking changes

* Invalid inputs to the d, p, q and r functions (non-finite or `NaN`
  parameters, `sigma <= 0`, `alpha <= 0`, probabilities outside `[0, 1]`)
  now return `NaN`, as R's and NIMBLE's built-in distributions do, instead
  of stopping with an error. `n != 1` in the r functions is still an error.
* `runmcmc_neonorm()` now runs 4 chains by default (was 1) and returns an
  object of class `neonorm_fit`.

## New features

* `prior_inits_neonorm()` draws dispersed initial values from the prior.
  `runmcmc_neonorm()` uses it when `inits = NULL` and `nchains > 1`, so that
  R-hat can detect chains that have not mixed.
* `summarise_neonorm()` summarises draws with the rank-normalised R-hat and
  bulk and tail effective sample sizes (via the posterior package), and
  warns when R-hat > 1.01 or ESS < 100 per chain. `runmcmc_neonorm()` uses
  it for its `summary`.
* `runmcmc_neonorm(WAIC = TRUE)` and `configure_mcmc_neonorm(enableWAIC =
  TRUE)` compute WAIC for model comparison.
* `prior_predict_neonorm()` and `posterior_predict_neonorm()` simulate
  prior and posterior predictive data.
* `runmcmc_neonorm()` checks that the log-probability at the starting
  values is finite before compiling.
* New vignette, `vignette("bayesian-workflow")`.

## Documentation

* Every exported function documents the references its methods rely on.
* Added `inst/CITATION`.

# nimbleNeonorm 0.1.1

* MSNBurr and MSNBurr-IIa d, p, q and r nimbleFunctions with numerically
  stable log-scale computations, registration helpers for NIMBLE, and MCMC
  helpers.
