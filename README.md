# nimbleNeonorm

<!-- badges: start -->
<!-- badges: end -->

**nimbleNeonorm** makes the MSNBurr and MSNBurr-IIa neonormal distributions
available in [NIMBLE](https://r-nimble.org) models, so that data with
skewness can be modelled with Bayesian inference without writing custom
distributions by hand.

Both distributions have a location `mu`, a scale `sigma` and a shape `alpha`.
Their mode is `mu`, where the density equals that of a normal distribution
with standard deviation `sigma`; `alpha = 1` gives a symmetric distribution,
and other values skew it to one side (Choir, 2020; Ramadani et al., 2025).

The package provides

* density, distribution, quantile and random-generation functions
  (`dmsnburr()`, `pmsnburr()`, `qmsnburr()`, `rmsnburr()` and the `2a`
  versions) that work in R and inside compiled NIMBLE models, with
  derivatives for HMC/NUTS and numerically stable log-scale computations;
* helpers that follow the Bayesian workflow (Gelman et al., 2020): prior and
  posterior predictive simulation, dispersed initial values, convergence
  diagnostics (rank-normalised R-hat and effective sample size) and WAIC.

## Installation

```r
# install.packages("remotes")
remotes::install_github("rafiqzand/nimbleNeonorm")
```

NIMBLE compiles models to C++, so a working compiler is needed (Rtools on
Windows, Xcode command-line tools on macOS).

## Example

```r
library(nimbleNeonorm)
register_neonorm()

code <- nimbleCode({
  for (i in 1:N) y[i] ~ dmsnburr2a(mu, sigma, alpha)
  mu ~ dnorm(0, sd = 5)
  sigma ~ dexp(1)
  alpha ~ dlnorm(0, sdlog = 1)
})

set.seed(1)
y <- replicate(200, rmsnburr2a(1, mu = 0, sigma = 1, alpha = 0.3))

model <- nimbleModel(
  code,
  constants = list(N = length(y)),
  data = list(y = y),
  inits = list(mu = 0, sigma = 1, alpha = 1)
)

fit <- runmcmc_neonorm(model, niter = 4000, sampler = "slice", WAIC = TRUE)
fit                                   # summary with R-hat, ESS and WAIC
yrep <- posterior_predict_neonorm(fit, ndraws = 200)

cleanup_neonorm()
```

See `vignette("bayesian-workflow", package = "nimbleNeonorm")` for the full
workflow: prior predictive check, parameter recovery, convergence
diagnostics, posterior predictive check and model comparison.

## Citation

```r
citation("nimbleNeonorm")
```

## References

Choir, A. S. (2020). *Distribusi neo-normal baru dan karakteristiknya*
[Doctoral dissertation]. Institut Teknologi Sepuluh Nopember.

de Valpine, P., Turek, D., Paciorek, C. J., Anderson-Bergman, C., Temple
Lang, D., & Bodik, R. (2017). Programming with models: Writing statistical
algorithms for general model structures with NIMBLE. *Journal of
Computational and Graphical Statistics*, 26(2), 403–413.
https://doi.org/10.1080/10618600.2016.1172487

Gelman, A., Vehtari, A., Simpson, D., Margossian, C. C., Carpenter, B.,
Yao, Y., Kennedy, L., Gabry, J., Bürkner, P.-C., & Modrák, M. (2020).
*Bayesian workflow*. arXiv:2011.01808.
https://doi.org/10.48550/arXiv.2011.01808

Ramadani, E. P., Choir, A. S., Pravitasari, A. A., & Paraguison, J. (2025).
Adding MSNBURR-IIa distribution to MultiBUGS. *Jurnal Aplikasi Statistika &
Komputasi Statistik*, 17(2), 109–130.
https://doi.org/10.34123/jurnalasks.v17i2.804
