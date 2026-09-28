# Shared helpers for the test suite. testthat sources helper-*.R files before
# running the tests.

# Evaluate a scalar nimbleFunction `f` at every element of `x`.
vec <- function(f, x, ...) {
  vapply(x, function(xi) f(xi, ...), numeric(1))
}

# The two distribution families, so that property tests can loop over both.
families <- list(
  MSNBurr = list(
    d = dmsnburr, p = pmsnburr, q = qmsnburr, r = rmsnburr
  ),
  `MSNBurr-IIa` = list(
    d = dmsnburr2a, p = pmsnburr2a, q = qmsnburr2a, r = rmsnburr2a
  )
)

# Names of all user-defined distributions currently known to NIMBLE. The
# tryCatch() covers NIMBLE versions that error when none are registered.
user_distributions <- function() {
  tryCatch(
    nimble:::getAllDistributionsInfo(
      "namesVector",
      nimbleOnly = FALSE,
      userOnly = TRUE
    ),
    error = function(e) character(0)
  )
}

# Register distributions for the duration of the calling test only.
local_neonorm <- function(distributions = "all", env = parent.frame()) {
  cleanup_neonorm(verbose = FALSE)
  register_neonorm(distributions, verbose = FALSE)
  withr::defer(cleanup_neonorm(verbose = FALSE), envir = env)
}

# A small MSNBurr model with priors on all three parameters. Extra arguments
# (e.g. `buildDerivs = TRUE`) are passed to nimble::nimbleModel().
new_msnburr_model <- function(y = 0, ...) {
  code <- nimble::nimbleCode({
    for (i in 1:N) {
      y[i] ~ dmsnburr(mu, sigma, alpha)
    }
    mu ~ dnorm(0, sd = 10)
    sigma ~ dexp(1)
    alpha ~ dgamma(2, 1)
  })

  nimble::nimbleModel(
    code,
    constants = list(N = length(y)),
    data = list(y = y),
    inits = list(mu = 0, sigma = 1, alpha = 1),
    ...
  )
}
