# ---- configure_mcmc_neonorm -------------------------------------------------

test_that("configure_mcmc_neonorm() works with the basic samplers", {
  local_neonorm("msnburr")
  model <- new_msnburr_model(y = rep(0, 5))

  for (sampler in c("default", "slice", "RW")) {
    conf <- configure_mcmc_neonorm(model, sampler = sampler)
    expect_s4_class(conf, "MCMCconf")
  }
})

test_that("the slice and RW strategies put the requested sampler on alpha", {
  local_neonorm("msnburr")
  model <- new_msnburr_model()

  for (sampler in c("slice", "RW")) {
    conf <- configure_mcmc_neonorm(model, sampler = sampler)
    alpha_samplers <- conf$getSamplers("alpha")
    expect_length(alpha_samplers, 1)
    expect_equal(alpha_samplers[[1]]$name, sampler)
  }
})

test_that("a model without an alpha node gives a warning", {
  code <- nimble::nimbleCode({
    y ~ dnorm(mu, sd = 1)
    mu ~ dnorm(0, sd = 10)
  })
  model <- nimble::nimbleModel(code, data = list(y = 0), inits = list(mu = 0))

  expect_warning(configure_mcmc_neonorm(model, sampler = "slice"), "alpha")
})

test_that("unknown samplers and non-models are rejected", {
  local_neonorm("msnburr")
  model <- new_msnburr_model()

  expect_error(
    configure_mcmc_neonorm(model, sampler = "not_a_sampler"),
    "Unknown sampler"
  )
  expect_error(configure_mcmc_neonorm(list()), "NIMBLE model")
})

test_that("NUTS-based samplers ask for nimbleHMC when it is missing", {
  skip_if(
    requireNamespace("nimbleHMC", quietly = TRUE),
    "nimbleHMC is installed"
  )
  local_neonorm("msnburr")
  model <- new_msnburr_model(buildDerivs = TRUE)

  expect_error(configure_mcmc_neonorm(model, sampler = "NUTS"), "nimbleHMC")
  expect_error(
    configure_mcmc_neonorm(model, sampler = "mixed_nuts"),
    "nimbleHMC"
  )
})


# ---- runmcmc_neonorm: argument checks (no compilation) ----------------------

test_that("runmcmc_neonorm() validates its arguments before compiling", {
  local_neonorm("msnburr")
  model <- new_msnburr_model()

  expect_error(runmcmc_neonorm(model, niter = 10, nburnin = 10), "nburnin")
  expect_error(runmcmc_neonorm(model, niter = 0), "niter")
  expect_error(runmcmc_neonorm(model, niter = 10, thin = 0), "thin")
  expect_error(runmcmc_neonorm(model, niter = 10, nchains = 1.5), "nchains")
  expect_error(runmcmc_neonorm(model, niter = 10, setSeed = "a"), "setSeed")
})


# ---- runmcmc_neonorm: full runs (compile C++, so slow) ----------------------

test_that("runmcmc_neonorm() returns draws, config and summary", {
  skip_on_cran()
  local_neonorm("msnburr")
  model <- new_msnburr_model()

  fit <- runmcmc_neonorm(
    model,
    niter = 30,
    nburnin = 10,
    sampler = "slice",
    setSeed = 123
  )

  expect_named(fit, c("samples", "config", "summary"))
  expect_true(is.matrix(fit$samples))
  expect_equal(nrow(fit$samples), 20)
  expect_setequal(colnames(fit$samples), c("mu", "sigma", "alpha"))
  expect_equal(colnames(fit$summary), c("Mean", "SD", "2.5%", "50%", "97.5%"))
})

test_that("runmcmc_neonorm() discards burn-in before thinning", {
  skip_on_cran()
  local_neonorm("msnburr")
  model <- new_msnburr_model()

  fit <- runmcmc_neonorm(
    model,
    niter = 30,
    nburnin = 10,
    thin = 2,
    summary = FALSE
  )

  expect_equal(nrow(fit$samples), 10)
  expect_null(fit$summary)
})
