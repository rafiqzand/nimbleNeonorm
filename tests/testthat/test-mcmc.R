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

test_that("enableWAIC is passed to the configuration and then restored", {
  local_neonorm("msnburr")
  model <- new_msnburr_model()
  before <- nimble::nimbleOptions("MCMCenableWAIC")

  conf <- configure_mcmc_neonorm(model, enableWAIC = TRUE)
  expect_true(isTRUE(conf$enableWAIC))
  expect_identical(nimble::nimbleOptions("MCMCenableWAIC"), before)

  expect_error(configure_mcmc_neonorm(model, enableWAIC = NA), "enableWAIC")
})


# ---- prior_inits_neonorm ----------------------------------------------------

test_that("prior_inits_neonorm() gives distinct, valid inits per chain", {
  local_neonorm("msnburr")
  model <- new_msnburr_model(y = c(-0.5, 0, 0.8))
  before <- nimble::values(model, c("mu", "sigma", "alpha"))

  inits <- prior_inits_neonorm(model, nchains = 3, seed = 42)

  expect_length(inits, 3)
  for (ini in inits) {
    expect_setequal(names(ini), c("mu", "sigma", "alpha"))
    expect_gt(ini$sigma, 0)
    expect_gt(ini$alpha, 0)
  }
  expect_false(identical(inits[[1]], inits[[2]]))
  # Reproducible with the same seed
  expect_identical(inits, prior_inits_neonorm(model, nchains = 3, seed = 42))
  # The model is left untouched
  expect_equal(nimble::values(model, c("mu", "sigma", "alpha")), before)
})

test_that("prior_inits_neonorm() does not change the global RNG state", {
  local_neonorm("msnburr")
  model <- new_msnburr_model()
  set.seed(1)
  expected <- runif(1)
  set.seed(1)
  prior_inits_neonorm(model, nchains = 2, seed = 99)
  expect_identical(runif(1), expected)
})


# ---- summarise_neonorm ------------------------------------------------------

test_that("summarise_neonorm() accepts a matrix, a list and an mcmc.list", {
  set.seed(1)
  chains <- lapply(1:4, function(i) cbind(mu = rnorm(400), s = rexp(400)))
  cols <- c("Mean", "SD", "2.5%", "50%", "97.5%",
            "Rhat", "ESS_bulk", "ESS_tail")

  from_list <- summarise_neonorm(chains)
  expect_equal(colnames(from_list), cols)
  expect_equal(rownames(from_list), c("mu", "s"))

  from_coda <- summarise_neonorm(coda::mcmc.list(lapply(chains, coda::mcmc)))
  expect_equal(from_coda, from_list)

  from_matrix <- summarise_neonorm(chains[[1]], warn = FALSE)
  expect_equal(from_matrix["mu", "Mean"], mean(chains[[1]][, "mu"]))
})

test_that("summarise_neonorm() reports well-mixed chains without warning", {
  set.seed(2)
  chains <- lapply(1:4, function(i) cbind(mu = rnorm(1000)))
  expect_no_warning(tab <- summarise_neonorm(chains))
  expect_lt(tab["mu", "Rhat"], 1.01)
  expect_gt(tab["mu", "ESS_bulk"], 400)
})

test_that("summarise_neonorm() warns about chains that disagree", {
  set.seed(3)
  chains <- lapply(1:4, function(i) cbind(mu = rnorm(500, mean = i)))
  expect_warning(summarise_neonorm(chains, ess_per_chain = 0), "R-hat")
})

test_that("summarise_neonorm() rejects malformed input", {
  expect_error(summarise_neonorm(1:10), "matrix")
  expect_error(
    summarise_neonorm(list(cbind(a = 1:5), cbind(b = 1:5))),
    "same named columns"
  )
  expect_error(summarise_neonorm(cbind(a = 1:5), probs = 2), "probs")
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
  expect_error(runmcmc_neonorm(model, niter = 10, WAIC = "yes"), "WAIC")
  expect_error(runmcmc_neonorm(model, niter = 10, summary = NA), "summary")
})

test_that("runmcmc_neonorm() stops early on invalid starting values", {
  local_neonorm("msnburr")
  model <- new_msnburr_model()
  model$sigma <- -1

  expect_error(
    runmcmc_neonorm(model, niter = 10, nchains = 1),
    "log-probability"
  )
})


# ---- runmcmc_neonorm: full runs (compile C++, so slow) ----------------------

test_that("runmcmc_neonorm() returns a neonorm_fit with draws and summary", {
  skip_on_cran()
  local_neonorm("msnburr")
  model <- new_msnburr_model()

  # Too few iterations for reliable diagnostics; only the structure is tested.
  fit <- suppressWarnings(runmcmc_neonorm(
    model,
    niter = 30,
    nburnin = 10,
    nchains = 1,
    sampler = "slice",
    setSeed = 123
  ))

  expect_s3_class(fit, "neonorm_fit")
  expect_true(all(
    c("samples", "summary", "config", "model", "compiled", "settings") %in%
      names(fit)
  ))
  expect_true(is.matrix(fit$samples))
  expect_equal(nrow(fit$samples), 20)
  expect_setequal(colnames(fit$samples), c("mu", "sigma", "alpha"))
  expect_equal(
    colnames(fit$summary),
    c("Mean", "SD", "2.5%", "50%", "97.5%", "Rhat", "ESS_bulk", "ESS_tail")
  )
  expect_output(print(fit), "nimbleNeonorm MCMC fit")
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
    nchains = 1,
    summary = FALSE
  )

  expect_equal(nrow(fit$samples), 10)
  expect_null(fit$summary)
})

test_that("NUTS compiles the derivatives and runs", {
  skip_on_cran()
  skip_if_not_installed("nimbleHMC")
  local_neonorm("msnburr")
  model <- new_msnburr_model(
    y = c(-1.2, -0.3, 0.1, 0.4, 0.9, 1.6),
    buildDerivs = TRUE
  )

  fit <- suppressWarnings(runmcmc_neonorm(
    model,
    niter = 100,
    nchains = 1,
    sampler = "NUTS",
    summary = FALSE
  ))

  expect_equal(nrow(fit$samples), 50)
  expect_true(all(is.finite(fit$samples)))
  expect_true(all(fit$samples[, "sigma"] > 0))
  expect_true(all(fit$samples[, "alpha"] > 0))
})

test_that("several chains start from dispersed inits and give WAIC", {
  skip_on_cran()
  local_neonorm("msnburr")
  model <- new_msnburr_model(y = c(-1.2, -0.3, 0.1, 0.4, 0.9, 1.6))

  fit <- suppressWarnings(runmcmc_neonorm(
    model,
    niter = 200,
    nchains = 2,
    sampler = "slice",
    WAIC = TRUE
  ))

  expect_s3_class(fit$samples, "mcmc.list")
  expect_length(fit$samples, 2)
  expect_length(fit$inits, 2)
  expect_false(identical(fit$inits[[1]], fit$inits[[2]]))
  expect_true(is.finite(fit$WAIC$WAIC))
})
