# [REVISI-06] FILE TEST BARU untuk R/predictive.R.

# ---- prior_predict_neonorm (uncompiled, no C++) -----------------------------

test_that("prior_predict_neonorm() returns parameter and data draws", {
  local_neonorm("msnburr")
  model <- new_msnburr_model(y = c(-0.5, 0, 0.8))

  prior <- prior_predict_neonorm(model, nsim = 20, seed = 1)

  expect_named(prior, c("parameters", "data"))
  expect_equal(dim(prior$parameters), c(20, 3))
  expect_setequal(colnames(prior$parameters), c("mu", "sigma", "alpha"))
  expect_equal(dim(prior$data), c(20, 3))
  expect_equal(colnames(prior$data), c("y[1]", "y[2]", "y[3]"))
  expect_true(all(is.finite(prior$data)))
  expect_true(all(prior$parameters[, "sigma"] > 0))
})

test_that("prior_predict_neonorm() restores the model and is reproducible", {
  local_neonorm("msnburr")
  y <- c(-0.5, 0, 0.8)
  model <- new_msnburr_model(y = y)
  lp_before <- model$calculate()

  first <- prior_predict_neonorm(model, nsim = 5, seed = 7)
  second <- prior_predict_neonorm(model, nsim = 5, seed = 7)

  expect_identical(first, second)
  expect_equal(as.numeric(model$y), y)
  expect_equal(model$calculate(), lp_before)
})

test_that("prior_predict_neonorm() validates its arguments", {
  local_neonorm("msnburr")
  model <- new_msnburr_model()
  expect_error(prior_predict_neonorm(model, nsim = 0), "nsim")
  expect_error(prior_predict_neonorm(model, seed = "a"), "seed")
  expect_error(prior_predict_neonorm(list()), "NIMBLE model")
})


# ---- posterior_predict_neonorm (compiles C++, so slow) ----------------------

test_that("posterior_predict_neonorm() simulates replicated data", {
  skip_on_cran()
  local_neonorm("msnburr")
  y <- c(-1.2, -0.3, 0.1, 0.4, 0.9, 1.6)
  model <- new_msnburr_model(y = y)

  fit <- suppressWarnings(runmcmc_neonorm(
    model, niter = 200, nchains = 2, sampler = "slice"
  ))
  yrep <- posterior_predict_neonorm(fit, ndraws = 50, seed = 1)

  expect_equal(dim(yrep), c(50, length(y)))
  expect_true(all(is.finite(yrep)))
  # Replicates differ from each other and the observed data are restored
  expect_gt(stats::sd(yrep[, 1]), 0)
  expect_equal(as.numeric(fit$compiled$model$y), y)

  expect_error(posterior_predict_neonorm(fit, ndraws = 1e6), "ndraws")
  expect_error(posterior_predict_neonorm(list()), "runmcmc_neonorm")
})
