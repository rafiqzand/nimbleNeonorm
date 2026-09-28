test_that("register_msnburr() registers dmsnburr with the right metadata", {
  local_neonorm("msnburr")

  expect_true("dmsnburr" %in% user_distributions())
  info <- nimble::getDistributionInfo("dmsnburr")
  expect_equal(info$densityName, "dmsnburr")
  expect_equal(info$simulateName, "rmsnburr")
  expect_true(isTRUE(info$pqAvail))
  expect_false(isTRUE(info$discrete))
})

test_that("register_msnburr2a() registers dmsnburr2a with the right metadata", {
  local_neonorm("msnburr2a")

  expect_true("dmsnburr2a" %in% user_distributions())
  info <- nimble::getDistributionInfo("dmsnburr2a")
  expect_equal(info$densityName, "dmsnburr2a")
  expect_equal(info$simulateName, "rmsnburr2a")
  expect_true(isTRUE(info$pqAvail))
  expect_false(isTRUE(info$discrete))
})

test_that("register_neonorm('all') registers both and returns their names", {
  cleanup_neonorm(verbose = FALSE)
  withr::defer(cleanup_neonorm(verbose = FALSE))

  registered <- register_neonorm("all", verbose = FALSE)

  expect_equal(registered, c("MSNBurr", "MSNBurr-IIa"))
  expect_true(all(c("dmsnburr", "dmsnburr2a") %in% user_distributions()))
})

test_that("distribution names are matched case-insensitively", {
  cleanup_neonorm(verbose = FALSE)
  withr::defer(cleanup_neonorm(verbose = FALSE))

  expect_equal(register_neonorm("MSNBurr-IIa", verbose = FALSE), "MSNBurr-IIa")
  expect_equal(register_neonorm("MSNBURR", verbose = FALSE), "MSNBurr")
})

test_that("registering twice does not fail", {
  local_neonorm("msnburr")
  expect_no_error(register_msnburr(verbose = FALSE))
})

test_that("unknown distribution names are rejected", {
  expect_error(
    register_neonorm("invalid_distribution", verbose = FALSE),
    "Unknown distribution"
  )
})

test_that("cleanup_neonorm() removes the distributions from NIMBLE", {
  cleanup_neonorm(verbose = FALSE)
  register_neonorm("all", verbose = FALSE)

  cleanup_neonorm(verbose = FALSE)

  expect_false(any(c("dmsnburr", "dmsnburr2a") %in% user_distributions()))
  expect_length(nimbleNeonorm:::.nimble_state$registered, 0)
})

test_that("cleanup_neonorm() is a no-op when nothing is registered", {
  cleanup_neonorm(verbose = FALSE)
  expect_message(cleanup_neonorm(), "No neo-normal distributions")
})

test_that("a model using dmsnburr can be built and evaluated", {
  local_neonorm("msnburr")

  code <- nimble::nimbleCode({
    y ~ dmsnburr(mu, sigma, alpha)
  })
  model <- nimble::nimbleModel(
    code,
    constants = list(mu = 0, sigma = 1, alpha = 1),
    data = list(y = 0)
  )

  # At the mode the density is dnorm(0), whatever alpha is.
  expect_equal(model$calculate(), dnorm(0, log = TRUE), tolerance = 1e-10)
})
