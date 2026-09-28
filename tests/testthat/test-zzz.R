test_that("the startup message points users to register_neonorm()", {
  expect_message(
    nimbleNeonorm:::.onAttach(libname = "", pkgname = "nimbleNeonorm"),
    "register_neonorm()",
    fixed = TRUE
  )
})

test_that("attaching the package does not register any distribution", {
  # Everything registered by earlier tests has been cleaned up, so the
  # package's own bookkeeping should be empty.
  expect_length(nimbleNeonorm:::.nimble_state$registered, 0)
})
