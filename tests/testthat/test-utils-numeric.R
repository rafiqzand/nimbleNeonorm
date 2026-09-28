# ---- log1pexp_nimble --------------------------------------------------------

test_that("log1pexp_nimble() matches log1p(exp(x)) in every region", {
  x <- c(-100, -50, -37, -20, -10, -1, -1e-8, 0, 1e-8, 1,
         10, 18, 20, 30, 33.3, 40, 100)
  expect_equal(vec(log1pexp_nimble, x), log1p(exp(x)), tolerance = 1e-12)
})

test_that("log1pexp_nimble() stays finite for extreme inputs", {
  x <- c(-1000, -100, -50, -37, 0, 18, 33.3, 100, 1000)
  expect_true(all(is.finite(vec(log1pexp_nimble, x))))
})

test_that("log1pexp_nimble() is monotonically increasing across cut-offs", {
  x <- seq(-100, 100, length.out = 1001)
  expect_true(all(diff(vec(log1pexp_nimble, x)) >= 0))
})


# ---- expm1_nimble -----------------------------------------------------------

test_that("expm1_nimble() matches expm1() on both sides of the cut-off", {
  x <- c(-0.7, -0.5, -0.1, -1e-8, -1e-12, 0, 1e-12, 1e-8,
         0.1, 0.5, 0.7, 1, 5, 10)
  expect_equal(vec(expm1_nimble, x), expm1(x), tolerance = 1e-12)
})

test_that("expm1_nimble() is zero at zero", {
  expect_equal(expm1_nimble(0), 0)
})

test_that("expm1_nimble() stays finite for finite inputs", {
  x <- c(-10, -1, -0.7, 0, 0.7, 1, 10, 100)
  expect_true(all(is.finite(vec(expm1_nimble, x))))
})


# ---- log1mexp_nimble --------------------------------------------------------

test_that("log1mexp_nimble() agrees with a 1024-bit Rmpfr reference", {
  skip_if_not_installed("Rmpfr")

  x <- c(1e-20, 1e-15, 1e-12, 1e-10, 1e-8, 1e-6, 1e-3, 0.1,
         log(2), 1, 2, 10, 100)
  expected <- vapply(x, function(a) {
    a <- Rmpfr::mpfr(a, precBits = 1024)
    as.numeric(log(-expm1(-a)))
  }, numeric(1))

  expect_equal(vec(log1mexp_nimble, x), expected, tolerance = 1e-14)
})

test_that("log1mexp_nimble() returns -Inf at zero", {
  expect_equal(log1mexp_nimble(0), -Inf)
})

test_that("log1mexp_nimble() rejects negative input", {
  expect_error(log1mexp_nimble(-1), "x must be >= 0")
})


# ---- log_omega_msnburr ------------------------------------------------------

test_that("log_omega_msnburr() matches the closed form", {
  alpha <- c(0.01, 0.5, 1, 2, 10, 1e6)
  expected <- (alpha + 1) * log1p(1 / alpha) - 0.5 * log(2 * pi)
  expect_equal(vec(log_omega_msnburr, alpha), expected, tolerance = 1e-12)
})
