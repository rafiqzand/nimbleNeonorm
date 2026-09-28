# Properties shared by MSNBurr and MSNBurr-IIa. Each test runs once per
# family; see `families` in helper-neonorm.R.

for (family in names(families)) {
  d <- families[[family]]$d
  p <- families[[family]]$p
  q <- families[[family]]$q
  r <- families[[family]]$r


  # ---- Parameter validation --------------------------------------------------

  test_that(paste(family, "density rejects invalid parameters"), {
    expect_error(d(0, 0, 0, 1), "sigma must be > 0")
    expect_error(d(0, 0, -1, 1), "sigma must be > 0")
    expect_error(d(0, 0, 1, 0), "alpha must be > 0")
    expect_error(d(0, 0, 1, -1), "alpha must be > 0")
    expect_error(d(0, Inf, 1, 1), "mu must be finite")
    expect_error(d(0, -Inf, 1, 1), "mu must be finite")
    expect_error(d(0, 0, Inf, 1), "sigma must be finite")
    expect_error(d(0, 0, -Inf, 1), "sigma must be finite")
    expect_error(d(0, 0, 1, Inf), "alpha must be finite")
    expect_error(d(0, 0, 1, -Inf), "alpha must be finite")
  })

  test_that(paste(family, "quantile rejects invalid probabilities"), {
    expect_error(q(-0.1, 0, 1, 1), "p must be between 0 and 1")
    expect_error(q(1.1, 0, 1, 1), "p must be between 0 and 1")
    expect_error(q(0.1, 0, 1, 1, log.p = 1), "log\\(p\\) must be <= 0")
  })

  test_that(paste(family, "random generator validates n and parameters"), {
    expect_error(r(0, 0, 1, 1), "only supports n = 1")
    expect_error(r(2, 0, 1, 1), "only supports n = 1")
    expect_error(r(1, 0, 0, 1), "sigma must be > 0")
    expect_error(r(1, 0, 1, 0), "alpha must be > 0")
    expect_error(r(1, Inf, 1, 1), "mu must be finite")
  })


  # ---- Density ---------------------------------------------------------------

  test_that(paste(family, "density is finite and positive"), {
    y <- vec(d, c(-5, -2, -1, 0, 1, 2, 5), mu = 0, sigma = 1, alpha = 1)
    expect_true(all(is.finite(y)))
    expect_true(all(y > 0))
  })

  test_that(paste(family, "density and log-density agree"), {
    x <- c(-10, -3, -1, 0, 1, 3, 10)
    dens <- vec(d, x, mu = 0.5, sigma = 2, alpha = 1.5, log = 0)
    log_dens <- vec(d, x, mu = 0.5, sigma = 2, alpha = 1.5, log = 1)
    expect_equal(log(dens), log_dens, tolerance = 1e-12)
  })

  test_that(paste(family, "density is zero at +/-Inf"), {
    expect_equal(d(-Inf, 0, 1, 1, log = 0), 0)
    expect_equal(d(Inf, 0, 1, 1, log = 0), 0)
    expect_equal(d(-Inf, 0, 1, 1, log = 1), -Inf)
    expect_equal(d(Inf, 0, 1, 1, log = 1), -Inf)
  })

  test_that(paste(family, "density at the mode is 1 / (sigma sqrt(2 pi))"), {
    for (alpha in c(0.1, 1, 5)) {
      expect_equal(d(1.5, 1.5, 2, alpha), dnorm(0, sd = 2), tolerance = 1e-12)
    }
  })

  test_that(paste(family, "density integrates to one"), {
    total <- integrate(
      function(x) vec(d, x, mu = 0, sigma = 1, alpha = 1.5),
      -Inf, Inf
    )$value
    expect_lt(abs(total - 1), 1e-5)
  })

  test_that(paste(family, "log-density is never NaN for extreme x"), {
    x <- c(-1e6, -1e4, -1000, -100, -50, 50, 100, 1000, 1e4, 1e6)
    for (alpha in c(0.01, 0.1, 1, 10, 100)) {
      log_dens <- vec(d, x, mu = 0, sigma = 1, alpha = alpha, log = 1)
      expect_false(any(is.nan(log_dens)))
    }
  })

  test_that(paste(family, "density obeys the location-scale property"), {
    x <- c(-5, -2, -1, 0, 1, 2, 5)
    mu <- 2.3
    sigma <- 1.7
    alpha <- 1.4
    expect_equal(
      vec(d, x, mu, sigma, alpha),
      vec(d, (x - mu) / sigma, 0, 1, alpha) / sigma,
      tolerance = 1e-12
    )
  })


  # ---- Distribution function -------------------------------------------------

  test_that(paste(family, "CDF lies in [0, 1] and is non-decreasing"), {
    x <- seq(-10, 10, length.out = 101)
    probs <- vec(p, x, mu = 0, sigma = 1, alpha = 1)
    expect_true(all(is.finite(probs)))
    expect_true(all(probs >= 0 & probs <= 1))
    expect_true(all(diff(probs) >= 0))
  })

  test_that(paste(family, "lower and upper tails sum to one"), {
    x <- c(-8, -4, -2, -1, 0, 1, 2, 4, 8)
    lower <- vec(p, x, 0, 1, 1.5, lower.tail = 1)
    upper <- vec(p, x, 0, 1, 1.5, lower.tail = 0)
    expect_equal(lower + upper, rep(1, length(x)), tolerance = 1e-12)
  })

  test_that(paste(family, "log.p = 1 gives the log of the probability"), {
    x <- c(-8, -4, -2, -1, 0, 1, 2, 4, 8)
    for (tail in c(1, 0)) {
      expect_equal(
        log(vec(p, x, 0, 1, 1.5, lower.tail = tail, log.p = 0)),
        vec(p, x, 0, 1, 1.5, lower.tail = tail, log.p = 1),
        tolerance = 1e-12
      )
    }
  })

  test_that(paste(family, "CDF has the right limits at +/-Inf"), {
    expect_equal(p(-Inf, 0, 1, 1, lower.tail = 1, log.p = 0), 0)
    expect_equal(p(Inf, 0, 1, 1, lower.tail = 1, log.p = 0), 1)
    expect_equal(p(-Inf, 0, 1, 1, lower.tail = 0, log.p = 0), 1)
    expect_equal(p(Inf, 0, 1, 1, lower.tail = 0, log.p = 0), 0)
    expect_equal(p(-Inf, 0, 1, 1, lower.tail = 1, log.p = 1), -Inf)
    expect_equal(p(Inf, 0, 1, 1, lower.tail = 1, log.p = 1), 0)
  })

  test_that(paste(family, "CDF derivative equals the density"), {
    x <- c(-5, -2, -1, 0, 1, 2, 5)
    h <- 1e-6
    slope <- (vec(p, x + h, 0, 1, 1.5) - vec(p, x - h, 0, 1, 1.5)) / (2 * h)
    expect_equal(slope, vec(d, x, 0, 1, 1.5), tolerance = 1e-7)
  })

  test_that(paste(family, "CDF obeys the location-scale property"), {
    x <- c(-5, -2, -1, 0, 1, 2, 5)
    mu <- -1.2
    sigma <- 2.4
    alpha <- 1.7
    expect_equal(
      vec(p, x, mu, sigma, alpha),
      vec(p, (x - mu) / sigma, 0, 1, alpha),
      tolerance = 1e-12
    )
  })


  # ---- Quantile function -----------------------------------------------------

  test_that(paste(family, "quantile has the right limits"), {
    expect_equal(q(0, 0, 1, 1, lower.tail = 1), -Inf)
    expect_equal(q(1, 0, 1, 1, lower.tail = 1), Inf)
    expect_equal(q(0, 0, 1, 1, lower.tail = 0), Inf)
    expect_equal(q(1, 0, 1, 1, lower.tail = 0), -Inf)
    expect_equal(q(-Inf, 0, 1, 1, lower.tail = 1, log.p = 1), -Inf)
    expect_equal(q(0, 0, 1, 1, lower.tail = 1, log.p = 1), Inf)
    expect_equal(q(-Inf, 0, 1, 1, lower.tail = 0, log.p = 1), Inf)
    expect_equal(q(0, 0, 1, 1, lower.tail = 0, log.p = 1), -Inf)
  })

  test_that(paste(family, "quantile inverts the CDF"), {
    x <- c(-5, -2, -1, -0.5, 0, 0.5, 1, 2, 5)
    probs <- vec(p, x, 0.3, 1.7, 1.4)
    expect_equal(vec(q, probs, 0.3, 1.7, 1.4), x, tolerance = 1e-10)
  })

  test_that(paste(family, "CDF inverts the quantile"), {
    probs <- c(1e-6, 1e-4, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75,
               0.9, 0.99, 0.999, 0.9999)
    x <- vec(q, probs, -0.5, 2, 2)
    expect_equal(vec(p, x, -0.5, 2, 2), probs, tolerance = 1e-10)
  })

  test_that(paste(family, "quantile agrees across tails and scales"), {
    probs <- c(0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999)
    expect_equal(
      vec(q, 1 - probs, 0, 1.5, 1.8, lower.tail = 0),
      vec(q, probs, 0, 1.5, 1.8, lower.tail = 1),
      tolerance = 1e-10
    )

    probs <- c(1e-6, 1e-4, 0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)
    expect_equal(
      vec(q, log(probs), 0, 1.2, 1.5, log.p = 1),
      vec(q, probs, 0, 1.2, 1.5, log.p = 0),
      tolerance = 1e-10
    )
  })

  test_that(paste(family, "quantile is increasing and finite in the tails"), {
    probs <- c(1e-12, 1e-10, 1e-8, 1e-6, seq(0.001, 0.999, length.out = 101),
               1 - 1e-8, 1 - 1e-10, 1 - 1e-12)
    x <- vec(q, probs, 0, 1, 1.5)
    expect_true(all(is.finite(x)))
    expect_true(all(diff(x) > 0))
  })

  test_that(paste(family, "quantile obeys the location-scale property"), {
    probs <- c(0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999)
    mu <- 3.2
    sigma <- 1.8
    alpha <- 1.6
    expect_equal(
      vec(q, probs, mu, sigma, alpha),
      mu + sigma * vec(q, probs, 0, 1, alpha),
      tolerance = 1e-10
    )
  })


  # ---- Random generation -----------------------------------------------------

  test_that(paste(family, "random draws are finite"), {
    set.seed(123)
    x <- replicate(1000, r(1, 0, 1, 1))
    expect_length(x, 1000)
    expect_true(all(is.finite(x)))
  })

  test_that(paste(family, "random draws follow the CDF"), {
    set.seed(12345)
    x <- replicate(10000, r(1, 0.5, 1.5, 1.4))
    at <- c(-2, 0, 0.5, 1, 3)
    empirical <- vapply(at, function(a) mean(x <= a), numeric(1))
    expect_equal(empirical, vec(p, at, 0.5, 1.5, 1.4), tolerance = 0.03)
  })
}


# ---- Relation between the two families ---------------------------------------

test_that("MSNBurr-IIa is the mirror image of MSNBurr about mu", {
  x <- c(-4, -1, 0, 0.7, 2, 6)
  mu <- 0.7
  expect_equal(
    vec(dmsnburr2a, x, mu, 1.3, 2.5),
    vec(dmsnburr, 2 * mu - x, mu, 1.3, 2.5),
    tolerance = 1e-12
  )
  expect_equal(
    vec(pmsnburr2a, x, mu, 1.3, 2.5),
    vec(pmsnburr, 2 * mu - x, mu, 1.3, 2.5, lower.tail = 0),
    tolerance = 1e-12
  )
})
