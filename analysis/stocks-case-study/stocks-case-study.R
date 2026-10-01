# Case study: daily stock returns of ADRO and PTBA (2010-2019)
# Bayesian workflow with nimbleNeonorm 0.2.0
#
# Selection criteria (following Choir, 2020, Section 5.2):
#   1. stocks that stayed in the LQ45 index throughout 2010-2019 (16 firms);
#   2. no autocorrelation in the returns (Ljung-Box test, lag 10, p > 0.05).
# Return: r_t = ln(P_t) - ln(P_{t-1}), with P the daily closing price.
#
# Data (one CSV per stock and source, columns `date` and `close`):
#   analysis/stocks-case-study/data/google/ADRO.csv   (main source, as in Choir)
#   analysis/stocks-case-study/data/google/PTBA.csv
#   analysis/stocks-case-study/data/yahoo/ADRO.csv    (robustness check)
#   analysis/stocks-case-study/data/yahoo/PTBA.csv
#
# Usage (from the package root, with nimbleNeonorm installed):
#   Rscript analysis/stocks-case-study/stocks-case-study.R google
#   Rscript analysis/stocks-case-study/stocks-case-study.R yahoo
# or, inside R: SOURCE <- "google"; source("analysis/stocks-case-study/stocks-case-study.R")
#
# Results go to analysis/stocks-case-study/output/<source>/<ticker>/.

library(nimbleNeonorm)
register_neonorm(verbose = FALSE)

# ---- Settings ---------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (!exists("SOURCE")) SOURCE <- if (length(args) > 0) args[1] else "google"
stopifnot(SOURCE %in% c("google", "yahoo"))

TICKERS <- c("ADRO", "PTBA")
BASE_DIR <- "analysis/stocks-case-study"

TICKERS <- TICKERS[file.exists(file.path(BASE_DIR, "data", SOURCE,
                                         paste0(TICKERS, ".csv")))]


skew_moment <- function(x) mean((x - mean(x))^3) / sd(x)^3
kurt_excess <- function(x) mean((x - mean(x))^4) / sd(x)^4 - 3
skew_bowley <- function(x) {
  q <- quantile(x, c(0.25, 0.5, 0.75), names = FALSE)
  (q[3] + q[1] - 2 * q[2]) / (q[3] - q[1])
}


# ---- 0. Data ----------------------------------------------------------------

get_returns <- function(ticker, source) {
  csv <- file.path(BASE_DIR, "data", source, paste0(ticker, ".csv"))
  if (!file.exists(csv)) stop("Missing data file: ", csv)

  px <- read.csv(csv)
  names(px) <- tolower(names(px))              # accepts Date/Close headers
  px$date <- as.Date(sub(" .*$", "", px$date),
                     tryFormats = c("%Y-%m-%d", "%m/%d/%Y"))
  if (anyNA(px$date)) stop("Unparseable dates in ", csv)
  px$close <- suppressWarnings(as.numeric(px$close))
  px <- px[order(px$date), ]
  px <- px[!is.na(px$date) & !is.na(px$close) & px$close > 0, ]
  diff(log(px$close))
}


# ---- 2. Models and priors ---------------------------------------------------
# Returns are standardised (mean 0, sd 1), so the weakly informative priors
# below are on the same scale for every stock:
#   mu ~ N(0, 5^2), sigma ~ Exp(1), alpha ~ LogNormal(0, 1) (median 1).

model_code <- list(
  normal = nimbleCode({
    for (i in 1:N) y[i] ~ dnorm(mu, sd = sigma)
    mu ~ dnorm(0, sd = 5)
    sigma ~ dexp(1)
  }),
  msnburr = nimbleCode({
    for (i in 1:N) y[i] ~ dmsnburr(mu, sigma, alpha)
    mu ~ dnorm(0, sd = 5)
    sigma ~ dexp(1)
    alpha ~ dlnorm(0, sdlog = 1)
  }),
  msnburr2a = nimbleCode({
    for (i in 1:N) y[i] ~ dmsnburr2a(mu, sigma, alpha)
    mu ~ dnorm(0, sd = 5)
    sigma ~ dexp(1)
    alpha ~ dlnorm(0, sdlog = 1)
  }),
  # Optional benchmark: Student's t. Choir (2020) found that it beats MSNBurr
  # on stock returns. Not tested with the package; see `include_t` below.
  student_t = nimbleCode({
    for (i in 1:N) y[i] ~ dt_nonstandard(mu, sigma, nu)
    mu ~ dnorm(0, sd = 5)
    sigma ~ dexp(1)
    nu ~ dgamma(2, 0.1)
  })
)

make_model <- function(family, y) {
  inits <- list(mu = 0, sigma = 1)
  if (family %in% c("msnburr", "msnburr2a")) inits$alpha <- 1
  if (family == "student_t") inits$nu <- 10
  nimbleModel(model_code[[family]], constants = list(N = length(y)),
              data = list(y = y), inits = inits)
}


# ---- Main function: the whole workflow for one stock ------------------------

run_case <- function(ticker, source, include_t = FALSE,
                     niter = 4000, nburnin = 2000) {
  out <- file.path(BASE_DIR, "output", source, ticker)
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  cat("\n=====================", ticker, "/", source, "=====================\n")

  r <- get_returns(ticker, source)

  # -- 1. Data characteristics and autocorrelation test ----------------------
  lb <- Box.test(r, lag = 10, type = "Ljung-Box")
  desc <- c(
    n = length(r), mean = mean(r), median = median(r), var = var(r),
    min = min(r), max = max(r),
    q1 = unname(quantile(r, 0.25)), q3 = unname(quantile(r, 0.75)),
    skewness = skew_moment(r), bowley = skew_bowley(r),
    excess_kurtosis = kurt_excess(r), n_zero = sum(r == 0),
    ljung_box_X2 = unname(lb$statistic), ljung_box_p = lb$p.value
  )
  print(signif(desc, 4))
  print(shapiro.test(sample(r, min(length(r), 5000))))
  write.csv(data.frame(statistic = names(desc), value = unname(desc)),
            file.path(out, "desc.csv"), row.names = FALSE)

  png(file.path(out, "01-data-characteristics.png"), 1800, 600, res = 150)
  op <- par(mfrow = c(1, 3))
  plot(r, type = "l", main = paste("Daily returns", ticker),
       xlab = "day", ylab = "return")
  hist(r, breaks = 80, freq = FALSE, main = "Histogram", xlab = "return")
  lines(density(r), lwd = 2)
  qqnorm(r, main = "Normal Q-Q plot"); qqline(r)
  par(op); dev.off()

  # Standardise (mean and sd) so the priors match the scale of the data.
  center <- mean(r)
  scale <- sd(r)
  y <- (r - center) / scale

  # -- 3. Prior predictive check ----------------------------------------------
  prior <- prior_predict_neonorm(make_model("msnburr", y), nsim = 200, seed = 1)
  prior_range <- apply(prior$data, 1, function(d) diff(range(d)))
  cat("\nRange of the prior predictive data:\n"); print(summary(prior_range))
  cat("Range of the observed data (standardised):",
      round(diff(range(y)), 2), "\n")
  png(file.path(out, "02-prior-predictive.png"), 1000, 700, res = 150)
  hist(prior$data[1:20, ], breaks = 100, main = "Prior predictive data",
       xlab = "y (standardised, 20 data sets)")
  dev.off()

  # -- 4. MCMC ----------------------------------------------------------------
  families <- c("normal", "msnburr", "msnburr2a", if (include_t) "student_t")
  set.seed(2026)
  fits <- lapply(families, function(f) {
    runmcmc_neonorm(
      make_model(f, y), niter = niter, nburnin = nburnin, WAIC = TRUE,
      sampler = if (f %in% c("msnburr", "msnburr2a")) "slice" else "default"
    )
  })
  names(fits) <- families
  saveRDS(fits, file.path(out, "fits.rds"))

  # -- 5. Convergence diagnostics ---------------------------------------------
  # Criteria: R-hat < 1.01 and ESS > 100 per chain (Vehtari et al., 2021).
  for (m in families) {
    cat("\n==", ticker, m, "(standardised scale) ==\n")
    print(round(fits[[m]]$summary, 3))
    write.csv(round(fits[[m]]$summary, 4),
              file.path(out, paste0("convergence-", m, ".csv")))
  }

  png(file.path(out, "03-traceplot-msnburr.png"), 1600, 500, res = 150)
  op <- par(mfrow = c(1, 3))
  for (p in c("mu", "sigma", "alpha")) {
    coda::traceplot(fits$msnburr$samples[, p], main = p)
  }
  par(op); dev.off()

  png(file.path(out, "03-traceplot-msnburr2a.png"), 1600, 500, res = 150)
  op <- par(mfrow = c(1, 3))
  for (p in c("mu", "sigma", "alpha")) {
    coda::traceplot(fits$msnburr2a$samples[, p], main = p)
  }
  par(op); dev.off()

  # Back to the original scale for mu and sigma; alpha is unaffected by the
  # standardisation.
  back <- function(fit) {
    d <- do.call(rbind, lapply(fit$samples, as.matrix))
    cbind(mu = center + scale * d[, "mu"], sigma = scale * d[, "sigma"],
          alpha = d[, "alpha"])
  }
  for (m in c("msnburr", "msnburr2a")) {
    tab <- t(apply(back(fits[[m]]), 2, function(x) {
      c(Mean = mean(x), quantile(x, c(0.025, 0.975)))
    }))
    cat("\n", ticker, m, "- original scale:\n"); print(round(tab, 5))
    write.csv(round(tab, 6), file.path(out, paste0("original-scale-", m, ".csv")))
  }

  # -- 6. Posterior predictive check ------------------------------------------
  # Test statistics: skewness, kurtosis, median and IQR. A posterior
  # predictive p-value near 0 or 1 means the model cannot reproduce that
  # statistic (Gelman et al., 2013, Ch. 6).
  stats_fun <- list(skewness = skew_moment, kurtosis = kurt_excess,
                    median = median, IQR = IQR)
  ppc <- list()
  png(file.path(out, "04-posterior-predictive.png"), 600 * length(families),
      600, res = 150)
  op <- par(mfrow = c(1, length(families)))
  for (m in families) {
    yrep <- posterior_predict_neonorm(fits[[m]], ndraws = 500, seed = 1)
    plot(density(y), lwd = 2, main = m, xlab = "y (standardised)",
         xlim = quantile(y, c(0.002, 0.998)))
    for (i in 1:50) lines(density(yrep[i, ]), col = grDevices::gray(0.6, 0.4))
    lines(density(y), lwd = 2)
    ppc[[m]] <- vapply(stats_fun, function(f) {
      mean(apply(yrep, 1, f) >= f(y))
    }, numeric(1))
  }
  par(op); dev.off()
  ppc <- do.call(rbind, ppc)
  cat("\nPosterior predictive p-values:\n"); print(round(ppc, 3))
  write.csv(round(ppc, 4), file.path(out, "ppc.csv"))

  # -- 7. Model comparison (WAIC; smaller is better) --------------------------
  waic <- vapply(fits, function(f) f$WAIC$WAIC, numeric(1))
  cat("\nWAIC:\n"); print(sort(round(waic, 2)))
  cat("Difference from the best model:\n")
  print(round(sort(waic) - min(waic), 2))
  write.csv(data.frame(model = names(waic), WAIC = waic,
                       delta = waic - min(waic)),
            file.path(out, "waic.csv"), row.names = FALSE)

  invisible(list(desc = desc, ppc = ppc, waic = waic))
}


# ---- Run --------------------------------------------------------------------
# Set include_t = TRUE in run_case() to add Student's t as a benchmark.

res <- lapply(setNames(TICKERS, TICKERS), run_case, source = SOURCE)
saveRDS(res, file.path(BASE_DIR, paste0("results-", SOURCE, ".rds")))
cat("\nDone. Results are in", file.path(BASE_DIR, "output", SOURCE), "\n")
