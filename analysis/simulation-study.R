# Simulation study: parameter recovery of the MSNBurr and MSNBurr-IIa models
# fitted with nimbleNeonorm.
#
# For each distribution family, sample size n and true shape alpha, the
# script simulates N_REP data sets with mu = 0 and sigma = 1, fits the model
# with four chains, and records the posterior mean, the 95% credible
# interval, R-hat and the effective sample sizes. It then summarises bias,
# mean squared error (MSE), coverage of the 95% interval and convergence.
#
# The model and the MCMC are compiled once per (family, n) and reused for
# every data set, so only the first fit of each combination pays for C++
# compilation.
#
# Design (revised after the thesis defence): a single large sample size,
# n = 10000, so that 2 families x 3 shapes = 6 scenarios. Each replication is
# saved to CHECKPOINT_DIR as soon as it finishes, so an interrupted run (for
# example a Google Colab session that disconnects) resumes where it stopped
# when the script is run again. Set PILOT <- TRUE first to time one
# replication per scenario before starting the full run.
#
# References
#   Cook, S. R., Gelman, A., & Rubin, D. B. (2006). Validation of software
#     for Bayesian models using posterior quantiles. Journal of
#     Computational and Graphical Statistics, 15(3), 675-692.
#   Gelman, A., et al. (2020). Bayesian workflow. arXiv:2011.01808.
#   Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
#     studies to evaluate statistical methods. Statistics in Medicine,
#     38(11), 2074-2102.
#   Vehtari, A., et al. (2021). Rank-normalization, folding, and
#     localization: An improved R-hat for assessing convergence of MCMC.
#     Bayesian Analysis, 16(2), 667-718.

library(nimbleNeonorm)
register_neonorm(verbose = FALSE)

# ---- Settings ---------------------------------------------------------------

N_REP <- 50                       # data sets per scenario; 100 lowers the
# Monte Carlo SE of a 0.95 coverage from
# about 0.031 to about 0.022
FAMILIES <- c("msnburr", "msnburr2a")
SAMPLE_SIZES <- 10000
ALPHAS <- c(0.5, 1, 3)
TRUE_MU <- 0
TRUE_SIGMA <- 1
N_ITER <- 4000                    # with n = 10000 the posterior is narrow and
N_BURNIN <- 2000                   # the chains settle quickly; R-hat and ESS
N_CHAINS <- 4                     # are still checked for every replication
OUT_DIR <- "analysis"
CHECKPOINT_DIR <- file.path(OUT_DIR, "simulation-checkpoints-4000")
PILOT <- FALSE                    # TRUE: one replication per scenario, timed
SCENARIOS_TO_RUN <- NULL          # e.g. 1:3 to split the work across machines
STOP_AT <- NA                     # e.g. "2026-10-11 15:15": no new replication
# starts after this time

# Command-line use (one process per scenario, see run-simulation-parallel.R):
#   Rscript analysis/simulation-study.R <scenario> [stop_at]
#   Rscript analysis/simulation-study.R summary
args <- commandArgs(trailingOnly = TRUE)
SUMMARY_ONLY <- length(args) > 0 && args[1] == "summary"
if (length(args) > 0 && !SUMMARY_ONLY) SCENARIOS_TO_RUN <- as.integer(args[1])
if (length(args) > 1) STOP_AT <- args[2]
if (SUMMARY_ONLY) SCENARIOS_TO_RUN <- integer(0)
stop_time <- if (is.na(STOP_AT)) Inf else as.POSIXct(STOP_AT)

dir.create(CHECKPOINT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Models -----------------------------------------------------------------

model_code <- list(
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
  })
)

simulators <- list(msnburr = rmsnburr, msnburr2a = rmsnburr2a)

# ---- Run --------------------------------------------------------------------

results <- list()
scenario_id <- 0

for (family in FAMILIES) {
  for (n in SAMPLE_SIZES) {
    ids_here <- scenario_id + seq_along(ALPHAS)
    if (!is.null(SCENARIOS_TO_RUN) && !any(ids_here %in% SCENARIOS_TO_RUN)) {
      scenario_id <- scenario_id + length(ALPHAS)
      next
    }
    model <- nimbleModel(
      model_code[[family]],
      constants = list(N = n),
      data = list(y = rep(0, n)),
      inits = list(mu = 0, sigma = 1, alpha = 1)
    )
    conf <- configure_mcmc_neonorm(model, sampler = "slice")
    mcmc <- buildMCMC(conf)
    cmodel <- compileNimble(model)
    cmcmc <- compileNimble(mcmc, project = model)

    for (alpha in ALPHAS) {
      scenario_id <- scenario_id + 1
      if (!is.null(SCENARIOS_TO_RUN) && !scenario_id %in% SCENARIOS_TO_RUN) {
        next
      }
      message(sprintf("Scenario %d: %s, n = %d, alpha = %g",
                      scenario_id, family, n, alpha))
      truth <- c(mu = TRUE_MU, sigma = TRUE_SIGMA, alpha = alpha)

      for (rep in seq_len(if (PILOT) 1 else N_REP)) {
        ckpt <- file.path(CHECKPOINT_DIR, sprintf(
          "s%02d-%s-n%d-a%g-r%03d.rds", scenario_id, family, n, alpha, rep))
        if (!PILOT && file.exists(ckpt)) {
          results[[length(results) + 1]] <- readRDS(ckpt)
          next
        }
        if (Sys.time() > stop_time) {
          message("  STOP_AT reached; stopping before rep ", rep)
          break
        }
        t_start <- Sys.time()
        seed <- scenario_id * 100000 + rep

        set.seed(seed)
        y <- replicate(
          n, simulators[[family]](1, TRUE_MU, TRUE_SIGMA, alpha)
        )

        # New data in both the R model (used for initial values) and the
        # compiled model (used by the MCMC).
        model$setData(list(y = y))
        cmodel$y <- y

        inits <- prior_inits_neonorm(model, nchains = N_CHAINS, seed = seed)
        samples <- runMCMC(
          cmcmc,
          niter = N_ITER,
          nburnin = N_BURNIN,
          nchains = N_CHAINS,
          inits = inits,
          setSeed = seed + seq_len(N_CHAINS),
          progressBar = FALSE
        )
        tab <- summarise_neonorm(samples, warn = FALSE)

        res <- data.frame(
          family = family,
          n = n,
          alpha_true = alpha,
          rep = rep,
          parameter = names(truth),
          true = unname(truth),
          mean = tab[names(truth), "Mean"],
          lower = tab[names(truth), "2.5%"],
          upper = tab[names(truth), "97.5%"],
          rhat = tab[names(truth), "Rhat"],
          ess_bulk = tab[names(truth), "ESS_bulk"],
          ess_tail = tab[names(truth), "ESS_tail"],
          row.names = NULL
        )
        secs <- as.numeric(Sys.time() - t_start, units = "secs")
        message(sprintf("  rep %d done in %.0f s", rep, secs))
        if (PILOT) {
          message(sprintf("  estimated time for this scenario: %.1f hours",
                          secs * N_REP / 3600))
        } else {
          saveRDS(res, ckpt)
        }
        results[[length(results) + 1]] <- res
      }
    }
  }
}

if (PILOT) {
  message("Pilot finished; set PILOT <- FALSE for the full run.")
  stop("Pilot run: no performance summary written.", call. = FALSE)
}
if (length(args) > 0 && !SUMMARY_ONLY) {
  message("Scenario ", SCENARIOS_TO_RUN, " finished; run with 'summary' ",
          "after all scenarios are done.")
  quit(save = "no")
}

# Collect every finished replication (including those from earlier sessions
# or from other machines that wrote to the same checkpoint folder).
draws <- do.call(rbind, lapply(
  list.files(CHECKPOINT_DIR, pattern = "[.]rds$", full.names = TRUE),
  readRDS
))

# ---- Performance measures (Morris et al., 2019) ------------------------------

performance <- do.call(rbind, lapply(
  split(draws, draws[c("family", "n", "alpha_true", "parameter")], drop = TRUE),
  function(d) {
    err <- d$mean - d$true
    covered <- d$lower <= d$true & d$true <= d$upper
    coverage <- mean(covered)
    data.frame(
      family = d$family[1],
      n = d$n[1],
      alpha_true = d$alpha_true[1],
      parameter = d$parameter[1],
      true = d$true[1],
      bias = mean(err),
      # Monte Carlo standard error of the bias (Morris et al., 2019)
      bias_mcse = stats::sd(err) / sqrt(nrow(d)),
      mse = mean(err^2),
      coverage = coverage,
      # Monte Carlo standard error of the coverage (Morris et al., 2019)
      coverage_mcse = sqrt(coverage * (1 - coverage) / nrow(d)),
      prop_rhat_ok = mean(d$rhat < 1.01, na.rm = TRUE),
      median_ess_bulk = stats::median(d$ess_bulk, na.rm = TRUE),
      n_rep = nrow(d)
    )
  }
))
performance <- performance[order(performance$family, performance$n,
                                 performance$alpha_true,
                                 performance$parameter), ]
rownames(performance) <- NULL

utils::write.csv(draws, file.path(OUT_DIR, "simulation-draws.csv"),
                 row.names = FALSE)
utils::write.csv(performance, file.path(OUT_DIR, "simulation-performance.csv"),
                 row.names = FALSE)

print(performance, digits = 3)
cleanup_neonorm(verbose = FALSE)
utils::sessionInfo()
