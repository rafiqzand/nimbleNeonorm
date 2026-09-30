# Convenience wrappers around NIMBLE's MCMC machinery for models that use the
# neonormal distributions, organised around the Bayesian workflow of
# Gelman et al. (2020): dispersed initial values, several chains,
# convergence diagnostics and model comparison.


.supported_samplers <- c("default", "slice", "rw", "nuts", "hmc", "mixed_nuts")

# Stochastic, non-data nodes called `alpha` (scalar or indexed): the shape
# parameter targeted by the "slice", "rw" and "mixed_nuts" strategies.
.alpha_nodes <- function(model) {
  nodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  nodes[grepl("^alpha(\\[.*\\])?$", nodes)]
}

.require_nimbleHMC <- function(sampler) {
  if (!requireNamespace("nimbleHMC", quietly = TRUE)) {
    stop(
      "Package 'nimbleHMC' is required for sampler = '", sampler, "'. ",
      "Install it with install.packages(\"nimbleHMC\").",
      call. = FALSE
    )
  }
}

.check_model <- function(model) {
  if (!methods::is(model, "modelBaseClass")) {
    stop(
      "`model` must be an uncompiled NIMBLE model created by ",
      "nimble::nimbleModel().",
      call. = FALSE
    )
  }
}

.is_single_number <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x)
}

# [REVISI-04] Helper baru untuk memvalidasi argumen logis (WAIC, summary).
.check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
}

# [REVISI-08] Helper baru: menambahkan satu sampler NUTS gabungan untuk
# [REVISI-08] `target`, memakai objek sampler (bukan nama teks).
# Add one joint NUTS sampler for `target`. The sampler object is passed
# directly, so nimbleHMC does not need to be attached for NIMBLE to find it.
.add_nuts <- function(conf, target) {
  conf$addSampler(target = target, type = nimbleHMC::sampler_NUTS)
  invisible(conf)
}

# Error unless `x` is a single whole number >= `min`.
.check_count <- function(x, name, min) {
  ok <- .is_single_number(x) && x >= min && x == round(x)
  if (!ok) {
    stop("`", name, "` must be a single whole number >= ", min, ".",
         call. = FALSE)
  }
}

# [REVISI-04] Helper baru: menyalakan opsi WAIC NIMBLE sementara lalu
# [REVISI-04] mengembalikannya. Dipakai lewat nimbleOptions (bukan argumen
# [REVISI-04] configureMCMC) agar berlaku seragam di semua jalur sampler.
# Evaluate `code` with NIMBLE's MCMCenableWAIC option temporarily set.
.with_waic_option <- function(enable, code) {
  old <- nimble::nimbleOptions("MCMCenableWAIC")
  nimble::nimbleOptions(MCMCenableWAIC = enable)
  on.exit(nimble::nimbleOptions(MCMCenableWAIC = old), add = TRUE)
  code
}


# ---- Configuration ----------------------------------------------------------

# [REVISI-04] Argumen baru `enableWAIC`. [REVISI-05] Ditambah @references
# [REVISI-05] untuk tiap strategi sampler (alasan ilmiah pemilihannya).
#' Configure an MCMC for a neonormal NIMBLE model
#'
#' Builds an MCMC configuration for a NIMBLE model, optionally replacing the
#' sampler of the shape parameter `alpha`, which often mixes poorly under
#' NIMBLE's default random-walk sampler.
#'
#' @param model An uncompiled NIMBLE model created by [nimble::nimbleModel()].
#'   For `"nuts"`, `"hmc"` and `"mixed_nuts"` it must have been created with
#'   `buildDerivs = TRUE`.
#' @param params Character vector of nodes to monitor. If `NULL` (default),
#'   all stochastic non-data nodes are monitored.
#' @param sampler Sampling strategy (case-insensitive):
#'   * `"default"`: NIMBLE's default samplers (for continuous scalar nodes,
#'     adaptive random-walk Metropolis-Hastings).
#'   * `"slice"`: univariate slice sampler on `alpha`, defaults elsewhere.
#'   * `"RW"`: adaptive random-walk Metropolis on `alpha`, defaults elsewhere.
#'   * `"NUTS"`: No-U-Turn Sampler, jointly on all continuous stochastic
#'     nodes (NIMBLE's defaults for any discrete ones), using
#'     [nimbleHMC::sampler_NUTS()].
#'   * `"HMC"`: alias for `"NUTS"`.
#'   * `"mixed_nuts"`: NUTS on `alpha` only, defaults elsewhere.
#' @param thin Thinning interval for the monitored nodes.
#' @param enableWAIC Logical; prepare the MCMC to compute WAIC
#'   (Watanabe, 2010) while it runs?
#'
#' @details
#' The `alpha`-specific strategies look for stochastic nodes literally named
#' `alpha` (e.g. `alpha` or `alpha[1]`). If there are none, a warning is
#' given and NIMBLE's defaults are kept.
#'
#' Why these strategies. Random-walk Metropolis (Metropolis et al., 1953;
#' Hastings, 1970) needs its proposal scale tuned; NIMBLE adapts it during
#' burn-in, but a shape parameter whose posterior is strongly skewed can
#' still mix slowly (Roberts et al., 1997). The slice sampler (Neal, 2003)
#' adapts its step to the local shape of the posterior without a proposal
#' scale, which suits such parameters. NUTS (Hoffman & Gelman, 2014) uses
#' gradients of the log-posterior, which NIMBLE obtains by automatic
#' differentiation; this is why the density, distribution and quantile
#' functions are built with `buildDerivs = TRUE`.
#'
#' NUTS-based strategies need the \pkg{nimbleHMC} package.
#'
#' @return An `MCMCconf` object, to be passed to [nimble::buildMCMC()].
#'
#' @references
#' Hastings, W. K. (1970). Monte Carlo sampling methods using Markov chains
#' and their applications. *Biometrika*, 57(1), 97--109.
#' \doi{10.1093/biomet/57.1.97}
#'
#' Hoffman, M. D., & Gelman, A. (2014). The No-U-Turn Sampler: Adaptively
#' setting path lengths in Hamiltonian Monte Carlo. *Journal of Machine
#' Learning Research*, 15, 1593--1623.
#'
#' Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H., &
#' Teller, E. (1953). Equation of state calculations by fast computing
#' machines. *The Journal of Chemical Physics*, 21(6), 1087--1092.
#' \doi{10.1063/1.1699114}
#'
#' Neal, R. M. (2003). Slice sampling. *The Annals of Statistics*, 31(3),
#' 705--767. \doi{10.1214/aos/1056562461}
#'
#' Roberts, G. O., Gelman, A., & Gilks, W. R. (1997). Weak convergence and
#' optimal scaling of random walk Metropolis algorithms. *The Annals of
#' Applied Probability*, 7(1), 110--120. \doi{10.1214/aoap/1034625254}
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
#' widely applicable information criterion in singular learning theory.
#' *Journal of Machine Learning Research*, 11, 3571--3594.
#'
#' @seealso [runmcmc_neonorm()], [nimble::configureMCMC()]
#'
#' @export
configure_mcmc_neonorm <- function(model,
                                   params = NULL,
                                   sampler = "default",
                                   thin = 1,
                                   enableWAIC = FALSE) {
  .check_model(model)
  .check_count(thin, "thin", min = 1)
  .check_flag(enableWAIC, "enableWAIC")

  sampler <- tolower(sampler)
  if (length(sampler) != 1L || !sampler %in% .supported_samplers) {
    stop(
      "Unknown sampler '", paste(sampler, collapse = ", "), "'. ",
      "Supported samplers are: ",
      paste(.supported_samplers, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (is.null(params)) {
    params <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  }

  # [REVISI-08] NUTS tidak lagi lewat nimbleHMC::configureHMC()/addHMC(),
  # [REVISI-08] karena keduanya mencari sampler "NUTS" berdasarkan nama di
  # [REVISI-08] search path, yang gagal jika nimbleHMC hanya dimuat lewat
  # [REVISI-08] requireNamespace() (tidak di-attach). Objek sampler
  # [REVISI-08] nimbleHMC::sampler_NUTS diberikan langsung ke addSampler().
  if (sampler %in% c("nuts", "hmc")) {
    .require_nimbleHMC(sampler)
    # NUTS jointly on all continuous stochastic nodes; NIMBLE's defaults for
    # any discrete ones, which NUTS cannot sample.
    nodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
    continuous <- nodes[!model$isDiscrete(nodes)]
    discrete <- setdiff(nodes, continuous)
    conf <- .with_waic_option(
      enableWAIC,
      nimble::configureMCMC(
        model,
        nodes = if (length(discrete) > 0L) discrete else NULL,
        monitors = params,
        thin = thin
      )
    )
    if (length(continuous) > 0L) {
      .add_nuts(conf, continuous)
    }
    return(conf)
  }

  conf <- .with_waic_option(
    enableWAIC,
    nimble::configureMCMC(model, monitors = params, thin = thin)
  )
  if (sampler == "default") {
    return(conf)
  }

  if (sampler == "mixed_nuts") {
    .require_nimbleHMC(sampler)
  }

  target <- .alpha_nodes(model)
  if (length(target) == 0L) {
    warning(
      "No stochastic node named 'alpha' found; ",
      "keeping NIMBLE's default samplers.",
      call. = FALSE
    )
    return(conf)
  }

  conf$removeSamplers(target)
  if (sampler == "mixed_nuts") {
    .add_nuts(conf, target)
  } else {
    type <- if (sampler == "slice") "slice" else "RW"
    for (node in target) {
      conf$addSampler(target = node, type = type)
    }
  }

  conf
}


# ---- Initial values ---------------------------------------------------------

# [REVISI-02] FUNGSI BARU. Alasan: sebelumnya, untuk nchains > 1 tanpa
# [REVISI-02] `inits`, rantai ke-2 dst. melanjutkan dari akhir rantai
# [REVISI-02] sebelumnya, sehingga titik awal TIDAK tersebar (overdispersed).
# [REVISI-02] Padahal R-hat mensyaratkan titik awal tersebar (Gelman &
# [REVISI-02] Rubin, 1992; Vehtari dkk., 2021). Fungsi ini membangkitkan
# [REVISI-02] nilai awal tiap rantai dari prior dan memastikan log-probabilitas
# [REVISI-02] model berhingga.
#' Draw dispersed initial values from the prior
#'
#' Simulates one set of initial values per chain from the prior of a NIMBLE
#' model, keeping only draws at which the model's log-probability is finite.
#'
#' The \eqn{\hat{R}} convergence diagnostic compares chains started from
#' points that are over-dispersed relative to the posterior (Gelman &
#' Rubin, 1992; Vehtari et al., 2021). Draws from the prior provide such
#' starting points when the prior is proper and not extremely vague; with
#' very vague priors (e.g. a normal with standard deviation 1000), supply
#' `inits` yourself or use weakly informative priors, as recommended by
#' Gelman et al. (2020).
#'
#' The model is left as it was: its values and log-probabilities are
#' restored before the function returns.
#'
#' @param model An uncompiled NIMBLE model created by [nimble::nimbleModel()].
#' @param nchains Number of sets of initial values to draw.
#' @param seed Optional single number. If given, the draws are reproducible
#'   and the global random-number state is restored afterwards.
#' @param max_tries Maximum number of prior draws per chain before giving up
#'   on finding one with a finite log-probability.
#'
#' @return A list of length `nchains`; each element is a named list of
#'   initial values, suitable for the `inits` argument of
#'   [runmcmc_neonorm()] or [nimble::runMCMC()].
#'
#' @references
#' Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation
#' using multiple sequences. *Statistical Science*, 7(4), 457--472.
#' \doi{10.1214/ss/1177011136}
#'
#' Gelman, A., Vehtari, A., Simpson, D., Margossian, C. C., Carpenter, B.,
#' Yao, Y., Kennedy, L., Gabry, J., Bürkner, P.-C., & Modrák, M. (2020).
#' *Bayesian workflow*. arXiv:2011.01808.
#' \doi{10.48550/arXiv.2011.01808}
#'
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C.
#' (2021). Rank-normalization, folding, and localization: An improved
#' \eqn{\hat{R}} for assessing convergence of MCMC. *Bayesian Analysis*,
#' 16(2), 667--718. \doi{10.1214/20-BA1221}
#'
#' @seealso [runmcmc_neonorm()]
#'
#' @examples
#' \donttest{
#' register_msnburr(verbose = FALSE)
#' code <- nimble::nimbleCode({
#'   for (i in 1:N) y[i] ~ dmsnburr(mu, sigma, alpha)
#'   mu ~ dnorm(0, sd = 5)
#'   sigma ~ dexp(1)
#'   alpha ~ dlnorm(0, sdlog = 1)
#' })
#' model <- nimble::nimbleModel(
#'   code,
#'   constants = list(N = 3),
#'   data = list(y = c(-0.2, 0.4, 1.1)),
#'   inits = list(mu = 0, sigma = 1, alpha = 1)
#' )
#' prior_inits_neonorm(model, nchains = 2, seed = 1)
#' cleanup_neonorm(verbose = FALSE)
#' }
#'
#' @export
prior_inits_neonorm <- function(model,
                                nchains = 4,
                                seed = NULL,
                                max_tries = 100) {
  .check_model(model)
  .check_count(nchains, "nchains", min = 1)
  .check_count(max_tries, "max_tries", min = 1)
  if (!is.null(seed) && !.is_single_number(seed)) {
    stop("`seed` must be NULL or a single number.", call. = FALSE)
  }

  param_nodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  if (length(param_nodes) == 0L) {
    stop("The model has no stochastic non-data nodes to initialise.",
         call. = FALSE)
  }
  # All non-data nodes in topological order: simulating a deterministic node
  # recalculates it, so dependent priors see the new parent values.
  sim_nodes <- model$getNodeNames(includeData = FALSE)
  vars <- unique(model$getVarNames(nodes = param_nodes))

  # Restore the model exactly as we found it.
  all_nodes <- model$getNodeNames()
  saved <- nimble::values(model, all_nodes)
  on.exit({
    nimble::values(model, all_nodes) <- saved
    model$calculate()
  }, add = TRUE)

  draw_one <- function(chain) {
    for (i in seq_len(max_tries)) {
      model$simulate(sim_nodes)
      if (is.finite(model$calculate())) {
        return(stats::setNames(lapply(vars, function(v) model[[v]]), vars))
      }
    }
    stop(
      "Could not draw initial values with a finite log-probability for ",
      "chain ", chain, " after ", max_tries, " prior draws. ",
      "Check the priors or supply `inits` yourself.",
      call. = FALSE
    )
  }

  draw_all <- function() lapply(seq_len(nchains), draw_one)
  if (is.null(seed)) draw_all() else withr::with_seed(seed, draw_all())
}


# ---- Posterior summary and convergence diagnostics --------------------------

# Coerce a matrix, a list of matrices or a coda::mcmc.list to a list of
# draw matrices (iterations x variables), one per chain.
.as_chain_list <- function(samples) {
  if (inherits(samples, "mcmc.list")) {
    samples <- lapply(samples, as.matrix)
  } else if (is.matrix(samples)) {
    samples <- list(samples)
  }
  ok <- is.list(samples) && length(samples) > 0L &&
    all(vapply(samples, is.matrix, logical(1)))
  if (!ok) {
    stop(
      "`samples` must be a matrix of draws, a list of such matrices ",
      "(one per chain) or a coda::mcmc.list.",
      call. = FALSE
    )
  }
  vars <- colnames(samples[[1]])
  same <- all(vapply(samples, function(m) {
    identical(colnames(m), vars) && nrow(m) == nrow(samples[[1]])
  }, logical(1)))
  if (is.null(vars) || !same) {
    stop("All chains must have the same named columns and length.",
         call. = FALSE)
  }
  samples
}

# [REVISI-03] FUNGSI BARU (diekspor). Menggantikan .summarise_samples() lama
# [REVISI-03] yang hanya berisi Mean, SD, dan kuantil. Sekarang ditambah
# [REVISI-03] R-hat rank-normalized, ESS bulk, dan ESS tail dari paket
# [REVISI-03] `posterior`, sesuai Vehtari dkk. (2021), plus peringatan
# [REVISI-03] otomatis jika R-hat > 1.01 atau ESS < 100 per rantai.
#' Posterior summary with convergence diagnostics
#'
#' Summarises MCMC draws with the posterior mean, standard deviation and
#' quantiles, together with the convergence diagnostics recommended by
#' Vehtari et al. (2021): the rank-normalised split-\eqn{\hat{R}}, and the
#' bulk and tail effective sample sizes (ESS).
#'
#' Following Vehtari et al. (2021), a warning is given when any
#' \eqn{\hat{R}} exceeds `rhat_threshold` (default 1.01) or when any bulk or
#' tail ESS is below `ess_per_chain` (default 100) times the number of
#' chains. The diagnostics are computed with [posterior::rhat()],
#' [posterior::ess_bulk()] and [posterior::ess_tail()]; nodes whose draws are
#' constant get `NA`.
#'
#' @param samples Draws as a matrix (iterations x nodes, one chain), a list
#'   of such matrices (one per chain) or a `coda::mcmc.list`.
#' @param probs Quantiles to report.
#' @param warn Logical; warn when the diagnostics indicate problems?
#' @param rhat_threshold,ess_per_chain Thresholds used for the warnings.
#'
#' @return A numeric matrix with one row per node and columns `Mean`, `SD`,
#'   one column per quantile (e.g. `2.5%`), `Rhat`, `ESS_bulk` and
#'   `ESS_tail`. Mean, SD and quantiles are computed from all chains pooled.
#'
#' @references
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C.
#' (2021). Rank-normalization, folding, and localization: An improved
#' \eqn{\hat{R}} for assessing convergence of MCMC. *Bayesian Analysis*,
#' 16(2), 667--718. \doi{10.1214/20-BA1221}
#'
#' Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation
#' using multiple sequences. *Statistical Science*, 7(4), 457--472.
#' \doi{10.1214/ss/1177011136}
#'
#' @seealso [runmcmc_neonorm()], which calls it when `summary = TRUE`.
#'
#' @examples
#' set.seed(1)
#' chains <- lapply(1:4, function(i) {
#'   cbind(mu = rnorm(500), sigma = rexp(500))
#' })
#' summarise_neonorm(chains)
#'
#' @export
summarise_neonorm <- function(samples,
                              probs = c(0.025, 0.5, 0.975),
                              warn = TRUE,
                              rhat_threshold = 1.01,
                              ess_per_chain = 100) {
  chains <- .as_chain_list(samples)
  .check_flag(warn, "warn")
  if (!is.numeric(probs) || anyNA(probs) || any(probs < 0 | probs > 1)) {
    stop("`probs` must be numbers between 0 and 1.", call. = FALSE)
  }

  vars <- colnames(chains[[1]])
  nchains <- length(chains)

  rows <- lapply(vars, function(v) {
    # iterations x chains, the layout expected by the posterior package
    draws <- vapply(chains, function(ch) as.numeric(ch[, v]),
                    numeric(nrow(chains[[1]])))
    draws <- matrix(draws, ncol = nchains)
    pooled <- as.vector(draws)
    c(
      Mean = mean(pooled),
      SD = stats::sd(pooled),
      stats::setNames(
        stats::quantile(pooled, probs, names = FALSE),
        paste0(probs * 100, "%")
      ),
      Rhat = .safe_diag(posterior::rhat, draws),
      ESS_bulk = .safe_diag(posterior::ess_bulk, draws),
      ESS_tail = .safe_diag(posterior::ess_tail, draws)
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- vars

  if (warn) {
    .warn_diagnostics(out, nchains, rhat_threshold, ess_per_chain)
  }
  out
}

# Constant or too-short chains make the posterior diagnostics return NA or
# warn; report NA quietly instead.
.safe_diag <- function(f, draws) {
  value <- tryCatch(suppressWarnings(f(draws)), error = function(e) NA_real_)
  if (length(value) != 1L || !is.finite(value)) NA_real_ else value
}

.warn_diagnostics <- function(tab, nchains, rhat_threshold, ess_per_chain) {
  high_rhat <- rownames(tab)[!is.na(tab[, "Rhat"]) &
                               tab[, "Rhat"] > rhat_threshold]
  min_ess <- ess_per_chain * nchains
  low_ess <- rownames(tab)[
    (!is.na(tab[, "ESS_bulk"]) & tab[, "ESS_bulk"] < min_ess) |
      (!is.na(tab[, "ESS_tail"]) & tab[, "ESS_tail"] < min_ess)
  ]
  if (length(high_rhat) > 0L) {
    warning(
      "R-hat > ", rhat_threshold, " for: ",
      paste(high_rhat, collapse = ", "),
      ". The chains have not converged; run them longer or reparameterise.",
      call. = FALSE
    )
  }
  if (length(low_ess) > 0L) {
    warning(
      "Bulk or tail ESS < ", min_ess, " (", ess_per_chain,
      " per chain) for: ", paste(low_ess, collapse = ", "),
      ". Posterior summaries may be unreliable; run more iterations.",
      call. = FALSE
    )
  }
  invisible(NULL)
}


# ---- Running the MCMC -------------------------------------------------------

# [REVISI-02] Perubahan di runmcmc_neonorm():
# [REVISI-02]  (a) default nchains 1 -> 4 (Vehtari dkk., 2021, menyarankan
# [REVISI-02]      minimal 4 rantai agar R-hat bermakna);
# [REVISI-02]  (b) jika inits = NULL dan nchains > 1, nilai awal diambil dari
# [REVISI-02]      prior_inits_neonorm() (titik awal tersebar);
# [REVISI-02]  (c) cek log-probabilitas awal berhingga SEBELUM kompilasi,
# [REVISI-02]      dengan pesan galat yang jelas (menggantikan peran nimStop).
# [REVISI-03] Ringkasan memakai summarise_neonorm() (R-hat & ESS).
# [REVISI-04] Argumen baru WAIC; hasil disimpan di fit$WAIC.
# [REVISI-06] Objek hasil kini berkelas "neonorm_fit" dan menyimpan model
# [REVISI-06] terkompilasi, agar posterior_predict_neonorm() bisa dipakai.
#' Compile and run an MCMC for a neonormal NIMBLE model
#'
#' Configures (with [configure_mcmc_neonorm()]), builds, compiles and runs an
#' MCMC, and returns the posterior draws with a summary that includes
#' convergence diagnostics and, optionally, WAIC.
#'
#' @inheritParams configure_mcmc_neonorm
#' @param niter Total number of iterations per chain, including burn-in.
#' @param nburnin Number of initial iterations to discard. Must be smaller
#'   than `niter`. Defaults to `floor(niter / 2)`.
#' @param nchains Number of chains. At least four are recommended so that
#'   \eqn{\hat{R}} is informative (Vehtari et al., 2021).
#' @param summary Logical; also return a posterior summary with convergence
#'   diagnostics (see [summarise_neonorm()])?
#' @param setSeed Controls reproducibility. `TRUE` (default) uses seeds
#'   `123, 124, ...` for chains `1, 2, ...`; a number `s` uses seeds
#'   `s, s + 1, ...`; `FALSE` sets no seed.
#' @param inits Optional initial values, passed to [nimble::runMCMC()]: a
#'   list, a list of lists (one per chain) or a function. If `NULL` and
#'   `nchains > 1`, dispersed initial values are drawn from the prior with
#'   [prior_inits_neonorm()]; if `NULL` and `nchains = 1`, the chain starts
#'   from the values currently stored in `model`.
#' @param WAIC Logical; compute the widely applicable information criterion
#'   (Watanabe, 2010) for model comparison? Lower is better.
#'
#' @details
#' Burn-in is discarded before thinning, so each chain returns
#' `floor((niter - nburnin) / thin)` draws.
#'
#' Before compiling, the function checks that the model's log-probability
#' at the starting values is finite, and stops with an informative message
#' otherwise.
#'
#' The model is compiled on every call. For repeated runs on the same model
#' (e.g. a simulation study), build and compile once with
#' [configure_mcmc_neonorm()], [nimble::buildMCMC()] and
#' [nimble::compileNimble()], then call [nimble::runMCMC()] directly, or
#' reuse `fit$compiled$mcmc`.
#'
#' @return An object of class `neonorm_fit`: a list with elements
#'   * `samples`: a matrix of draws (one column per monitored node) for a
#'     single chain, or a `coda::mcmc.list` for several chains;
#'   * `summary`: only if `summary = TRUE`, the matrix returned by
#'     [summarise_neonorm()];
#'   * `WAIC`: only if `WAIC = TRUE`, NIMBLE's WAIC result (a list with
#'     elements `WAIC`, `lppd` and `pWAIC`);
#'   * `config`: the `MCMCconf` object used;
#'   * `inits`: the initial values used (`NULL` if taken from `model`);
#'   * `model`: the uncompiled model;
#'   * `compiled`: a list with the compiled `model` and `mcmc`;
#'   * `settings`: the run settings (`niter`, `nburnin`, `thin`, `nchains`,
#'     `sampler`).
#'
#' @references
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C.
#' (2021). Rank-normalization, folding, and localization: An improved
#' \eqn{\hat{R}} for assessing convergence of MCMC. *Bayesian Analysis*,
#' 16(2), 667--718. \doi{10.1214/20-BA1221}
#'
#' Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. *Statistics
#' and Computing*, 27(5), 1413--1432. \doi{10.1007/s11222-016-9696-4}
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
#' widely applicable information criterion in singular learning theory.
#' *Journal of Machine Learning Research*, 11, 3571--3594.
#'
#' @seealso [configure_mcmc_neonorm()], [prior_inits_neonorm()],
#'   [summarise_neonorm()], [posterior_predict_neonorm()],
#'   [nimble::runMCMC()]
#'
#' @examples
#' \donttest{
#' register_msnburr(verbose = FALSE)
#'
#' code <- nimble::nimbleCode({
#'   for (i in 1:N) {
#'     y[i] ~ dmsnburr(mu, sigma, alpha)
#'   }
#'   mu ~ dnorm(0, sd = 5)
#'   sigma ~ dexp(1)
#'   alpha ~ dlnorm(0, sdlog = 1)
#' })
#'
#' set.seed(1)
#' y <- replicate(50, rmsnburr(1, mu = 0, sigma = 1, alpha = 3))
#' model <- nimble::nimbleModel(
#'   code,
#'   constants = list(N = length(y)),
#'   data = list(y = y),
#'   inits = list(mu = 0, sigma = 1, alpha = 1)
#' )
#'
#' fit <- runmcmc_neonorm(model, niter = 2000, sampler = "slice", WAIC = TRUE)
#' fit
#'
#' cleanup_neonorm(verbose = FALSE)
#' }
#'
#' @export
runmcmc_neonorm <- function(model,
                            niter = 10000,
                            nburnin = NULL,
                            thin = 1,
                            nchains = 4,
                            params = NULL,
                            sampler = "default",
                            summary = TRUE,
                            setSeed = TRUE,
                            inits = NULL,
                            WAIC = FALSE) {
  .check_model(model)
  .check_count(niter, "niter", min = 1)
  if (is.null(nburnin)) {
    nburnin <- floor(niter / 2)
  }
  .check_count(nburnin, "nburnin", min = 0)
  if (nburnin >= niter) {
    stop("`nburnin` must be smaller than `niter`.", call. = FALSE)
  }
  .check_count(thin, "thin", min = 1)
  .check_count(nchains, "nchains", min = 1)
  .check_flag(summary, "summary")
  .check_flag(WAIC, "WAIC")

  if (isFALSE(setSeed)) {
    seeds <- FALSE
  } else if (isTRUE(setSeed) || .is_single_number(setSeed)) {
    base_seed <- if (isTRUE(setSeed)) 123 else setSeed
    # One distinct seed per chain; a single seed would make chains that
    # start from the same values identical.
    seeds <- base_seed + seq_len(nchains) - 1
  } else {
    stop("`setSeed` must be TRUE, FALSE or a single number.", call. = FALSE)
  }

  # Starting values: dispersed draws from the prior for several chains, the
  # model's current values for a single chain.
  if (is.null(inits)) {
    if (nchains > 1L) {
      inits_seed <- if (isFALSE(seeds)) NULL else seeds[1] + 10000
      inits <- prior_inits_neonorm(model, nchains = nchains,
                                   seed = inits_seed)
    } else {
      lp <- model$calculate()
      if (!is.finite(lp)) {
        stop(
          "The model's log-probability at its current values is ", lp,
          ". Check the initial values (e.g. sigma > 0 and alpha > 0) ",
          "or supply `inits`.",
          call. = FALSE
        )
      }
    }
  }

  conf <- configure_mcmc_neonorm(
    model,
    params = params,
    sampler = sampler,
    thin = thin,
    enableWAIC = WAIC
  )
  mcmc <- nimble::buildMCMC(conf)
  cmodel <- nimble::compileNimble(model)
  cmcmc <- nimble::compileNimble(mcmc, project = model)

  run_args <- list(
    cmcmc,
    niter = niter,
    nburnin = nburnin,
    nchains = nchains,
    setSeed = seeds,
    WAIC = WAIC
  )
  if (!is.null(inits)) {
    run_args$inits <- inits
  }
  # Thinning is taken from `conf`. Multiple chains come back as a list of
  # matrices; with WAIC = TRUE the draws and WAIC come back in a list.
  out <- do.call(nimble::runMCMC, run_args)
  samples <- if (WAIC) out$samples else out

  result <- list(samples = samples)

  if (summary) {
    result$summary <- summarise_neonorm(samples)
  }
  if (WAIC) {
    result$WAIC <- out$WAIC
  }

  if (nchains > 1L) {
    result$samples <- coda::mcmc.list(
      lapply(samples, coda::mcmc, thin = thin)
    )
  }

  result$config <- conf
  result$inits <- inits
  result$model <- model
  result$compiled <- list(model = cmodel, mcmc = cmcmc)
  result$settings <- list(
    niter = niter, nburnin = nburnin, thin = thin,
    nchains = nchains, sampler = tolower(sampler)
  )

  structure(result, class = "neonorm_fit")
}


# [REVISI-06] METHOD BARU: print() yang ringkas untuk objek neonorm_fit,
# [REVISI-06] agar model terkompilasi di dalamnya tidak ikut tercetak.
#' @export
print.neonorm_fit <- function(x, digits = 3, ...) {
  s <- x$settings
  cat(
    "nimbleNeonorm MCMC fit: ", s$nchains, " chain(s), ",
    s$niter, " iterations (", s$nburnin, " burn-in, thin = ", s$thin,
    "), sampler = \"", s$sampler, "\"\n\n",
    sep = ""
  )
  if (!is.null(x$summary)) {
    print(round(x$summary, digits))
  }
  if (!is.null(x$WAIC)) {
    cat("\nWAIC = ", format(x$WAIC$WAIC, digits = digits + 3),
        "  (pWAIC = ", format(x$WAIC$pWAIC, digits = digits), ")\n",
        sep = "")
  }
  invisible(x)
}
