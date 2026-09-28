# Convenience wrappers around NIMBLE's MCMC machinery for models that use the
# neo-normal distributions.


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

# Error unless `x` is a single whole number >= `min`.
.check_count <- function(x, name, min) {
  ok <- .is_single_number(x) && x >= min && x == round(x)
  if (!ok) {
    stop("`", name, "` must be a single whole number >= ", min, ".",
         call. = FALSE)
  }
}

# One row per monitored node: mean, SD and the 2.5/50/97.5% quantiles.
.summarise_samples <- function(samples) {
  t(apply(samples, 2, function(x) {
    q <- stats::quantile(x, c(0.025, 0.5, 0.975), names = FALSE)
    c(Mean = mean(x), SD = stats::sd(x),
      `2.5%` = q[1], `50%` = q[2], `97.5%` = q[3])
  }))
}


#' Configure an MCMC for a neo-normal NIMBLE model
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
#'   * `"default"`: NIMBLE's default samplers.
#'   * `"slice"`: slice sampler on `alpha`, defaults elsewhere.
#'   * `"RW"`: adaptive random-walk Metropolis on `alpha`, defaults elsewhere.
#'   * `"NUTS"`: NUTS on all continuous nodes, via [nimbleHMC::configureHMC()].
#'   * `"HMC"`: alias for `"NUTS"`.
#'   * `"mixed_nuts"`: NUTS on `alpha` only (via [nimbleHMC::addHMC()]),
#'     defaults elsewhere.
#' @param thin Thinning interval for the monitored nodes.
#'
#' @details
#' The `alpha`-specific strategies look for stochastic nodes literally named
#' `alpha` (e.g. `alpha` or `alpha[1]`). If there are none, a warning is
#' given and NIMBLE's defaults are kept.
#'
#' NUTS-based strategies need the \pkg{nimbleHMC} package.
#'
#' @return An `MCMCconf` object, to be passed to [nimble::buildMCMC()].
#'
#' @seealso [runmcmc_neonorm()], [nimble::configureMCMC()]
#'
#' @export
configure_mcmc_neonorm <- function(model,
                                   params = NULL,
                                   sampler = "default",
                                   thin = 1) {
  .check_model(model)
  .check_count(thin, "thin", min = 1)

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

  if (sampler %in% c("nuts", "hmc")) {
    .require_nimbleHMC(sampler)
    return(nimbleHMC::configureHMC(model, monitors = params, thin = thin))
  }

  conf <- nimble::configureMCMC(model, monitors = params, thin = thin)
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
    nimbleHMC::addHMC(conf, target = target, type = "NUTS")
  } else {
    type <- if (sampler == "slice") "slice" else "RW"
    for (node in target) {
      conf$addSampler(target = node, type = type)
    }
  }

  conf
}


#' Compile and run an MCMC for a neo-normal NIMBLE model
#'
#' Configures (with [configure_mcmc_neonorm()]), builds, compiles and runs an
#' MCMC, and returns the posterior draws with a short summary.
#'
#' @inheritParams configure_mcmc_neonorm
#' @param niter Total number of iterations per chain, including burn-in.
#' @param nburnin Number of initial iterations to discard. Must be smaller
#'   than `niter`. Defaults to `floor(niter / 2)`.
#' @param nchains Number of chains.
#' @param summary Logical; also return a posterior summary?
#' @param setSeed Controls reproducibility. `TRUE` (default) uses seeds
#'   `123, 124, ...` for chains `1, 2, ...`; a number `s` uses seeds
#'   `s, s + 1, ...`; `FALSE` sets no seed.
#' @param inits Optional initial values, passed to [nimble::runMCMC()]: a
#'   list, a list of lists (one per chain) or a function. If `NULL`, the
#'   first chain starts from the values currently stored in `model` and each
#'   later chain continues from where the previous one ended; supply `inits`
#'   for chains with dispersed starting values.
#'
#' @details
#' Burn-in is discarded before thinning, so each chain returns
#' `floor((niter - nburnin) / thin)` draws.
#'
#' The model is compiled on every call. For repeated runs on the same model,
#' build and compile once with [configure_mcmc_neonorm()],
#' [nimble::buildMCMC()] and [nimble::compileNimble()], then call
#' [nimble::runMCMC()] directly.
#'
#' @return A list with elements
#'   * `samples`: a matrix of draws (one column per monitored node) for a
#'     single chain, or a `coda::mcmc.list` for several chains;
#'   * `config`: the `MCMCconf` object used;
#'   * `summary`: only if `summary = TRUE`, a matrix with columns `Mean`,
#'     `SD`, `2.5%`, `50%` and `97.5%`, computed from all chains pooled.
#'
#' @seealso [configure_mcmc_neonorm()], [nimble::runMCMC()]
#'
#' @examples
#' \donttest{
#' register_msnburr(verbose = FALSE)
#'
#' code <- nimble::nimbleCode({
#'   for (i in 1:N) {
#'     y[i] ~ dmsnburr(mu, sigma, alpha)
#'   }
#'   mu ~ dnorm(0, sd = 10)
#'   sigma ~ dexp(1)
#'   alpha ~ dgamma(2, 1)
#' })
#'
#' y <- c(-0.4, 0.2, 1.3, 0.8, -1.1)
#' model <- nimble::nimbleModel(
#'   code,
#'   constants = list(N = length(y)),
#'   data = list(y = y),
#'   inits = list(mu = 0, sigma = 1, alpha = 1)
#' )
#'
#' fit <- runmcmc_neonorm(model, niter = 2000, sampler = "slice")
#' fit$summary
#'
#' cleanup_neonorm(verbose = FALSE)
#' }
#'
#' @export
runmcmc_neonorm <- function(model,
                            niter = 10000,
                            nburnin = NULL,
                            thin = 1,
                            nchains = 1,
                            params = NULL,
                            sampler = "default",
                            summary = TRUE,
                            setSeed = TRUE,
                            inits = NULL) {
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

  conf <- configure_mcmc_neonorm(
    model,
    params = params,
    sampler = sampler,
    thin = thin
  )
  mcmc <- nimble::buildMCMC(conf)
  nimble::compileNimble(model)
  cmcmc <- nimble::compileNimble(mcmc, project = model)

  run_args <- list(
    cmcmc,
    niter = niter,
    nburnin = nburnin,
    nchains = nchains,
    setSeed = seeds
  )
  if (!is.null(inits)) {
    run_args$inits <- inits
  }
  # Thinning is taken from `conf`. Multiple chains come back as a list of
  # matrices.
  samples <- do.call(nimble::runMCMC, run_args)

  result <- list(samples = samples, config = conf)

  if (summary) {
    pooled <- if (nchains == 1L) samples else do.call(rbind, samples)
    result$summary <- .summarise_samples(pooled)
  }

  if (nchains > 1L) {
    result$samples <- coda::mcmc.list(
      lapply(samples, coda::mcmc, thin = thin)
    )
  }

  result
}
