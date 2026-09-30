# [REVISI-06] FILE BARU. Prior dan posterior predictive simulation, dua
# [REVISI-06] langkah Bayesian workflow (Gelman dkk., 2020; Gabry dkk., 2019)
# [REVISI-06] yang sebelumnya belum difasilitasi paket ini.

# Prior and posterior predictive simulation for NIMBLE models, two checks of
# the Bayesian workflow (Gelman et al., 2020; Gabry et al., 2019).


# Pooled posterior draws (iterations x nodes) from a matrix or mcmc.list.
.pool_draws <- function(samples) {
  do.call(rbind, .as_chain_list(samples))
}

.check_seed <- function(seed) {
  if (!is.null(seed) && !.is_single_number(seed)) {
    stop("`seed` must be NULL or a single number.", call. = FALSE)
  }
}


#' Prior predictive simulation
#'
#' Simulates parameters from the prior and then data from the likelihood,
#' `nsim` times, to check whether the prior implies plausible data before any
#' data are used (Gelman et al., 2020, Section 2.4; Gabry et al., 2019).
#'
#' The simulation runs on the uncompiled model, so it needs no C++
#' compilation; for large models or many simulations it can be slow. The
#' model's values, including its data, are restored before the function
#' returns.
#'
#' @param model An uncompiled NIMBLE model created by [nimble::nimbleModel()]
#'   with its data set.
#' @param nsim Number of prior predictive draws.
#' @param seed Optional single number for reproducible draws; the global
#'   random-number state is restored afterwards.
#'
#' @return A list with two matrices, each with `nsim` rows:
#'   * `parameters`: the simulated values of the stochastic non-data nodes;
#'   * `data`: the simulated data, one column per data node, in the same
#'     order as `model$getNodeNames(dataOnly = TRUE)`.
#'
#' @references
#' Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., & Gelman, A. (2019).
#' Visualization in Bayesian workflow. *Journal of the Royal Statistical
#' Society: Series A*, 182(2), 389--402. \doi{10.1111/rssa.12378}
#'
#' Gelman, A., Vehtari, A., Simpson, D., Margossian, C. C., Carpenter, B.,
#' Yao, Y., Kennedy, L., Gabry, J., Bürkner, P.-C., & Modrák, M. (2020).
#' *Bayesian workflow*. arXiv:2011.01808.
#' \doi{10.48550/arXiv.2011.01808}
#'
#' @seealso [posterior_predict_neonorm()]
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
#'   constants = list(N = 20),
#'   data = list(y = rep(0, 20)),
#'   inits = list(mu = 0, sigma = 1, alpha = 1)
#' )
#' prior <- prior_predict_neonorm(model, nsim = 50, seed = 1)
#' summary(as.vector(prior$data))
#' cleanup_neonorm(verbose = FALSE)
#' }
#'
#' @export
prior_predict_neonorm <- function(model, nsim = 500, seed = NULL) {
  .check_model(model)
  .check_count(nsim, "nsim", min = 1)
  .check_seed(seed)

  data_nodes <- model$getNodeNames(dataOnly = TRUE)
  if (length(data_nodes) == 0L) {
    stop("The model has no data nodes to simulate.", call. = FALSE)
  }
  param_nodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  all_nodes <- model$getNodeNames()

  saved <- nimble::values(model, all_nodes)
  on.exit({
    nimble::values(model, all_nodes) <- saved
    model$calculate()
  }, add = TRUE)

  simulate_all <- function() {
    params <- vector("list", nsim)
    data <- vector("list", nsim)
    for (i in seq_len(nsim)) {
      # All nodes in topological order, data included.
      model$simulate(all_nodes, includeData = TRUE)
      params[[i]] <- nimble::values(model, param_nodes)
      data[[i]] <- nimble::values(model, data_nodes)
    }
    list(parameters = do.call(rbind, params), data = do.call(rbind, data))
  }

  out <- if (is.null(seed)) {
    simulate_all()
  } else {
    withr::with_seed(seed, simulate_all())
  }
  colnames(out$parameters) <- model$expandNodeNames(
    param_nodes, returnScalarComponents = TRUE
  )
  colnames(out$data) <- model$expandNodeNames(
    data_nodes, returnScalarComponents = TRUE
  )
  out
}


#' Posterior predictive simulation
#'
#' For `ndraws` posterior draws from a [runmcmc_neonorm()] fit, sets the
#' parameters to the drawn values and simulates replicated data from the
#' likelihood. Comparing replicated with observed data is the posterior
#' predictive check of the Bayesian workflow (Gelman et al., 2013,
#' Chapter 6; Gabry et al., 2019).
#'
#' The simulation uses the compiled model stored in the fit, so it is fast.
#' All stochastic non-data nodes must have been monitored (the default of
#' [runmcmc_neonorm()]), otherwise the replicated data would not condition on
#' them. The compiled model's values, including its data, are restored
#' before the function returns.
#'
#' @param fit An object returned by [runmcmc_neonorm()].
#' @param ndraws Number of posterior draws to use, taken evenly spaced from
#'   the pooled draws of all chains. At most the number of available draws.
#' @param seed Optional single number for reproducible draws; the global
#'   random-number state is restored afterwards.
#'
#' @return A matrix with `ndraws` rows and one column per data node: each
#'   row is one replicated data set \eqn{y^{rep}}. The columns are in the
#'   order of `fit$model$getNodeNames(dataOnly = TRUE)`, expanded to scalar
#'   elements. The result can be passed to, e.g., `bayesplot::ppc_dens_overlay()`.
#'
#' @references
#' Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., & Gelman, A. (2019).
#' Visualization in Bayesian workflow. *Journal of the Royal Statistical
#' Society: Series A*, 182(2), 389--402. \doi{10.1111/rssa.12378}
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., &
#' Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.). CRC Press.
#'
#' @seealso [prior_predict_neonorm()], [runmcmc_neonorm()]
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
#' set.seed(1)
#' y <- replicate(50, rmsnburr(1, mu = 0, sigma = 1, alpha = 3))
#' model <- nimble::nimbleModel(
#'   code,
#'   constants = list(N = length(y)),
#'   data = list(y = y),
#'   inits = list(mu = 0, sigma = 1, alpha = 1)
#' )
#' fit <- runmcmc_neonorm(model, niter = 2000, sampler = "slice")
#' yrep <- posterior_predict_neonorm(fit, ndraws = 100, seed = 1)
#' dim(yrep)
#' cleanup_neonorm(verbose = FALSE)
#' }
#'
#' @export
posterior_predict_neonorm <- function(fit, ndraws = 500, seed = NULL) {
  if (!inherits(fit, "neonorm_fit")) {
    stop("`fit` must be an object returned by runmcmc_neonorm().",
         call. = FALSE)
  }
  .check_count(ndraws, "ndraws", min = 1)
  .check_seed(seed)

  model <- fit$model
  cmodel <- fit$compiled$model
  draws <- .pool_draws(fit$samples)

  needed <- model$expandNodeNames(
    model$getNodeNames(stochOnly = TRUE, includeData = FALSE),
    returnScalarComponents = TRUE
  )
  missing_nodes <- setdiff(needed, colnames(draws))
  if (length(missing_nodes) > 0L) {
    stop(
      "Posterior predictive simulation needs draws of every stochastic ",
      "non-data node, but these were not monitored: ",
      paste(missing_nodes, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (ndraws > nrow(draws)) {
    stop("`ndraws` (", ndraws, ") is larger than the number of posterior ",
         "draws (", nrow(draws), ").", call. = FALSE)
  }

  # Only data and deterministic nodes are simulated; the parameters come
  # from the posterior draws.
  data_nodes <- model$getNodeNames(dataOnly = TRUE)
  deps <- model$getDependencies(needed, self = FALSE, downstream = TRUE)
  sim_nodes <- deps[model$isData(deps) | model$isDeterm(deps)]

  all_nodes <- model$getNodeNames()
  saved <- nimble::values(cmodel, all_nodes)
  on.exit({
    nimble::values(cmodel, all_nodes) <- saved
    cmodel$calculate()
  }, add = TRUE)

  rows <- unique(round(seq(1, nrow(draws), length.out = ndraws)))
  simulate_all <- function() {
    reps <- lapply(rows, function(i) {
      nimble::values(cmodel, needed) <- draws[i, needed]
      cmodel$simulate(sim_nodes, includeData = TRUE)
      nimble::values(cmodel, data_nodes)
    })
    do.call(rbind, reps)
  }

  yrep <- if (is.null(seed)) {
    simulate_all()
  } else {
    withr::with_seed(seed, simulate_all())
  }
  colnames(yrep) <- model$expandNodeNames(
    data_nodes, returnScalarComponents = TRUE
  )
  yrep
}
