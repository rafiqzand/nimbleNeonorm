# MSNBurr-IIa distribution: density, distribution, quantile and random
# generation, written as nimbleFunctions so they can be used both from R and
# inside NIMBLE models (after `register_msnburr2a()`).

#' @include utils-numeric.R
NULL
#' The MSNBurr-IIa distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the MSNBurr-IIa distribution with location `mu`, scale `sigma` and
#' shape `alpha`. All four are nimbleFunctions: they can be called from R
#' and, once registered with [register_msnburr2a()], used in NIMBLE model
#' code as `y ~ dmsnburr2a(mu, sigma, alpha)`.
#'
#' @details
#' MSNBurr-IIa is the mirror image of [MSNBurr] about \eqn{\mu}: if
#' \eqn{f_{I}} is the MSNBurr density, then
#' \eqn{f_{IIa}(x) = f_{I}(2\mu - x)}.
#'
#' With \eqn{z = (x - \mu)/\sigma} and
#' \eqn{\omega = (1 + 1/\alpha)^{\alpha + 1} / \sqrt{2\pi}}, the density is
#' \deqn{
#'   f(x) = \frac{\omega}{\sigma} e^{\omega z}
#'   \left(1 + \frac{e^{\omega z}}{\alpha}\right)^{-(\alpha + 1)},
#' }
#' the survival function is
#' \deqn{
#'   S(x) = 1 - F(x) = \left(1 + \frac{e^{\omega z}}{\alpha}\right)^{-\alpha},
#' }
#' and the quantile function is
#' \deqn{
#'   Q(p) = \mu + \frac{\sigma}{\omega}
#'   \left[\log\alpha + \log\left((1 - p)^{-1/\alpha} - 1\right)\right].
#' }
#'
#' All computations are carried out on the log scale with
#' [log1pexp_nimble()] and [log1mexp_nimble()], so tail probabilities and
#' log-densities stay accurate far into the tails.
#'
#' Following NIMBLE conventions, each function works on scalars and the
#' logical flags `log`, `lower.tail` and `log.p` are integers (`0` = `FALSE`,
#' non-zero = `TRUE`). `rmsnburr2a()` returns a single draw and only accepts
#' `n = 1`. Invalid parameters (non-finite values, `sigma <= 0` or
#' `alpha <= 0`) raise an error.
#'
#' @inheritParams MSNBurr
#'
#' @return A numeric scalar: the density (`dmsnburr2a()`), probability
#'   (`pmsnburr2a()`), quantile (`qmsnburr2a()`) or random draw
#'   (`rmsnburr2a()`).
#'
#' @seealso [MSNBurr] for the mirror-image distribution,
#'   [register_msnburr2a()] to use it in NIMBLE models.
#'
#' @examples
#' dmsnburr2a(0, mu = 0, sigma = 1, alpha = 1)
#' pmsnburr2a(1, mu = 0, sigma = 1, alpha = 2)
#' qmsnburr2a(0.5, mu = 0, sigma = 1, alpha = 2)
#' set.seed(1)
#' rmsnburr2a(1, mu = 0, sigma = 1, alpha = 2)
#'
#' @name MSNBurr2a
NULL


# ---- Density ----------------------------------------------------------------

#' @rdname MSNBurr2a
#' @export
dmsnburr2a <- nimble::nimbleFunction(
  run = function(x = double(0),
                 mu = double(0),
                 sigma = double(0),
                 alpha = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))

    if (mu == Inf | mu == -Inf) nimStop("mu must be finite")
    if (sigma == Inf | sigma == -Inf) nimStop("sigma must be finite")
    if (alpha == Inf | alpha == -Inf) nimStop("alpha must be finite")
    if (sigma <= 0) nimStop("sigma must be > 0")
    if (alpha <= 0) nimStop("alpha must be > 0")

    # Handle infinite x before standardising, where Inf - Inf would give NaN.
    if (x == Inf | x == -Inf) {
      if (log != 0) return(-Inf)
      return(0)
    }

    log_omega <- log_omega_msnburr(alpha)
    omega <- exp(log_omega)
    z <- (x - mu) / sigma

    # log f = log(omega) - log(sigma) + omega * z
    #         - (alpha + 1) * log(1 + exp(omega * z) / alpha)
    eta <- omega * z - log(alpha)
    log_density <- log_omega - log(sigma) + omega * z -
      (alpha + 1) * log1pexp_nimble(eta)

    if (log != 0) return(log_density)
    return(exp(log_density))
  },
  buildDerivs = TRUE
)


# ---- Distribution function --------------------------------------------------

#' @rdname MSNBurr2a
#' @export
pmsnburr2a <- nimble::nimbleFunction(
  run = function(q = double(0),
                 mu = double(0),
                 sigma = double(0),
                 alpha = double(0),
                 lower.tail = integer(0, default = 1),
                 log.p = integer(0, default = 0)) {
    returnType(double(0))

    if (mu == Inf | mu == -Inf) nimStop("mu must be finite")
    if (sigma == Inf | sigma == -Inf) nimStop("sigma must be finite")
    if (alpha == Inf | alpha == -Inf) nimStop("alpha must be finite")
    if (sigma <= 0) nimStop("sigma must be > 0")
    if (alpha <= 0) nimStop("alpha must be > 0")

    # Limits at the ends of the support.
    if (q == -Inf) {
      if (lower.tail != 0) {
        if (log.p != 0) return(-Inf)
        return(0)
      }
      if (log.p != 0) return(0)
      return(1)
    }
    if (q == Inf) {
      if (lower.tail != 0) {
        if (log.p != 0) return(0)
        return(1)
      }
      if (log.p != 0) return(-Inf)
      return(0)
    }

    omega <- exp(log_omega_msnburr(alpha))

    # log S = -alpha * log(1 + exp(omega * z) / alpha)
    eta <- omega * ((q - mu) / sigma) - log(alpha)
    log_sf <- -alpha * log1pexp_nimble(eta)

    if (lower.tail == 0) {
      if (log.p != 0) return(log_sf)
      return(exp(log_sf))
    }

    # log F = log(1 - S) from log S without forming 1 - S, which would
    # cancel when S is close to one.
    log_cdf <- log1mexp_nimble(-log_sf)
    if (log.p != 0) return(log_cdf)
    return(exp(log_cdf))
  },
  buildDerivs = TRUE
)


# ---- Quantile function ------------------------------------------------------

#' @rdname MSNBurr2a
#' @export
qmsnburr2a <- nimble::nimbleFunction(
  run = function(p = double(0),
                 mu = double(0),
                 sigma = double(0),
                 alpha = double(0),
                 lower.tail = integer(0, default = 1),
                 log.p = integer(0, default = 0)) {
    returnType(double(0))

    if (p != p) return(NaN)

    if (mu != mu) nimStop("mu must not be NaN")
    if (sigma != sigma) nimStop("sigma must not be NaN")
    if (alpha != alpha) nimStop("alpha must not be NaN")
    if (mu == Inf | mu == -Inf) nimStop("mu must be finite")
    if (sigma == Inf | sigma == -Inf) nimStop("sigma must be finite")
    if (alpha == Inf | alpha == -Inf) nimStop("alpha must be finite")
    if (sigma <= 0) nimStop("sigma must be > 0")
    if (alpha <= 0) nimStop("alpha must be > 0")

    # Convert the input, whatever its tail and scale, to log S (upper tail),
    # returning early at the ends of the support.
    log_sf <- 0
    if (log.p != 0) {
      if (p > 0) nimStop("log(p) must be <= 0")
      if (lower.tail != 0) {
        if (p == -Inf) return(-Inf)
        log_sf <- log1mexp_nimble(-p)
      } else {
        if (p == -Inf) return(Inf)
        log_sf <- p
      }
    } else {
      if (p < 0 | p > 1) nimStop("p must be between 0 and 1")
      if (lower.tail != 0) {
        if (p == 0) return(-Inf)
        if (p == 1) return(Inf)
        log_sf <- log1mexp_nimble(-log(p)) # log(1 - p)
      } else {
        if (p == 0) return(Inf)
        if (p == 1) return(-Inf)
        log_sf <- log(p)
      }
    }

    omega <- exp(log_omega_msnburr(alpha))

    # With a = -log(S) / alpha > 0, log(S^(-1/alpha) - 1) = log(exp(a) - 1),
    # evaluated stably as a + log(1 - exp(-a)).
    a <- -log_sf / alpha
    log_term <- a + log1mexp_nimble(a)

    return(mu + (sigma / omega) * (log(alpha) + log_term))
  },
  buildDerivs = TRUE
)


# ---- Random generation ------------------------------------------------------

#' @rdname MSNBurr2a
#' @export
rmsnburr2a <- nimble::nimbleFunction(
  run = function(n = integer(0),
                 mu = double(0),
                 sigma = double(0),
                 alpha = double(0)) {
    returnType(double(0))

    # NIMBLE simulation functions return one draw per call.
    if (n != 1) nimStop("rmsnburr2a only supports n = 1")

    if (mu != mu) nimStop("mu must not be NaN")
    if (sigma != sigma) nimStop("sigma must not be NaN")
    if (alpha != alpha) nimStop("alpha must not be NaN")
    if (mu == Inf | mu == -Inf) nimStop("mu must be finite")
    if (sigma == Inf | sigma == -Inf) nimStop("sigma must be finite")
    if (alpha == Inf | alpha == -Inf) nimStop("alpha must be finite")
    if (sigma <= 0) nimStop("sigma must be > 0")
    if (alpha <= 0) nimStop("alpha must be > 0")

    omega <- exp(log_omega_msnburr(alpha))

    # Inverse-transform sampling with S = 1 - U. log1p(-u) keeps precision
    # for small u; runif() never returns 0 or 1, so the result is finite.
    u <- runif(1, 0, 1)
    a <- -log1p(-u) / alpha
    log_term <- a + log1mexp_nimble(a)

    return(mu + (sigma / omega) * (log(alpha) + log_term))
  }
)
