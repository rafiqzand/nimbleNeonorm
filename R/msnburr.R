# MSNBurr distribution: density, distribution, quantile and random
# generation, written as nimbleFunctions so they can be used both from R and
# inside NIMBLE models (after `register_msnburr()`).

#' @include utils-numeric.R
NULL
#' The MSNBurr distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the MSNBurr distribution with location `mu`, scale `sigma` and shape
#' `alpha`. All four are nimbleFunctions: they can be called from R and, once
#' registered with [register_msnburr()], used in NIMBLE model code as
#' `y ~ dmsnburr(mu, sigma, alpha)`.
#'
#' @details
#' With \eqn{z = (x - \mu)/\sigma} and
#' \eqn{\omega = (1 + 1/\alpha)^{\alpha + 1} / \sqrt{2\pi}}, the density is
#' \deqn{
#'   f(x) = \frac{\omega}{\sigma} e^{-\omega z}
#'   \left(1 + \frac{e^{-\omega z}}{\alpha}\right)^{-(\alpha + 1)},
#' }
#' the distribution function is
#' \deqn{
#'   F(x) = \left(1 + \frac{e^{-\omega z}}{\alpha}\right)^{-\alpha},
#' }
#' and the quantile function is
#' \deqn{
#'   Q(p) = \mu - \frac{\sigma}{\omega}
#'   \left[\log\alpha + \log\left(p^{-1/\alpha} - 1\right)\right].
#' }
#' The mode is \eqn{\mu}, and for every \eqn{\alpha} the density there is
#' \eqn{1/(\sigma\sqrt{2\pi})}, the same as a normal density with standard
#' deviation \eqn{\sigma}.
#'
#' All computations are carried out on the log scale with
#' [log1pexp_nimble()] and [log1mexp_nimble()], so tail probabilities and
#' log-densities stay accurate far into the tails.
#'
#' Following NIMBLE conventions, each function works on scalars and the
#' logical flags `log`, `lower.tail` and `log.p` are integers (`0` = `FALSE`,
#' non-zero = `TRUE`). `rmsnburr()` returns a single draw and only accepts
#' `n = 1`. Invalid parameters (non-finite values, `sigma <= 0` or
#' `alpha <= 0`) raise an error.
#'
#' @param x,q Numeric scalar, the point at which to evaluate the density or
#'   distribution function.
#' @param p Numeric scalar, a probability (or log-probability if `log.p` is
#'   non-zero).
#' @param n Number of draws. Must be `1`.
#' @param mu Finite location parameter.
#' @param sigma Finite, positive scale parameter.
#' @param alpha Finite, positive shape parameter.
#' @param log,log.p Integer flag; if non-zero, densities/probabilities are
#'   given (or taken) on the log scale.
#' @param lower.tail Integer flag; if non-zero (default) probabilities are
#'   \eqn{P(X \le x)}, otherwise \eqn{P(X > x)}.
#'
#' @return A numeric scalar: the density (`dmsnburr()`), probability
#'   (`pmsnburr()`), quantile (`qmsnburr()`) or random draw (`rmsnburr()`).
#'
#' @seealso [MSNBurr2a] for the mirror-image distribution,
#'   [register_msnburr()] to use it in NIMBLE models.
#'
#' @examples
#' dmsnburr(0, mu = 0, sigma = 1, alpha = 1)
#' pmsnburr(1, mu = 0, sigma = 1, alpha = 2)
#' qmsnburr(0.5, mu = 0, sigma = 1, alpha = 2)
#' set.seed(1)
#' rmsnburr(1, mu = 0, sigma = 1, alpha = 2)
#'
#' @name MSNBurr
NULL


# ---- Density ----------------------------------------------------------------

#' @rdname MSNBurr
#' @export
dmsnburr <- nimble::nimbleFunction(
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

    # log f = log(omega) - log(sigma) - omega * z
    #         - (alpha + 1) * log(1 + exp(-omega * z) / alpha)
    eta <- -omega * z - log(alpha)
    log_density <- log_omega - log(sigma) - omega * z -
      (alpha + 1) * log1pexp_nimble(eta)

    if (log != 0) return(log_density)
    return(exp(log_density))
  },
  buildDerivs = TRUE
)


# ---- Distribution function --------------------------------------------------

#' @rdname MSNBurr
#' @export
pmsnburr <- nimble::nimbleFunction(
  run = function(q = double(0),
                 mu = double(0),
                 sigma = double(0),
                 alpha = double(0),
                 lower.tail = integer(0, default = 1),
                 log.p = integer(0, default = 0)) {
    returnType(double(0))

    if (q != q) return(NaN)

    if (mu != mu) nimStop("mu must not be NaN")
    if (sigma != sigma) nimStop("sigma must not be NaN")
    if (alpha != alpha) nimStop("alpha must not be NaN")
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

    # log F = -alpha * log(1 + exp(-omega * z) / alpha)
    eta <- -omega * ((q - mu) / sigma) - log(alpha)
    log_cdf <- -alpha * log1pexp_nimble(eta)

    if (lower.tail != 0) {
      if (log.p != 0) return(log_cdf)
      return(exp(log_cdf))
    }

    # log(1 - F) from log F without forming 1 - F, which would cancel when
    # F is close to one.
    log_sf <- log1mexp_nimble(-log_cdf)
    if (log.p != 0) return(log_sf)
    return(exp(log_sf))
  },
  buildDerivs = TRUE
)


# ---- Quantile function ------------------------------------------------------

#' @rdname MSNBurr
#' @export
qmsnburr <- nimble::nimbleFunction(
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

    # Convert the input, whatever its tail and scale, to log F (lower tail),
    # returning early at the ends of the support.
    log_cdf <- 0
    if (log.p != 0) {
      if (p > 0) nimStop("log(p) must be <= 0")
      if (lower.tail != 0) {
        if (p == -Inf) return(-Inf)
        log_cdf <- p
      } else {
        if (p == -Inf) return(Inf)
        log_cdf <- log1mexp_nimble(-p)
      }
    } else {
      if (p < 0 | p > 1) nimStop("p must be between 0 and 1")
      if (lower.tail != 0) {
        if (p == 0) return(-Inf)
        if (p == 1) return(Inf)
        log_cdf <- log(p)
      } else {
        if (p == 0) return(Inf)
        if (p == 1) return(-Inf)
        log_cdf <- log1mexp_nimble(-log(p)) # log(1 - p)
      }
    }

    omega <- exp(log_omega_msnburr(alpha))

    # With a = -log(F) / alpha > 0, log(F^(-1/alpha) - 1) = log(exp(a) - 1),
    # evaluated stably as a + log(1 - exp(-a)).
    a <- -log_cdf / alpha
    log_term <- a + log1mexp_nimble(a)

    return(mu - (sigma / omega) * (log(alpha) + log_term))
  },
  buildDerivs = TRUE
)


# ---- Random generation ------------------------------------------------------

#' @rdname MSNBurr
#' @export
rmsnburr <- nimble::nimbleFunction(
  run = function(n = integer(0),
                 mu = double(0),
                 sigma = double(0),
                 alpha = double(0)) {
    returnType(double(0))

    # NIMBLE simulation functions return one draw per call.
    if (n != 1) nimStop("rmsnburr only supports n = 1")

    if (mu != mu) nimStop("mu must not be NaN")
    if (sigma != sigma) nimStop("sigma must not be NaN")
    if (alpha != alpha) nimStop("alpha must not be NaN")
    if (mu == Inf | mu == -Inf) nimStop("mu must be finite")
    if (sigma == Inf | sigma == -Inf) nimStop("sigma must be finite")
    if (alpha == Inf | alpha == -Inf) nimStop("alpha must be finite")
    if (sigma <= 0) nimStop("sigma must be > 0")
    if (alpha <= 0) nimStop("alpha must be > 0")

    omega <- exp(log_omega_msnburr(alpha))

    # Inverse-transform sampling, X = Q(U). runif() never returns 0 or 1,
    # so log(u) is finite.
    u <- runif(1, 0, 1)
    a <- -log(u) / alpha
    log_term <- a + log1mexp_nimble(a)

    return(mu - (sigma / omega) * (log(alpha) + log_term))
  }
)
