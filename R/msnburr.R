# MSNBurr distribution: density, distribution, quantile and random
# generation, written as nimbleFunctions so they can be used both from R and
# inside NIMBLE models (after `register_msnburr()`).

#' @include utils-numeric.R
NULL

# [REVISI-01] Dokumentasi: kalimat "Invalid parameters ... raise an error"
# [REVISI-01] diganti karena sekarang parameter tidak valid menghasilkan NaN.
# [REVISI-05] Dokumentasi: ditambah paragraf sifat kemiringan (alpha < 1,
# [REVISI-05] = 1, > 1) dan blok @references (Choir 2020, Iriawan 2012,
# [REVISI-05] Burr 1942, Maechler 2012, Devroye 1986).
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
#' deviation \eqn{\sigma}. The distribution is symmetric (logistic) for
#' \eqn{\alpha = 1}, skewed to the left for \eqn{\alpha < 1} and to the
#' right for \eqn{\alpha > 1}; it inherits from the Burr II distribution a
#' stronger capacity for left skewness than for right skewness.
#'
#' All computations are carried out on the log scale with
#' [log1pexp_nimble()] and [log1mexp_nimble()], so tail probabilities and
#' log-densities stay accurate far into the tails. Random draws use
#' inverse-transform sampling, \eqn{X = Q(U)} with
#' \eqn{U \sim \mathrm{Uniform}(0, 1)}.
#'
#' Following NIMBLE conventions, each function works on scalars and the
#' logical flags `log`, `lower.tail` and `log.p` are integers (`0` = `FALSE`,
#' non-zero = `TRUE`). `rmsnburr()` returns a single draw and only accepts
#' `n = 1`.
#'
#' Invalid inputs return `NaN`, as R's and NIMBLE's built-in distributions
#' do: non-finite or `NaN` parameters, `sigma <= 0`, `alpha <= 0`, a
#' probability outside \eqn{[0, 1]} or a positive log-probability.
#' One exception: when the d, p and q functions are called from R without
#' compilation, a `NaN` argument stops with R's "missing value" error,
#' because their `NaN` check (`x != x`) only works in compiled code;
#' `is.na()` cannot be used there as it has no derivative support in
#' NIMBLE. Inside compiled models, `NaN` arguments give `NaN`. Inside
#' an MCMC, a `NaN` log-density makes the sampler reject the proposal
#' instead of stopping the run.
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
#'   (`pmsnburr()`), quantile (`qmsnburr()`) or random draw (`rmsnburr()`);
#'   `NaN` for invalid inputs.
#'
#' @references
#' Burr, I. W. (1942). Cumulative frequency functions. *The Annals of
#' Mathematical Statistics*, 13(2), 215--232. \doi{10.1214/aoms/1177731607}
#'
#' Choir, A. S. (2020). *Distribusi neo-normal baru dan karakteristiknya*
#' \[Doctoral dissertation\]. Institut Teknologi Sepuluh Nopember.
#'
#' Devroye, L. (1986). *Non-Uniform Random Variate Generation*. Springer.
#' \doi{10.1007/978-1-4613-8643-8}
#'
#' Iriawan, N. (2012). *Pemodelan dan Analisis Data-Driven*. ITS Press.
#'
#' Mächler, M. (2012). *Accurately computing* \eqn{\log(1 - \exp(-|a|))}.
#' Vignette of the \pkg{Rmpfr} package.
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
#' # Invalid parameters give NaN rather than an error
#' dmsnburr(0, mu = 0, sigma = -1, alpha = 1)
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

    # [REVISI-01] Sebelumnya: nimStop() untuk parameter tidak valid.
    # [REVISI-01] Sekarang: mengembalikan NaN seperti dnorm(0, sd = -1) di R
    # [REVISI-01] dan distribusi bawaan NIMBLE. Ditambah cek NaN (x != x)
    # [REVISI-01] yang sebelumnya belum ada di fungsi d.
    # Invalid inputs give NaN, as R's and NIMBLE's built-in densities do.
    # x != x is TRUE only for NaN in compiled code. is.na() cannot be used
    # here: NIMBLE's automatic differentiation (buildDerivs) does not
    # support it for double arguments.
    if (x != x | mu != mu | sigma != sigma | alpha != alpha) return(NaN)
    if (abs(mu) == Inf | sigma == Inf | alpha == Inf) return(NaN)
    if (sigma <= 0 | alpha <= 0) return(NaN)

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

    # [REVISI-01] nimStop() diganti return(NaN); cek NaN digabung satu baris.
    # x != x is TRUE only for NaN in compiled code. is.na() cannot be used
    # here: NIMBLE's automatic differentiation (buildDerivs) does not
    # support it for double arguments.
    if (q != q | mu != mu | sigma != sigma | alpha != alpha) return(NaN)
    if (abs(mu) == Inf | sigma == Inf | alpha == Inf) return(NaN)
    if (sigma <= 0 | alpha <= 0) return(NaN)

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

    # [REVISI-01] nimStop() diganti return(NaN), termasuk untuk p di luar
    # [REVISI-01] [0, 1] dan log(p) > 0 (sama seperti qnorm(1.1) di R).
    # x != x is TRUE only for NaN in compiled code. is.na() cannot be used
    # here: NIMBLE's automatic differentiation (buildDerivs) does not
    # support it for double arguments.
    if (p != p | mu != mu | sigma != sigma | alpha != alpha) return(NaN)
    if (abs(mu) == Inf | sigma == Inf | alpha == Inf) return(NaN)
    if (sigma <= 0 | alpha <= 0) return(NaN)

    # Convert the input, whatever its tail and scale, to log F (lower tail),
    # returning early at the ends of the support.
    log_cdf <- 0
    if (log.p != 0) {
      if (p > 0) return(NaN)
      if (lower.tail != 0) {
        if (p == -Inf) return(-Inf)
        log_cdf <- p
      } else {
        if (p == -Inf) return(Inf)
        log_cdf <- log1mexp_nimble(-p)
      }
    } else {
      if (p < 0 | p > 1) return(NaN)
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

    # NIMBLE simulation functions return one draw per call (NIMBLE User
    # Manual, "Creating user-defined distributions"). A wrong `n` is a
    # programming error, so it still stops.
    if (n != 1) nimStop("rmsnburr only supports n = 1")

    # [REVISI-01] Parameter tidak valid: nimStop() diganti return(NaN),
    # [REVISI-01] sama seperti rnorm(1, sd = -1) di R.
    if (is.na(mu) | is.na(sigma) | is.na(alpha)) return(NaN)
    if (abs(mu) == Inf | sigma == Inf | alpha == Inf) return(NaN)
    if (sigma <= 0 | alpha <= 0) return(NaN)

    omega <- exp(log_omega_msnburr(alpha))

    # Inverse-transform sampling, X = Q(U). runif() never returns 0 or 1,
    # so log(u) is finite.
    u <- runif(1, 0, 1)
    a <- -log(u) / alpha
    log_term <- a + log1mexp_nimble(a)

    return(mu - (sigma / omega) * (log(alpha) + log_term))
  }
)
