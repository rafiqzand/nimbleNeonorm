# MSNBurr-IIa distribution: density, distribution, quantile and random
# generation, written as nimbleFunctions so they can be used both from R and
# inside NIMBLE models (after `register_msnburr2a()`).

#' @include utils-numeric.R
NULL

# [REVISI-01] Dokumentasi: kalimat "Invalid parameters ... raise an error"
# [REVISI-01] diganti karena sekarang parameter tidak valid menghasilkan NaN.
# [REVISI-05] Dokumentasi: ditambah paragraf sifat kemiringan dan blok
# [REVISI-05] @references (Ramadani dkk. 2025, Choir 2020, Maechler 2012).
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
#' The distribution is symmetric for \eqn{\alpha = 1}, skewed to the right
#' for \eqn{\alpha < 1} and to the left for \eqn{\alpha > 1}; it inherits
#' from the Burr IIa distribution a stronger capacity for right skewness than
#' for left skewness (Ramadani et al., 2025).
#'
#' All computations are carried out on the log scale with
#' [log1pexp_nimble()] and [log1mexp_nimble()], so tail probabilities and
#' log-densities stay accurate far into the tails. Random draws use
#' inverse-transform sampling with \eqn{S = 1 - U}.
#'
#' Following NIMBLE conventions, each function works on scalars and the
#' logical flags `log`, `lower.tail` and `log.p` are integers (`0` = `FALSE`,
#' non-zero = `TRUE`). `rmsnburr2a()` returns a single draw and only accepts
#' `n = 1`.
#'
#' Invalid inputs return `NaN`, as R's and NIMBLE's built-in distributions
#' do: non-finite or `NaN` parameters, `sigma <= 0`, `alpha <= 0`, a
#' probability outside \eqn{[0, 1]} or a positive log-probability.
#' One exception: when the d, p and q functions are called from R without
#' compilation, a `NaN` argument stops with R's "missing value" error,
#' because their `NaN` check (`x != x`) only works in compiled code;
#' `is.na()` cannot be used there as it has no derivative support in
#' NIMBLE. Inside compiled models, `NaN` arguments give `NaN`.
#'
#' @inheritParams MSNBurr
#'
#' @return A numeric scalar: the density (`dmsnburr2a()`), probability
#'   (`pmsnburr2a()`), quantile (`qmsnburr2a()`) or random draw
#'   (`rmsnburr2a()`); `NaN` for invalid inputs.
#'
#' @references
#' Ramadani, E. P., Choir, A. S., Pravitasari, A. A., & Paraguison, J.
#' (2025). Adding MSNBURR-IIa distribution to MultiBUGS. *Jurnal Aplikasi
#' Statistika & Komputasi Statistik*, 17(2), 109--130.
#' \doi{10.34123/jurnalasks.v17i2.804}
#'
#' Choir, A. S. (2020). *Distribusi neo-normal baru dan karakteristiknya*
#' \[Doctoral dissertation\]. Institut Teknologi Sepuluh Nopember.
#'
#' Mächler, M. (2012). *Accurately computing* \eqn{\log(1 - \exp(-|a|))}.
#' Vignette of the \pkg{Rmpfr} package.
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

    # [REVISI-01] nimStop() diganti return(NaN); ditambah cek NaN.
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

    # [REVISI-01] nimStop() diganti return(NaN). Cek NaN pada q dan parameter
    # [REVISI-01] sebelumnya TIDAK ADA di fungsi ini (beda dengan pmsnburr);
    # [REVISI-01] sekarang disamakan.
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

    # [REVISI-01] nimStop() diganti return(NaN), termasuk untuk p di luar
    # [REVISI-01] [0, 1] dan log(p) > 0.
    # x != x is TRUE only for NaN in compiled code. is.na() cannot be used
    # here: NIMBLE's automatic differentiation (buildDerivs) does not
    # support it for double arguments.
    if (p != p | mu != mu | sigma != sigma | alpha != alpha) return(NaN)
    if (abs(mu) == Inf | sigma == Inf | alpha == Inf) return(NaN)
    if (sigma <= 0 | alpha <= 0) return(NaN)

    # Convert the input, whatever its tail and scale, to log S (upper tail),
    # returning early at the ends of the support.
    log_sf <- 0
    if (log.p != 0) {
      if (p > 0) return(NaN)
      if (lower.tail != 0) {
        if (p == -Inf) return(-Inf)
        log_sf <- log1mexp_nimble(-p)
      } else {
        if (p == -Inf) return(Inf)
        log_sf <- p
      }
    } else {
      if (p < 0 | p > 1) return(NaN)
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

    # [REVISI-01] Parameter tidak valid: nimStop() diganti return(NaN).
    if (is.na(mu) | is.na(sigma) | is.na(alpha)) return(NaN)
    if (abs(mu) == Inf | sigma == Inf | alpha == Inf) return(NaN)
    if (sigma <= 0 | alpha <= 0) return(NaN)

    omega <- exp(log_omega_msnburr(alpha))

    # Inverse-transform sampling with S = 1 - U. log1p(-u) keeps precision
    # for small u; runif() never returns 0 or 1, so the result is finite.
    u <- runif(1, 0, 1)
    a <- -log1p(-u) / alpha
    log_term <- a + log1mexp_nimble(a)

    return(mu + (sigma / omega) * (log(alpha) + log_term))
  }
)
