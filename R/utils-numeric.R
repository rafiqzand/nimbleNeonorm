# Numerically stable building blocks shared by the MSNBurr and MSNBurr-IIa
# nimbleFunctions.
#
# All helpers are exported (but kept out of the documentation index through
# `@keywords internal`). This is required: when NIMBLE compiles a model, it
# resolves nimbleFunctions called from other nimbleFunctions by name on the
# search path, so helpers that live only in the package namespace cannot be
# found and compilation fails with "Problem with type of arg2 in
# sizeBinaryCwise".

#' @importFrom nimble nimbleFunction nimStop
NULL


# ---- log1pexp_nimble --------------------------------------------------------

#' Numerically stable `log(1 + exp(x))`
#'
#' Evaluates \eqn{\log(1 + e^x)} without overflow for large `x` and without
#' loss of precision for very negative `x`, using the piecewise scheme of
#' Mächler (2012).
#'
#' | Region               | Expression used    |
#' |----------------------|--------------------|
#' | `x <= -37`           | `exp(x)`           |
#' | `-37 < x <= 18`      | `log1p(exp(x))`    |
#' | `18 < x <= 33.3`     | `x + exp(-x)`      |
#' | `x > 33.3`           | `x`                |
#'
#' @param x Numeric scalar.
#'
#' @return Numeric scalar, \eqn{\log(1 + e^x)}.
#'
#' @references
#' Mächler, M. (2012). *Accurately computing* \eqn{\log(1 - \exp(-|a|))}.
#' Vignette of the \pkg{Rmpfr} package.
#'
#' @keywords internal
#' @export
log1pexp_nimble <- nimble::nimbleFunction(
  run = function(x = double(0)) {
    returnType(double(0))

    # At these cut-offs each approximation equals log(1 + exp(x)) to full
    # double precision.
    if (x <= -37.0) {
      return(exp(x))
    }
    if (x <= 18.0) {
      return(log1p(exp(x)))
    }
    if (x <= 33.3) {
      return(x + exp(-x))
    }
    return(x)
  },
  buildDerivs = TRUE
)


# ---- expm1_nimble -----------------------------------------------------------

#' Numerically stable `exp(x) - 1`
#'
#' Evaluates \eqn{e^x - 1} without the cancellation that `exp(x) - 1`
#' suffers for `x` near zero. For \eqn{|x| \le 0.7} a 16th-order Taylor
#' polynomial (in Horner form) is used; otherwise `exp(x) - 1` is accurate.
#'
#' The NIMBLE DSL has no `expm1()`, hence this helper.
#'
#' @param x Numeric scalar.
#'
#' @return Numeric scalar, \eqn{e^x - 1}.
#'
#' @keywords internal
#' @export
expm1_nimble <- nimble::nimbleFunction(
  run = function(x = double(0)) {
    returnType(double(0))

    if (abs(x) <= 0.7) {
      # sum_{k = 1}^{16} x^k / k!, in Horner form.
      return(x * (1 + x * (1 / 2 + x * (1 / 6 + x * (1 / 24 +
        x * (1 / 120 + x * (1 / 720 + x * (1 / 5040 + x * (1 / 40320 +
        x * (1 / 362880 + x * (1 / 3628800 + x * (1 / 39916800 +
        x * (1 / 479001600 + x * (1 / 6227020800 + x * (1 / 87178291200 +
        x * (1 / 1307674368000 + x * (1 / 20922789888000)))))))))))))))))
    }
    return(exp(x) - 1)
  },
  buildDerivs = TRUE
)


# ---- log1mexp_nimble --------------------------------------------------------

#' Numerically stable `log(1 - exp(-x))`
#'
#' Evaluates \eqn{\log(1 - e^{-x})} for \eqn{x \ge 0} following Mächler
#' (2012): `log(-expm1(-x))` for \eqn{x \le \log 2} and `log1p(-exp(-x))`
#' otherwise.
#'
#' Its main use is turning a log-probability into the log of its
#' complement: if \eqn{\ell = \log F \le 0}, then
#' \eqn{\log(1 - F) = } `log1mexp_nimble(-l)`.
#'
#' @param x Non-negative numeric scalar.
#'
#' @return Numeric scalar, \eqn{\log(1 - e^{-x})}; `-Inf` when `x = 0`.
#'   Stops with an error when `x < 0`, where the result is not real.
#'
#' @references
#' Mächler, M. (2012). *Accurately computing* \eqn{\log(1 - \exp(-|a|))}.
#' Vignette of the \pkg{Rmpfr} package.
#'
#' @keywords internal
#' @export
log1mexp_nimble <- nimble::nimbleFunction(
  run = function(x = double(0)) {
    returnType(double(0))

    if (x < 0) {
      nimStop("log1mexp_nimble: x must be >= 0")
    }
    if (x == 0) {
      return(-Inf)
    }
    if (x <= 0.6931471805599453) { # log(2)
      return(log(-expm1_nimble(-x)))
    }
    return(log1p(-exp(-x)))
  },
  buildDerivs = TRUE
)


# ---- log_omega_msnburr ------------------------------------------------------

#' Log of the MSNBurr scaling constant
#'
#' Computes \eqn{\log\omega} where
#' \deqn{\omega = \frac{(1 + 1/\alpha)^{\alpha + 1}}{\sqrt{2\pi}},}
#' the constant shared by the MSNBurr and MSNBurr-IIa distributions. Working
#' on the log scale avoids overflow of \eqn{(1 + 1/\alpha)^{\alpha + 1}}, and
#' `log1p()` keeps precision when \eqn{1/\alpha} is small.
#'
#' @param alpha Positive numeric scalar, the shape parameter.
#'
#' @return Numeric scalar, \eqn{\log\omega}.
#'
#' @keywords internal
#' @export
log_omega_msnburr <- nimble::nimbleFunction(
  run = function(alpha = double(0)) {
    returnType(double(0))

    # 0.9189385332046727 = 0.5 * log(2 * pi)
    return((alpha + 1) * log1p(1 / alpha) - 0.9189385332046727)
  },
  buildDerivs = TRUE
)
