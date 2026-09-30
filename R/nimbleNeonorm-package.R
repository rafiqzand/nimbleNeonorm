# [REVISI-05] Dokumentasi tingkat paket: ditambah ringkasan alur kerja
# [REVISI-05] (fungsi per langkah Bayesian workflow) dan @references inti.
# [REVISI-05] Bagian kode (import & globalVariables) TIDAK berubah.
#' nimbleNeonorm: Neo-Normal Distributions for 'NIMBLE'
#'
#' Provides the MSNBurr and MSNBurr-IIa neo-normal distributions as
#' user-defined distributions for 'NIMBLE', together with helpers that follow
#' the Bayesian workflow of Gelman et al. (2020).
#'
#' @section Workflow:
#' 1. **Distributions**: [MSNBurr], [MSNBurr2a]; make them available in
#'    model code with [register_neonorm()].
#' 2. **Prior predictive check**: [prior_predict_neonorm()].
#' 3. **Fitting** with several chains from dispersed starting values:
#'    [runmcmc_neonorm()], [configure_mcmc_neonorm()],
#'    [prior_inits_neonorm()].
#' 4. **Convergence diagnostics** (\eqn{\hat{R}}, bulk and tail ESS):
#'    [summarise_neonorm()].
#' 5. **Posterior predictive check**: [posterior_predict_neonorm()].
#' 6. **Model comparison** with WAIC: `runmcmc_neonorm(..., WAIC = TRUE)`.
#'
#' See `vignette("bayesian-workflow", package = "nimbleNeonorm")`.
#'
#' @references
#' Choir, A. S. (2020). *Distribusi neo-normal baru dan karakteristiknya*
#' \[Doctoral dissertation\]. Institut Teknologi Sepuluh Nopember.
#'
#' de Valpine, P., Turek, D., Paciorek, C. J., Anderson-Bergman, C.,
#' Temple Lang, D., & Bodik, R. (2017). Programming with models: Writing
#' statistical algorithms for general model structures with NIMBLE.
#' *Journal of Computational and Graphical Statistics*, 26(2), 403--413.
#' \doi{10.1080/10618600.2016.1172487}
#'
#' Gelman, A., Vehtari, A., Simpson, D., Margossian, C. C., Carpenter, B.,
#' Yao, Y., Kennedy, L., Gabry, J., Bürkner, P.-C., & Modrák, M. (2020).
#' *Bayesian workflow*. arXiv:2011.01808.
#' \doi{10.48550/arXiv.2011.01808}
#'
#' Ramadani, E. P., Choir, A. S., Pravitasari, A. A., & Paraguison, J.
#' (2025). Adding MSNBURR-IIa distribution to MultiBUGS. *Jurnal Aplikasi
#' Statistika & Komputasi Statistik*, 17(2), 109--130.
#' \doi{10.34123/jurnalasks.v17i2.804}
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import methods
#' @importFrom stats runif
## usethis namespace: end
NULL

# `returnType()` only exists inside the NIMBLE DSL, so R CMD check cannot see
# a definition for it in the nimbleFunction bodies.
utils::globalVariables("returnType")
