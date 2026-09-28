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
