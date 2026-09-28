# Package hooks. Nothing is registered with NIMBLE on load: attaching the
# package must not change NIMBLE's global distribution registry behind the
# user's back, so registration is left to register_neonorm().

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "nimbleNeonorm ", utils::packageVersion(pkgname),
    ": call register_neonorm() before building models with ",
    "dmsnburr() or dmsnburr2a()."
  )
}
