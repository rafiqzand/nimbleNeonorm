# Registration of the distributions with NIMBLE's user-defined distribution
# registry, and the matching cleanup.


# Names of the distributions currently registered by this package, so that
# re-registration and cleanup only touch what we added.
.nimble_state <- new.env(parent = emptyenv())
.nimble_state$registered <- character(0)

# Density name -> display name, used in messages and return values.
.formal_names <- c(dmsnburr = "MSNBurr", dmsnburr2a = "MSNBurr-IIa")


# NIMBLE registration metadata for the requested distributions.
.make_distribution_registration <- function(distributions) {
  # Both distributions are continuous on the whole real line and provide
  # p and q functions, which lets NIMBLE handle truncation, e.g. T(y, 0, ).
  spec <- function(name) {
    list(
      BUGSdist = paste0(name, "(mu, sigma, alpha)"),
      Rdist = paste0(name, "(mu, sigma, alpha)"),
      discrete = FALSE,
      pqAvail = TRUE,
      range = c(-Inf, Inf)
    )
  }

  dist_list <- list()
  if ("msnburr" %in% distributions) {
    dist_list$dmsnburr <- spec("dmsnburr")
  }
  if ("msnburr2a" %in% distributions) {
    dist_list$dmsnburr2a <- spec("dmsnburr2a")
  }
  dist_list
}


#' Register the neo-normal distributions with NIMBLE
#'
#' Makes the MSNBurr and/or MSNBurr-IIa distributions available in NIMBLE
#' model code, e.g. `y[i] ~ dmsnburr(mu, sigma, alpha)`. Call this before
#' [nimble::nimbleModel()].
#'
#' `register_msnburr()` and `register_msnburr2a()` are shortcuts for
#' registering a single distribution.
#'
#' Calling `register_neonorm()` again first removes the distributions it
#' registered previously, so it is safe to call repeatedly (for example after
#' reloading the package during development). Use [cleanup_neonorm()] to
#' remove them when you are done.
#'
#' @param distributions Character vector of distributions to register:
#'   `"MSNBurr"`, `"MSNBurr-IIa"` (or `"MSNBurr2a"`), or `"all"`.
#'   Matching is case-insensitive.
#' @param verbose Logical; print a message listing what was registered?
#'
#' @return Invisibly, a character vector with the display names of the
#'   registered distributions, e.g. `c("MSNBurr", "MSNBurr-IIa")`.
#'
#' @seealso [MSNBurr], [MSNBurr2a], [cleanup_neonorm()],
#'   [nimble::registerDistributions()]
#'
#' @examples
#' \donttest{
#' register_neonorm()
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
#' model$calculate()
#'
#' cleanup_neonorm()
#' }
#'
#' @export
register_neonorm <- function(distributions = "all", verbose = TRUE) {
  distributions <- tolower(distributions)

  allowed <- c("all", "msnburr", "msnburr-iia", "msnburr2a")
  invalid <- setdiff(distributions, allowed)
  if (length(invalid) > 0) {
    stop(
      "Unknown distribution(s): ", paste(invalid, collapse = ", "), ". ",
      "Allowed values are: ", paste(allowed, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if ("all" %in% distributions) {
    distributions <- c("msnburr", "msnburr2a")
  } else {
    distributions[distributions == "msnburr-iia"] <- "msnburr2a"
    distributions <- unique(distributions)
  }

  # NIMBLE refuses to register a name twice, so drop our earlier
  # registration first.
  if (length(.nimble_state$registered) > 0) {
    nimble::deregisterDistributions(.nimble_state$registered)
    .nimble_state$registered <- character(0)
  }

  dist_list <- .make_distribution_registration(distributions)

  # The nimbleFunctions named in `Rdist` live in the package namespace.
  nimble::registerDistributions(
    dist_list,
    userEnv = asNamespace("nimbleNeonorm"),
    verbose = FALSE
  )
  .nimble_state$registered <- names(dist_list)

  registered_names <- unname(.formal_names[.nimble_state$registered])
  if (verbose) {
    message(
      "Successfully registered: ",
      paste(registered_names, collapse = ", ")
    )
  }

  invisible(registered_names)
}


#' @rdname register_neonorm
#' @export
register_msnburr <- function(verbose = TRUE) {
  register_neonorm(distributions = "msnburr", verbose = verbose)
}


#' @rdname register_neonorm
#' @export
register_msnburr2a <- function(verbose = TRUE) {
  register_neonorm(distributions = "msnburr2a", verbose = verbose)
}


#' Remove the neo-normal distributions from NIMBLE
#'
#' Deregisters every distribution added by [register_neonorm()] from NIMBLE's
#' registry of user-defined distributions. Distributions registered by other
#' packages or by the user are left untouched. Calling it when nothing is
#' registered does nothing.
#'
#' @param verbose Logical; print a message listing what was removed?
#'
#' @return Invisibly, `NULL`.
#'
#' @seealso [register_neonorm()], [nimble::deregisterDistributions()]
#'
#' @examples
#' \donttest{
#' register_neonorm(verbose = FALSE)
#' cleanup_neonorm()
#' }
#'
#' @export
cleanup_neonorm <- function(verbose = TRUE) {
  registered <- .nimble_state$registered

  if (length(registered) == 0) {
    if (verbose) {
      message("No neo-normal distributions are currently registered.")
    }
    return(invisible(NULL))
  }

  nimble::deregisterDistributions(registered)
  .nimble_state$registered <- character(0)

  if (verbose) {
    message(
      "Deregistered: ",
      paste(unname(.formal_names[registered]), collapse = ", ")
    )
  }

  invisible(NULL)
}
