#' Construct a classed `bacipair_abstention` condition
#'
#' `bacipair`'s estimators decline to produce a result when a *statistical*
#' assumption fails on the input panel (too few comparison donors, too few
#' pre-treatment years) -- as distinct from a malformed call, which stays a
#' plain validation error. Those declines are signalled as a classed
#' condition rather than a bare `stop()`, so that cross-member integration
#' code can recognise "assumptions not met" without depending on this
#' package's namespace, per the orchestra's refusal/abstention contract
#' (`ORCHESTRA_dev/integration/refusal_contract.R`): any class ending in
#' `_refusal` or `_abstention` is a decline, and `bacipair` produces a typed
#' object recording what was checked rather than refusing outright, so its
#' class ends in `_abstention`.
#'
#' @param message Character scalar. Human-readable diagnostic, matching the
#'   message the equivalent `stop()` used to carry.
#' @param reason Character scalar. Short machine-readable reason code
#'   (e.g. `"donor_pool_too_small"`, `"insufficient_pre_years"`).
#' @param data Named list of diagnostic values relevant to the decline
#'   (e.g. the observed donor count and the required minimum).
#' @return A condition object of class `c("bacipair_abstention",
#'   "orchestra_refusal", "error", "condition")`, produced by
#'   `structure()` around a standard R condition list (`message`, `call`).
#' @keywords internal
new_bacipair_abstention <- function(message, reason, data = list()) {
  structure(
    class = c("bacipair_abstention", "orchestra_refusal",
             "error", "condition"),
    list(
      message = message,
      call    = sys.call(-1L),
      reason  = reason,
      data    = data
    )
  )
}

#' Signal a `bacipair_abstention`
#'
#' Raises the condition built by [new_bacipair_abstention()]. `bacipair`'s
#' abstentions are unrecoverable for the calling estimator -- there is no
#' result to return -- so this is a thin wrapper around `stop()` on the
#' classed condition rather than a full `rlang::abort()`-style restart
#' protocol; the package carries no dependency on `rlang`.
#'
#' @inheritParams new_bacipair_abstention
#' @return Never returns; raises the condition.
#' @keywords internal
abstain <- function(message, reason, data = list()) {
  stop(new_bacipair_abstention(message, reason, data))
}

#' Is `x` a `bacipair` abstention?
#'
#' Producer-side predicate; complements the orchestra-wide,
#' producer-agnostic `is_orchestra_decline()` in
#' `ORCHESTRA_dev/integration/refusal_contract.R`, which recognises this
#' class (and every other member's) by the `_refusal`/`_abstention` naming
#' convention without depending on `bacipair`'s namespace.
#'
#' @param x Any object, typically a caught condition.
#' @return A length-1 logical.
#' @examples
#' sim <- simulate_staggered_panel(seed = 20260520L, n_pairs = 4L,
#'                                  n_proponents = 2L)
#' cond <- tryCatch(
#'   synth_control_baci(sim$panel),
#'   bacipair_abstention = function(e) e
#' )
#' is_bacipair_abstention(cond)
#' @export
is_bacipair_abstention <- function(x) {
  inherits(x, "bacipair_abstention")
}
