## Cached small simulator output reused across tests; cuts redundant
## simulator work in `R CMD check`. The seed is locked to keep tests
## deterministic. We pass an explicit small `n_proponents` so callers
## can stay below the default of 50 without tripping the n_proponents
## <= n_pairs guard.

cached_sim <- function(seed = 20260520L, n_pairs = 30L,
                       n_proponents = 8L) {
  key <- sprintf("sim_%d_%d_%d", seed, n_pairs, n_proponents)
  if (!exists(key, envir = .GlobalEnv, inherits = FALSE)) {
    sim <- bacipair::simulate_staggered_panel(
      seed = seed, n_pairs = n_pairs,
      n_proponents = n_proponents
    )
    assign(key, sim, envir = .GlobalEnv)
  }
  get(key, envir = .GlobalEnv)
}
