#' Simulate a staggered-treatment HIR-style BACI panel
#'
#' Generates a synthetic Before-After-Control-Impact (BACI) panel with
#' staggered treatment timing, two arms (`"credited"` and `"comparison"`),
#' per-pair static covariates (baseline cover, elevation, rainfall
#' climatology), per-year shared rainfall anomaly, and proponent
#' assignments drawn from a Zipf distribution. The data-generating process
#' is structured so that the panel-average ATT is recoverable from the
#' comparison arm; paired-DiD with the matched 1:1 comparator carries
#' residual bias when elevation / rainfall climatology covary with
#' treatment timing; and the donor-pool synthetic-control estimator should
#' recover the oracle within the kill-gate tolerance documented in
#' [synth_control_baci()].
#'
#' All values, distributions, and seeds are explicit. There is no hidden
#' RNG. Inputs are validated by `stopifnot()` guards.
#'
#' @param seed Integer scalar. The simulation RNG seed (default
#'   `20260520L`).
#' @param n_pairs Integer scalar (>= 2). Number of treated/comparison
#'   pairs (default 200).
#' @param t_min,t_max Integer scalars. Inclusive panel year range
#'   (defaults 2010 / 2025).
#' @param treatment_year_min,treatment_year_max Integer scalars. Range
#'   of staggered treatment years (defaults 2013 / 2022). Must lie inside
#'   `[t_min, t_max]`.
#' @param additionality_share Numeric in `[0, 1]`. Fraction of pairs with
#'   non-zero true additionality (default 0.18).
#' @param additionality_mean_pp,additionality_sd_pp Numeric. Mean / SD of
#'   the (untruncated) Normal asymptotic-additionality distribution among
#'   pairs with `has_real_additionality == TRUE`, in percentage points of
#'   woody cover (defaults 4.5 / 1.5). The draw is signed, not
#'   non-negative: a pair flagged as genuinely additional can realise a
#'   negative `true_additionality_max` when `additionality_mean_pp` is
#'   small relative to `additionality_sd_pp`, and that signed value flows
#'   unaltered into `oracle_panel_att`.
#' @param additionality_dynamics Numeric vector of non-negative weights
#'   indexed by years-since-treatment (1, 2, ...). The realised
#'   individual treatment effect in year `y` is
#'   `true_additionality_max * additionality_dynamics[y - treatment_year + 1]`.
#'   Saturates at the last element.
#' @param rainfall_sd_pp,noise_sd_pp Numeric. Rainfall-anomaly and
#'   idiosyncratic noise SDs in percentage points of woody cover.
#' @param site_baseline_mean_pp,site_baseline_sd_pp Numeric. Mean and SD
#'   of the site-baseline draw in percentage points of woody cover.
#' @param site_pair_sd_pp Numeric. Within-pair shift SD in percentage
#'   points of woody cover.
#' @param secular_trend_pp_per_year Numeric. Common secular trend in
#'   percentage points of woody cover per year (default 0.4).
#' @param elevation_mean_m,elevation_sd_m Numeric. Per-pair elevation
#'   distribution (defaults 350 / 200 m).
#' @param elevation_trend_coupling Numeric. Coupling between elevation
#'   and the secular trend (default `0.0015` pp/year per metre).
#' @param rainfall_climatology_mean_pp,rainfall_climatology_sd_pp Numeric.
#'   Per-pair rainfall-climatology offset distribution.
#' @param elevation_treatment_year_coupling Numeric. Strength of the
#'   selection-on-observables coupling between elevation and treatment
#'   year (default -2.5; higher elevation means earlier registration).
#' @param n_proponents Integer (1 <= n_proponents <= n_pairs). Number of
#'   distinct proponents to assign (default 50).
#' @param proponent_zipf_s Numeric > 0. Zipf shape parameter for proponent
#'   cluster sizes (default 1.2).
#'
#' @return A list with elements:
#'   \describe{
#'     \item{`panel`}{A `data.table` with columns including `pair_id`,
#'       `proponent_id`, `arm`, `year`, `treatment_year`,
#'       `post_treatment`, `treated_arm`, `post_idx`, `baseline_cover`,
#'       `elevation_m`, `rainfall_climatology`, `rainfall_sensitivity`,
#'       `site_pair_shift`, `rainfall_anomaly`,
#'       `true_additionality_max`, `has_real_additionality`, `ite`,
#'       `noise`, `woody_cover`.}
#'     \item{`pair_table`}{A `data.table` of per-pair covariates.}
#'     \item{`rainfall`}{A `data.table` of per-year shared rainfall
#'       anomaly.}
#'     \item{`oracle_panel_att`}{Numeric scalar. The oracle panel-average
#'       ATT (mean ITE over treated post-treatment observations).}
#'     \item{`oracle_unit_att`}{A `data.table` of per-unit oracle ATT
#'       (`pair_id`, `unit_att`).}
#'     \item{`parameters`}{Named list of every parameter used to generate
#'       the panel, plus `generated_at` ISO-8601 timestamp.}
#'   }
#'
#' @examples
#' sim <- simulate_staggered_panel(seed = 20260520L, n_pairs = 50L)
#' sim$oracle_panel_att
#' head(sim$panel)
#' @export
simulate_staggered_panel <- function(seed = 20260520L,
                                     n_pairs = 200L,
                                     t_min = 2010L,
                                     t_max = 2025L,
                                     treatment_year_min = 2013L,
                                     treatment_year_max = 2022L,
                                     additionality_share = 0.18,
                                     additionality_mean_pp = 4.5,
                                     additionality_sd_pp = 1.5,
                                     additionality_dynamics =
                                       c(0.20, 0.55, 0.85, 1.00, 1.00,
                                         1.00, 1.00, 1.00, 1.00, 1.00),
                                     rainfall_sd_pp = 1.4,
                                     site_baseline_mean_pp = 22.0,
                                     site_baseline_sd_pp = 4.0,
                                     site_pair_sd_pp = 1.0,
                                     noise_sd_pp = 0.8,
                                     secular_trend_pp_per_year = 0.4,
                                     elevation_mean_m = 350,
                                     elevation_sd_m = 200,
                                     elevation_trend_coupling = 0.0015,
                                     rainfall_climatology_mean_pp = 0.0,
                                     rainfall_climatology_sd_pp = 0.6,
                                     elevation_treatment_year_coupling = -2.5,
                                     n_proponents = 50L,
                                     proponent_zipf_s = 1.2) {
  stopifnot(
    is.numeric(seed), length(seed) == 1L, !is.na(seed),
    n_pairs > 1L,
    t_max > t_min,
    treatment_year_min >= t_min, treatment_year_max <= t_max,
    treatment_year_max > treatment_year_min,
    additionality_share >= 0, additionality_share <= 1,
    n_proponents >= 1L, n_proponents <= n_pairs,
    proponent_zipf_s > 0,
    length(additionality_dynamics) >= 1L,
    all(additionality_dynamics >= 0)
  )

  set.seed(seed)

  pair_id_vec <- seq_len(n_pairs)

  baseline_cover_vec <- stats::rnorm(
    n_pairs, mean = site_baseline_mean_pp, sd = site_baseline_sd_pp
  )
  elevation_vec <- pmax(0, stats::rnorm(
    n_pairs, mean = elevation_mean_m, sd = elevation_sd_m
  ))
  rainfall_climatology_vec <- stats::rnorm(
    n_pairs,
    mean = rainfall_climatology_mean_pp,
    sd   = rainfall_climatology_sd_pp
  )
  rainfall_sensitivity_vec <- stats::rnorm(n_pairs, mean = 1.0, sd = 0.25)

  treatment_window <- seq.int(treatment_year_min, treatment_year_max)
  elevation_z <- as.numeric(scale(elevation_vec))
  base_logits <- matrix(0, nrow = n_pairs, ncol = length(treatment_window))
  for (k in seq_along(treatment_window)) {
    base_logits[, k] <- elevation_treatment_year_coupling *
      elevation_z * (k - mean(seq_along(treatment_window))) /
      length(treatment_window)
  }
  base_logits <- base_logits + matrix(
    stats::rnorm(n_pairs * length(treatment_window), sd = 1.0),
    nrow = n_pairs
  )
  treatment_year_vec <- treatment_window[apply(base_logits, 1L, which.max)]

  has_real_add_vec <- stats::rbinom(
    n_pairs, size = 1L, prob = additionality_share
  ) == 1L
  true_additionality_max_vec <- ifelse(
    has_real_add_vec,
    stats::rnorm(n_pairs, mean = additionality_mean_pp,
                 sd = additionality_sd_pp),
    0
  )

  ranks <- seq_len(n_proponents)
  weights <- 1 / (ranks ^ proponent_zipf_s)
  weights <- weights / sum(weights)
  proponent_vec <- sample.int(n_proponents, size = n_pairs,
                              replace = TRUE, prob = weights)

  pair_dt <- data.table::data.table(
    pair_id              = pair_id_vec,
    proponent_id         = proponent_vec,
    baseline_cover       = baseline_cover_vec,
    elevation_m          = elevation_vec,
    rainfall_climatology = rainfall_climatology_vec,
    rainfall_sensitivity = rainfall_sensitivity_vec,
    treatment_year       = treatment_year_vec,
    true_additionality_max = true_additionality_max_vec,
    has_real_additionality = has_real_add_vec
  )

  years <- seq.int(t_min, t_max)
  rainfall_anomaly_vec <- stats::rnorm(
    length(years), mean = 0, sd = rainfall_sd_pp
  )
  rainfall_dt <- data.table::data.table(
    year             = years,
    rainfall_anomaly = rainfall_anomaly_vec
  )

  arms <- c("credited", "comparison")
  panel <- data.table::CJ(pair_id = pair_id_vec, year = years, arm = arms,
                          sorted = FALSE)
  panel <- merge(panel, pair_dt,    by = "pair_id", all.x = TRUE)
  panel <- merge(panel, rainfall_dt, by = "year",   all.x = TRUE)

  pair_shift_vec <- stats::rnorm(n_pairs, mean = 0, sd = site_pair_sd_pp)
  panel[, "site_pair_shift" := pair_shift_vec[get("pair_id")]]
  panel[, "noise" := stats::rnorm(.N, mean = 0, sd = noise_sd_pp)]

  panel[, "post_treatment" := as.integer(get("year") >= get("treatment_year"))]
  panel[, "treated_arm"    := as.integer(get("arm") == "credited")]
  panel[, "post_idx"       := pmin(
    pmax(get("year") - get("treatment_year") + 1L, 0L),
    length(additionality_dynamics)
  )]
  dynamics_lookup <- c(0, additionality_dynamics)
  panel[, "ite" := get("true_additionality_max") *
                   dynamics_lookup[get("post_idx") + 1L] *
                   get("treated_arm")]

  panel[, "woody_cover" :=
          get("baseline_cover") +
          get("site_pair_shift") +
          (secular_trend_pp_per_year +
             elevation_trend_coupling * get("elevation_m")) *
            (get("year") - t_min) +
          get("rainfall_sensitivity") * get("rainfall_anomaly") +
          get("rainfall_climatology") +
          get("ite") +
          get("noise")]

  oracle_panel_att <- panel[
    get("treated_arm") == 1L & get("post_treatment") == 1L,
    mean(get("ite"))
  ]
  oracle_unit_att <- panel[
    get("treated_arm") == 1L & get("post_treatment") == 1L,
    list(unit_att = mean(get("ite"))),
    by = "pair_id"
  ]

  data.table::setcolorder(panel, c(
    "pair_id", "proponent_id", "arm", "year",
    "treatment_year", "post_treatment", "treated_arm", "post_idx",
    "baseline_cover", "elevation_m", "rainfall_climatology",
    "rainfall_sensitivity",
    "site_pair_shift", "rainfall_anomaly",
    "true_additionality_max", "has_real_additionality", "ite",
    "noise", "woody_cover"
  ))

  list(
    panel             = panel[],
    pair_table        = pair_dt[],
    rainfall          = rainfall_dt[],
    oracle_panel_att  = oracle_panel_att,
    oracle_unit_att   = oracle_unit_att,
    parameters        = list(
      seed                              = seed,
      n_pairs                           = n_pairs,
      t_min                             = t_min,
      t_max                             = t_max,
      treatment_year_min                = treatment_year_min,
      treatment_year_max                = treatment_year_max,
      additionality_share               = additionality_share,
      additionality_mean_pp             = additionality_mean_pp,
      additionality_sd_pp               = additionality_sd_pp,
      additionality_dynamics            = additionality_dynamics,
      rainfall_sd_pp                    = rainfall_sd_pp,
      site_baseline_mean_pp             = site_baseline_mean_pp,
      site_baseline_sd_pp               = site_baseline_sd_pp,
      site_pair_sd_pp                   = site_pair_sd_pp,
      noise_sd_pp                       = noise_sd_pp,
      secular_trend_pp_per_year         = secular_trend_pp_per_year,
      elevation_mean_m                  = elevation_mean_m,
      elevation_sd_m                    = elevation_sd_m,
      elevation_trend_coupling          = elevation_trend_coupling,
      rainfall_climatology_mean_pp      = rainfall_climatology_mean_pp,
      rainfall_climatology_sd_pp        = rainfall_climatology_sd_pp,
      elevation_treatment_year_coupling = elevation_treatment_year_coupling,
      n_proponents                      = n_proponents,
      proponent_zipf_s                  = proponent_zipf_s,
      generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    )
  )
}
