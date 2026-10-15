# bacipair <img src="man/figures/logo.png" align="right" height="139" alt="bacipair logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/max578/bacipair/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/max578/bacipair/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/max578/bacipair/graph/badge.svg)](https://app.codecov.io/gh/max578/bacipair)
[![r-universe](https://max578.r-universe.dev/badges/bacipair)](https://max578.r-universe.dev/bacipair)
<!-- badges: end -->

`bacipair` estimates a treatment effect on a paired Before-After-Control-Impact
(BACI) panel with staggered treatment timing -- paired difference-in-differences
with cluster-robust standard errors, a simplex-weighted synthetic-control
comparator, a proponent-level cluster bootstrap, and a Goodman-Bacon
diagnostic -- and grades every estimator against a known-oracle simulator
before it is trusted on a real register.

## Related packages

- [did](https://CRAN.R-project.org/package=did) estimates group-time average
  treatment effects in difference-in-differences designs with more than two
  periods and staggered treatment timing (Callaway and Sant'Anna, 2021).
- [Synth](https://CRAN.R-project.org/package=Synth) builds a synthetic
  control unit from a weighted combination of comparison units for a
  comparative case study of one treated unit.
- [bacondecomp](https://CRAN.R-project.org/package=bacondecomp) decomposes a
  two-way fixed-effects difference-in-differences estimate under staggered
  timing into its 2x2 comparisons (Goodman-Bacon).

bacipair is built for paired BACI panels, where each impact site has its own
control. It puts paired difference-in-differences and synthetic control side
by side on the same panel, uses bacondecomp for its staggered-timing
diagnostic, and checks each estimator against a simulator with a known effect
before it is used on real data.

## Installation

```r
# install.packages("remotes")
remotes::install_github("max578/bacipair")
```

## Example

```r
library(bacipair)

sim <- simulate_staggered_panel(seed = 20260520L, n_pairs = 80L)
fit <- synth_control_baci(sim$panel)

cat(sprintf("ATT_hat = %.3f pp  (oracle = %.3f pp)\n",
            fit$panel_att_hat, sim$oracle_panel_att))
```

See `vignette("getting-started", package = "bacipair")` for the estimand and
estimator mathematics, the paired-DiD comparator, the cluster bootstrap, and
the Goodman-Bacon staggered-DiD diagnostic.

## Contributing

Bug reports and suggestions are welcome through the
[issue tracker](https://github.com/max578/bacipair/issues).

## Citation

If you use this package, please cite both the package and the source
methodology paper:

> Macintosh, A., et al. (2024). Australian human-induced native forest
> regeneration carbon offset projects have limited impact on changes in
> woody vegetation cover and carbon removals. *Communications Earth &
> Environment*, 5, Article 149. <https://doi.org/10.1038/s43247-024-01313-x>

## Licence

MIT
