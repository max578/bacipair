# bacipair

`bacipair` estimates a treatment effect on a paired Before-After-Control-Impact
(BACI) panel with staggered treatment timing -- paired difference-in-differences
with cluster-robust standard errors, a simplex-weighted synthetic-control
comparator, a proponent-level cluster bootstrap, and a Goodman-Bacon
diagnostic -- and grades every estimator against a known-oracle simulator
before it is trusted on a real register.

## Installation

```r
# install.packages("remotes")
remotes::install_local("path/to/bacipair")
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

## Citation

If you use this package, please cite both the package and the source
methodology paper:

> Macintosh, A., et al. (2024). Australian human-induced native forest
> regeneration carbon offset projects have limited impact on changes in
> woody vegetation cover and carbon removals. *Communications Earth &
> Environment*, 5, Article 149. <https://doi.org/10.1038/s43247-024-01313-x>

## Licence

MIT
