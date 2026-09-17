# InvaCost Damage Cost:Management Ratio

## Project shape

This repository accompanies the 2024 *Ecological Economics* paper on national invasive-species damage costs relative to management expenditure. It is an R analysis repository, not an R package: there is no dependency lockfile, build system, test suite, or lint configuration.

The analysis is driven by `scripts/ntlScaleCostDiffInvaCostGithub2024.R`:

1. It loads the `invacost` package dataset, retains high-reliability records, and bootstraps country-level damage and management summaries over 2000--2020. It derives the damage:management ratio, its rate of change, and damage/management proportions.
2. It joins those outcomes to `data/continent.countryINVACOST.csv` and the national socio-economic datasets. Country joins use `country` initially and then `cntry.code`; preserve these identifiers and their spelling when updating data.
3. It derives per-capita imports and scientific output using `pop2021.csv`, selects the most recent FAO agriculture value by country, imputes the eight predictor variables with `mice` (eight PMM imputations, seed 101), and applies the existing log/logit/scale transformations.
4. It analyses the ratio, proportion damaged, and ratio rate of change through phased boosted regression trees and exhaustive `lmer` model sets with regional random intercepts. The ratio is additionally assessed with spatial GLS using country centroids.

`scripts/new_lmer_AIC_tables3.r` supplies AICc/BIC/deviance-explained utilities for model lists. `scripts/r.squared.R` defines the `r.squared()` S3 generic and methods used to obtain marginal and conditional R-squared values for `lm`, `merMod`, and `lme` fits.

## Running and checking

Required packages are listed in the README: `invacost`, `lme4`, `dismo`, `gbm`, `boot`, `VIM`, `mice`, `performance`, `sjPlot`, `rworldmap`, `SpatialEpi`, `nlme`, and `rcompanion`. The script also uses spatial functionality associated with the commented `rgeos` import, so run it only in an R environment where its spatial dependencies are available.

There is no supported one-command full analysis run in the committed state. The main script reads helpers and all CSVs by bare filename, while they reside in separate `scripts/` and `data/` directories; it also changes the working directory late in the script to the author's machine-specific output path. Before a full rerun, make those paths explicit or stage the required files in the chosen working directory, and replace the hard-coded output directory. Do not commit generated CSV outputs unless intentionally adding a reproducible result artefact.

Use parse checks for a changed R file:

```sh
Rscript -e 'parse("scripts/ntlScaleCostDiffInvaCostGithub2024.R")'
Rscript -e 'parse("scripts/r.squared.R")'
```

There is no automated test suite or lint command, and consequently no single-test command. For a focused behavioural check of a helper after changing it, source that helper in `Rscript -e` and exercise the affected function with a minimal model object.

## Conventions that affect results

- Treat `ntlScaleCostDiffInvaCostGithub2024.R` as an ordered analysis: later sections depend on objects built earlier in the same R session. Do not extract or rerun a later block without preparing its upstream objects.
- Preserve the defined transformations and names (`gdp`, `cpi`, `ghsi`, `igs`, `agrL`, `VAag`, `govexpedu`, `stja`; responses `ratio`, `pdam`, and `r`). Model formula construction and BRT column-index selection depend on this column order.
- The grouped-region variable `reg2` deliberately combines NAM/CAR, EUR/ME, and ASIA/OC before mixed modelling to increase category sample sizes. Keep that recoding consistent across response datasets.
- Bootstrap and BRT output filenames are reused by distinct analysis blocks (for example, `BRT.boot.*.csv`). Direct reruns can overwrite earlier results; use an explicit output directory and distinct names when retaining outputs from multiple blocks.
- The main analysis contains computationally intensive resampling (country bootstrap iterations, `mice` with 500 iterations, and BRT bootstrap loops). Preserve seeds and iteration settings when reproducing published results; make deliberate, documented changes when reducing them for exploratory work.
