# JRSSA group-testing reproducibility

R code and numerical results for **Designing group testing studies under assay and specimen constraints: a Chlamydia planning study** (revision in progress).

This repository preserves the complete original `inst/` tree from [chikuang/gtDesign](https://github.com/chikuang/gtDesign) and provides executable analyses for the revision, including finite-sample bias, variance and RMSE under a common budget and estimator.

## Contents

| Location | Purpose |
|---|---|
| `inst/examples/` | Original R and R Markdown examples, including the supplied HTML. |
| `inst/paper_reproduce/` | Original paper table/figure scripts and supplied PNG figures. |
| `vendor/gtDesign/` | Unmodified package source used by the calculations. |
| `R/bias_simulation.R` | Joint MLE simulation and allocation construction. |
| `R/profile_audit.R` | Independent likelihood checks by profiling prevalence. |
| `R/validate_theory.R` | Nonidentifiability example and exact binomial bias calculation. |
| `R/revision_checks.R` | Revised rounding, precision bounds, specimen-budget and Stage II checks. |
| `R/summarize_bias.R` | Figures and LaTeX tables generated from the saved R results. |
| `results/bias/` | All 120,000 estimates, generated counts, allocations, Monte Carlo SEs and diagnostics. |
| `results/revision/` | Checks for the revised manuscript. |
| `results/original/` | Results from the preserved original scripts. |
| `figures/`, `docs/` | Generated scientific figures, tables and a technical note. |
| `upstream_manifest.csv` | SHA-256 inventory of the unchanged original `inst/` files. |

The source is `gtDesign` 0.1.1, revision `81341810777b8941aff5dd2937f5dbcb10feedb3`, retrieved on 10 November 2026. Details are in `vendor/gtDesign/SOURCE.txt`. Original filenames refer to earlier manuscript numbering, not the latest revision's table numbers.

## Run

Run from the repository root. The validated environment uses R 4.6.1, CVXR 1.9.3 and CLARABEL. Dependency versions are recorded in `dependencies.csv`; each analysis records session information. Solver/platform differences may cause small numerical differences.

```sh
# Install missing CRAN dependencies and the vendored package into .library/.
Rscript scripts/install_dependencies.R

# 10,000 replicates per cell, four worker processes.
Rscript run.R bias 10000 4

# Revised numerical checks and Stage II figures.
Rscript run.R revision

# Original table/figure scripts and the R example.
Rscript run.R original

# All analyses.
Rscript run.R all 10000 4
```

On Windows the simulation uses one process. Installation does not silently upgrade existing CRAN packages; consult `dependencies.csv` for the validated versions. Runs replace their generated files. Use a separate checkout for pilot runs if you want to preserve the supplied 10,000-replicate results.

Individual original scripts can be run as follows:

```sh
Rscript scripts/run_original.R inst/paper_reproduce/table1.R
Rscript scripts/run_original.R inst/examples/example_calc_Eopt.R
```

Each script runs in a fresh R process. Printed results, data frames, named ggplot objects and session information are saved in `results/original/<script-name>/`. The `.Rmd` example and original HTML are preserved; rendering requires R Markdown and Pandoc.

## Simulation specification

Independently over pool sizes, `Y_x ~ Binomial(n_x, pi_x)`, where

`pi_x = sensitivity - (sensitivity + specificity - 1) * (1 - prevalence)^x`.

| Setting | (Prevalence, sensitivity, specificity) | Maximum size | Budgets |
|---|---|---|---|
| Cost | (0.07, 0.93, 0.96) | 61 | C = 100, 500, 2000; cost(x) = 0.8 + 0.2x |
| Low prevalence | (0.022, 0.92, 0.965) | 15 | N = 500, 2000, 10000 assays |

These settings were specified before examining the Monte Carlo results. Within each cell, rounded D-optimal and uniform four-size designs spend **exactly the same budget**. Assay and specimen totals need not match when monetary cost is binding; both are reported. Uniform designs use sizes (1,21,41,61) or (1,6,10,15), with counts as close to equal as the integer budget allows. All allocations have rank-three information. D-optimal support is (1,10,61) or (1,7,15). For this simulation, nearby fully spent integer allocations on those supports are enumerated and the highest determinant is selected. This is distinct from the manuscript's local-search rounding comparison and does not certify global integer optimality.

The same joint constrained MLE is used throughout: prevalence in [1e-7,1-1e-7], sensitivity and specificity in [0.5,1]. These bounds approximate the closure of the manuscript model. Four common starting prevalences and data-derived accuracy starts are used. Feasible saturated three-size solutions are global binomial maxima; other fits use analytic-gradient optimisation, polishing and expanded-start checks. Estimation never receives the generating parameter. Nominal parameters equal generating parameters; design misspecification is not varied.

Seed: **20261110**, with independent L'Ecuyer-CMRG streams for the 12 cells. Boundary estimates are retained. Bias is the mean error; variance uses denominator R-1; RMSE uses mean squared error with denominator R. Monte Carlo SEs are reported for all three. The finite-R identity is `RMSE^2 = bias^2 + (R-1)/R * variance`.

There were 832 nonzero optimiser termination codes among 120,000 retained estimates. Expanded-start checks found likelihood differences below 2e-9. All projected gradients divided by assay count were below 5.3e-6. An independent prevalence-profile audit checked the first 50 datasets per cell and every large unnormalised-gradient case, finding no improvement above 1.8e-9. These numerical checks are not a symbolic global certificate for every dataset.

## Findings

At cost budget C=500, prevalence bias is 0.00307 under D versus 0.00838 under the uniform comparator; RMSE is 0.01930 versus 0.03815. Sensitivity RMSE is **smaller for the uniform design**: 0.05624 versus 0.07381. D-optimality targets generalized variance, not every coordinate's bias or RMSE.

At prevalence 0.022 and 500 assays, sensitivity bias is about -0.143 under both designs, with an accuracy boundary reached in about 88% of samples. The results do not establish that non-optimal design necessarily causes bias or that optimal design ensures unbiasedness.

With only one or two pool sizes, all three unknown parameters are not identified without external information. The technical note gives a likelihood ridge and proves that an estimator cannot be unbiased at two distinct parameter points with the same data distribution. This is separate from finite-sample bias in an identifiable design.

## Original code and revision corrections

Every file in the original `inst/` tree is preserved unchanged. Two distinctions matter:

1. The trailing `p_regret_grouped` plot in the original Stage II script uses an undefined `regret_group_tbl`. The runner reproduces the complete scenario calculations and preceding plots, records this limitation, and stops before that incomplete plot. `R/revision_checks.R` regenerates the revised scenario/regret figures from the original scenario table.
2. Original `table6.R` rescales an assay-optimised allocation to a specimen cap. The revised checks assess the specimen-specific optimisation problem. Original output is retained as provenance, not represented as the corrected result.

No participant-level data are required. Inputs are published planning parameters and simulated binomial outcomes. Code is distributed under GPL-3; see `LICENSE` and the retained package authorship metadata.
