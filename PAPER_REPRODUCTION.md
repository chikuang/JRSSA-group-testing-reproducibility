# Tables and figures in the revised paper

Numbering is taken from the compiled `AI_assisted/main_AI_assist.tex` and
`AI_assisted/suppl_AI_assist.tex` synchronized with Overleaf on **11 November 2026**.
LaTeX labels are included because numbers can change in later revisions.

Run `Rscript reproduce.R list` to list the 33 items, or e.g.
`Rscript reproduce.R main-table-06` to execute one item's command. Install the
archived package first with `Rscript scripts/install_dependencies.R`.

**Coverage matters:** “Computed” means the indicated numerical calculation is run
at the stated inputs, not that every historical support weight, integer tie or
machine timing is bit-for-bit identical. “Verification” treats printed allocations
as inputs and checks them; it does not reconstruct them by optimization.
“From saved simulation” rebuilds the display from the supplied 10000-replicate
results. `Rscript run.R bias 10000 4` regenerates those simulations. Schematics and
literature summaries are authored in LaTeX and do not have R calculations.

The original `inst/` files remain unchanged. Their old filename numbers and some
header comments do **not** describe the current manuscript. In particular,
`table1.R` currently sets M=150,q=0.8 (one block of current Table 6), and `table2.R`
sets M=61,q=0.8 (one block of current Table 7). The new `R/paper_designs.R` loops over
all relevant blocks. Original `table5.R` omits E-optimality and the coordinate
columns; `R/low_prevalence_tables.R` supplies these for current Table 2. Original
`table6.R` is a legacy assay-optimal-to-specimen rescaling and is **not** the
current specimen-constrained Table 3 calculation. Original `table7.R` and
`table8.R` both calculate assay-constrained allocations related to current Table 4.

## Main text

| Item | LaTeX label | R source | Output | Coverage |
|---|---|---|---|---|
| Main Table 1 | `tab:identification-comparison` | `R/validate_theory.R` | `results/bias/three_size_comparison.csv` | Computed |
| Main Table 2 | `tbl:app-approximate` | `R/low_prevalence_tables.R` | `results/paper/main_table_02.csv` | Computed |
| Main Table 3 | `tbl:app-exact-individual` | `R/low_prevalence_tables.R` | `results/paper/main_table_03_verification.csv` | Verification |
| Main Table 4 | `tbl:app-num_tests` | `R/low_prevalence_tables.R`; `inst/paper_reproduce/table7.R` | `results/paper/main_table_04_verification.csv; results/original/table7/table7.csv` | Verification |
| Main Table 5 | `tab:analytic-bounds` | `R/revision_checks.R` | `results/revision/analytic_checks.txt` | Computed |
| Main Table 6 | `tbl:1-2` | `R/paper_designs.R` | `results/paper/main-table-06/` | Computed |
| Main Table 7 | `tbl:multiobj_maximin_time` | `R/paper_designs.R` | `results/paper/main-table-07/` | Computed |
| Main Table 8 | `tab:rounding-comparison` | `R/revision_checks.R` | `results/revision/rounding_R.csv` | Computed |
| Main Table 9 | `tab:stage2_grouped` | `R/stage2_outputs.R` | `results/paper/main_table_09.csv` | Computed |
| Main Table 10 | `tab:finite-bias-main` | `R/bias_simulation.R`; `R/summarize_bias.R` | `docs/main_bias_rows.tex; results/bias/summary.csv` | From saved simulation |
| Main Figure 1 | `tikz:sequential` | `Overleaf: AI_assisted/main_AI_assist.tex` | `Overleaf Figure 1` | LaTeX schematic |
| Main Figure 2 | `fig:equivalence3` | `inst/paper_reproduce/figure2.R` | `results/original/figure2/plots.pdf` | Computed |
| Main Figure 3 | `fig:eff-MO-exact` | `inst/paper_reproduce/figure3.R` | `results/original/figure3/g_n.pdf; results/original/figure3/g_C.pdf` | Computed |
| Main Figure 4 | `fig:stage2_regret` | `R/stage2_outputs.R` | `figures/main_figure_04.pdf` | Computed |
## Supplementary text

| Item | LaTeX label | R source | Output | Coverage |
|---|---|---|---|---|
| Supplementary Table S1 | `tbl:MOD_rounding_combined` | `inst/paper_reproduce/table4.R` | `results/original/table4/` | Computed |
| Supplementary Table S2 | `tbl:MOD_rounding_budget-2` | `inst/paper_reproduce/table3.R` | `results/original/table3/` | Computed |
| Supplementary Table S3 | `tbl:multiobj_maximin_time-supp` | `R/paper_designs.R` | `results/paper/supp-table-S03/` | Computed |
| Supplementary Table S4 | `tbl:SOD_rounding_costfree-supp` | `R/supp_rounding.R` | `results/paper/supp-table-S04/allocations.csv` | Recomputed designs |
| Supplementary Table S5 | `tbl:MOD_rounding_combined-supp` | `inst/paper_reproduce/tableS2.R` | `results/original/tableS2/` | Computed |
| Supplementary Table S6 | `tbl:MOD_rounding_budget-2-supp` | `R/supp_rounding.R` | `results/paper/supp-table-S06/allocations.csv` | Recomputed designs |
| Supplementary Table S7 | `tab:stage2_scenarios_full` | `R/stage2_outputs.R` | `results/paper/supp_table_S07.csv` | Computed |
| Supplementary Table S8 | `tab:literature-boundary` | `Overleaf: AI_assisted/suppl_AI_assist.tex` | `Overleaf Table S8` | Literature summary |
| Supplementary Table S9 | `tab:all-rounding` | `R/revision_checks.R` | `results/revision/rounding_R.csv` | Computed |
| Supplementary Table S10 | `tab:bias-cost-full` | `R/bias_simulation.R`; `R/profile_audit.R`; `R/summarize_bias.R` | `docs/bias_cost_rows.tex` | From saved simulation |
| Supplementary Table S11 | `tab:bias-low-full` | `R/bias_simulation.R`; `R/profile_audit.R`; `R/summarize_bias.R` | `docs/bias_low_rows.tex` | From saved simulation |
| Supplementary Table S12 | `tab:bias-allocations` | `R/bias_simulation.R`; `R/profile_audit.R`; `R/summarize_bias.R` | `docs/allocation_rows.tex` | From saved simulation |
| Supplementary Table S13 | `tab:bias-diagnostics` | `R/bias_simulation.R`; `R/profile_audit.R`; `R/summarize_bias.R` | `docs/diagnostic_rows.tex` | From saved simulation |
| Supplementary Figure S1 | `fig:equivalence2-supp` | `R/equivalence_plot.R` | `figures/supp_figure_S01.pdf` | Computed |
| Supplementary Figure S2 | `fig:stage2_k_full` | `R/revision_checks.R` | `figures/stage2_scenario_rebuilt.pdf` | Computed |
| Supplementary Figure S3 | `fig:stage2_regret_full` | `R/revision_checks.R` | `figures/stage2_regret_rebuilt.pdf` | Computed |
| Supplementary Figure S4 | `fig:stage2_actual_full` | `R/revision_checks.R` | `figures/stage2_expected_rebuilt.pdf` | Computed |
| Supplementary Figure S5 | `fig:bias-cost` | `R/summarize_bias.R` | `figures/bias_cost.pdf` | From saved simulation |
| Supplementary Figure S6 | `fig:bias-low` | `R/summarize_bias.R` | `figures/bias_low.pdf` | From saved simulation |

## Commands and scope

### Main Table 1 — `main-table-01`

All five allocations, including the rank-one and rank-two cases.

```sh
Rscript R/validate_theory.R vendor/gtDesign results/bias 10000 4
```

### Main Table 2 — `main-table-02`

All five criteria and coordinate efficiencies. Numerical supports may differ near solver tolerance; full precision is retained. A and E share rounded weights in the manuscript but have distinct specificity efficiencies (about 0.36554 and 0.36547); the output keeps them separate. The c calculation uses a conditioned basis and the correspondingly transformed contrast.

```sh
Rscript R/low_prevalence_tables.R
```

### Main Table 3 — `main-table-03`

Checks all displayed specimen allocations and information. Counts are manuscript inputs, not reconstructed optima. Original table6.R solves an older rescaling problem.

```sh
Rscript R/low_prevalence_tables.R
```

### Main Table 4 — `main-table-04`

Checks displayed counts. table7.R recomputes alternatives; its D tie uses 3840 specimens and its E allocation can include size 6, unlike the printed allocation.

```sh
Rscript R/low_prevalence_tables.R
```

### Main Table 5 — `main-table-05`

Base-R evaluation of the proved bound, no optimization needed.

```sh
Rscript run.R revision
```

### Main Table 6 — `main-table-06`

All M=61,150 and q=0,0.2,0.8 blocks. Solver-dependent near ties and computation times can differ from printed values.

```sh
Rscript R/paper_designs.R main-table-06
```

### Main Table 7 — `main-table-07`

All M=61,150 and q=0,0.2,0.8 blocks. Solver-dependent near ties and computation times can differ from printed values.

```sh
Rscript R/paper_designs.R main-table-07
```

### Main Table 8 — `main-table-08`

The main table selects D and E from the full eight-case comparison.

```sh
Rscript run.R revision
```

### Main Table 9 — `main-table-09`

Uses gtDesign response probabilities and the same ceiling(N/x)+N*pi(x) approximation as the manuscript.

```sh
Rscript R/stage2_outputs.R
```

### Main Table 10 — `main-table-10`

Selects budget 500. To regenerate all simulated data first: Rscript run.R bias 10000 4.

```sh
Rscript R/summarize_bias.R
```

### Main Figure 1 — `main-figure-01`

Authored TikZ diagram in the manuscript; no simulated or computed R result.

### Main Figure 2 — `main-figure-02`

D-A-Ds at M=150,q=0.2. The original code calls the prevalence contrast c with c=(1,0,0), which is Ds here.

```sh
Rscript scripts/run_original.R inst/paper_reproduce/figure2.R
```

### Main Figure 3 — `main-figure-03`

Both panels; n=10,...,100 and C=100,...,1000. Dashed line is the recomputed approximate minimum efficiency.

```sh
Rscript scripts/run_original.R inst/paper_reproduce/figure3.R
```

### Main Figure 4 — `main-figure-04`

Recreates the grouped plot without the undefined regret_group_tbl in the preserved original script.

```sh
Rscript R/stage2_outputs.R
```

### Supplementary Table S1 — `supp-table-S01`

Uses the original parameter settings. Integer ties and numerical tolerances can select alternative counts.

```sh
Rscript scripts/run_original.R inst/paper_reproduce/table4.R
```

### Supplementary Table S2 — `supp-table-S02`

Uses the original parameter settings. Integer ties and numerical tolerances can select alternative counts.

```sh
Rscript scripts/run_original.R inst/paper_reproduce/table3.R
```

### Supplementary Table S3 — `supp-table-S03`

All M=61,150 and q=0,0.2,0.6,0.8; D-Ds and D-Ds-c criteria.

```sh
Rscript R/paper_designs.R supp-table-S03
```

### Supplementary Table S4 — `supp-table-S04`

Calls package rounding for all printed budgets. Near-tied supports and counts may differ from the historical printed allocations; no global integer-optimality claim.

```sh
Rscript R/supp_rounding.R supp-table-S04
```

### Supplementary Table S5 — `supp-table-S05`

Uses the original parameter settings. Integer ties and numerical tolerances can select alternative counts.

```sh
Rscript scripts/run_original.R inst/paper_reproduce/tableS2.R
```

### Supplementary Table S6 — `supp-table-S06`

Calls package rounding for all printed budgets. Near-tied supports and counts may differ from the historical printed allocations; no global integer-optimality claim.

```sh
Rscript R/supp_rounding.R supp-table-S06
```

### Supplementary Table S7 — `supp-table-S07`

All original parameter scenarios, with the package response probability.

```sh
Rscript R/stage2_outputs.R
```

### Supplementary Table S8 — `supp-table-S08`

Authored comparison of cited studies; no R numerical calculation.

### Supplementary Table S9 — `supp-table-S09`

All exact counts, expenditures and eight comparisons.

```sh
Rscript run.R revision
```

### Supplementary Table S10 — `supp-table-S10`

Summarizes the saved 10000-replicate cells. Full rerun: Rscript run.R bias 10000 4.

```sh
Rscript R/summarize_bias.R
```

### Supplementary Table S11 — `supp-table-S11`

Summarizes the saved 10000-replicate cells. Full rerun: Rscript run.R bias 10000 4.

```sh
Rscript R/summarize_bias.R
```

### Supplementary Table S12 — `supp-table-S12`

Summarizes the saved 10000-replicate cells. Full rerun: Rscript run.R bias 10000 4.

```sh
Rscript R/summarize_bias.R
```

### Supplementary Table S13 — `supp-table-S13`

Summarizes the saved 10000-replicate cells. Full rerun: Rscript run.R bias 10000 4.

```sh
Rscript R/summarize_bias.R
```

### Supplementary Figure S1 — `supp-figure-S01`

D-A at M=61,q=0.8, using maximin_design_workflow and its equivalence plot.

```sh
Rscript R/equivalence_plot.R
```

### Supplementary Figure S2 — `supp-figure-S02`

Scenario plot corresponding to Supplementary Table S7.

```sh
Rscript run.R revision
```

### Supplementary Figure S3 — `supp-figure-S03`

Scenario plot corresponding to Supplementary Table S7.

```sh
Rscript run.R revision
```

### Supplementary Figure S4 — `supp-figure-S04`

Scenario plot corresponding to Supplementary Table S7.

```sh
Rscript run.R revision
```

### Supplementary Figure S5 — `supp-figure-S05`

All three parameters and three error measures across budgets; bars show Monte Carlo uncertainty.

```sh
Rscript R/summarize_bias.R
```

### Supplementary Figure S6 — `supp-figure-S06`

All three parameters and three error measures across budgets; bars show Monte Carlo uncertainty.

```sh
Rscript R/summarize_bias.R
```
