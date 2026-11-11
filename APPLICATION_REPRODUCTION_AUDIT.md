# Application-table reproduction audit

Checked against the preserved `JRSSA_rev2.tex` and the current
`AI_assisted/main_AI_assist.tex` downloaded from Overleaf on 11 November 2026.
The table labels are `tbl:app-exact-individual` and `tbl:app-num_tests`.
The original reproduction scripts call them Tables 6–7; in the preserved
`JRSSA_rev2.tex`, intervening additions place them at Tables 7–8. They are
Tables 3–4 in the current AI_assisted manuscript. The stable labels identify
the comparison independently of these numbering changes.

## What the original code reproduces

The author's original `inst/paper_reproduce/table6.R` and `table7.R` remain
unchanged. With the installed gtDesign 0.1.1, R 4.6.1, CVXR 1.9.3 and CLARABEL,
seven of the ten criterion rows reproduce the original counts directly.
The remaining three are explained by solver accuracy and a D-optimal tie.
With the E solver settings and D tie rule below, all ten rows match the original
manuscript exactly, including assay and specimen totals.

- **Original `table6.R` (specimen-capped rescaling):** D, A, c and Ds match at defaults. E also matches after
  setting `ridge=0` and `tol_gap_abs=tol_gap_rel=tol_feas=1e-12`.
- **Original `table7.R` (500 assays):** A, c and Ds match at defaults. The same E settings recover
  its printed E allocation. D returns an equally efficient permutation; choosing
  the one using the fewest specimens recovers the printed allocation.

The default E solve leaves weight 0.0000625115 at pool size 6. With
`fix_zero_floor=TRUE`, rounding transfers an assay to this numerical support.
The higher-accuracy solution has support 1, 7, 15, and this extra assay disappears.
`E_solver_weights.csv` and `E_rounding_trace.csv` show the two paths. This is a
support/rounding sensitivity, not a difference in the probability model.

For D with 500 assays at sizes 1, 7, 15, the three permutations of counts
(167,167,166) have exactly the same determinant: on three nonsingular support
points, `det(J) = det(H)^2 * n1*n2*n3`, where the rows of H are the regressors.
Their specimen totals are 3826, 3834 and 3840. The original manuscript prints
3826. The package's enumeration order can select 3840 without changing D-efficiency.

These findings qualify the earlier repository wording that Tables 3–4 were only
verified. Differences from the revised specimen table did not show that the
original code failed to reproduce the original manuscript.

## Why the specimen table changed in the revision

Original `table6.R` computes approximate designs with q=0 (information per assay),
then `round_gt_design_subject_budget` converts those weights to an assay count
and selects a completion under the 6632-specimen cap. This reproduces the old table.

The revised Main Table 3 instead computes the approximate design with q=1,
so the cost is pool size and the objective uses information per specimen.
It then calls `round_gt_design_budget` with C=6632 and q=1. The two procedures
solve different optimization problems. For example, the revised D allocation
uses 2842 assays at sizes 1, 4, 5, 15, whereas the original uses 865 assays at
sizes 1, 7, 8, 15. Both use 6632 specimens. Neither comparison imposes a separate
500-assay cap.

## Complete construction of current Main Tables 3–4

`R/application_exact_tables.R` starts from the nominal parameter vector,
candidate sizes and resource limits. It does not read printed allocations.
It calls the public gtDesign optimizers and budget-rounding function for all
five criteria in each table. Every generated count matches the current manuscript.

Numerical settings are explicit:

- D uses a nonsingular congruence transformation, zero ridge, support threshold
  1e-6 and solver tolerances 1e-11. Congruence preserves the D optimizer.
- A uses an orthogonal rotation from the regressor SVD, support threshold 1e-6
  and tolerances 1e-9. Orthogonality preserves the trace criterion on the original
  parameter scale; arbitrary coordinate rescaling would not.
- E uses the original regressors, zero ridge, threshold 1e-6 and tolerances 1e-12.
- c and Ds use a conditioned basis with the corresponding transformed contrast,
  threshold 1e-6 and tolerances 1e-10.
- Rounding uses `n_index=2` and `fix_zero_floor=FALSE`. It keeps the floor counts
  fixed and searches for feasible local completions. This is not a global
  integer-optimality certificate.
- For D at 500 assays, candidate losses within 1e-9 in negative log determinant
  are ties; the candidate with the fewest specimens is selected. All three tied
  candidates and the selected one are saved. Here their losses differ by less
  than 2e-13 numerically and are equal by the determinant formula above.

All ten approximate solves have status `optimal`. Directional scores are checked
on all 15 candidate sizes using the original parameter scale; the largest relative
excess is 5.75e-6 or less. Every implemented design has rank-three information,
integer counts and exact expenditure of its stated budget. All ten efficiencies
against their own approximate benchmarks exceed 0.99999.

## Run and inspect

From the repository root, after dependency installation:

```sh
# Recompute both current tables, including approximate designs and rounding.
Rscript reproduce.R main-table-03
Rscript reproduce.R main-table-04

# Rerun the original scripts, then compare to the original manuscript.
Rscript scripts/run_original.R inst/paper_reproduce/table6.R
Rscript scripts/run_original.R inst/paper_reproduce/table7.R
Rscript R/audit_original_application.R
```

Current outputs in `results/paper/`:

- `main_table_03.csv`, `main_table_04.csv`: counts, resources, efficiencies, local
  standard errors and diagnostics; corresponding `_rows.tex` files are generated.
- `application_weights.csv`, `application_certificates.csv`: complete approximate
  weights and full-candidate checks.
- `application_designs.rds`: optimizer objects, floor allocations, candidate
  completions and selected exact designs.
- `main_table_04_D_ties.csv`: all tied D candidates and selection.
- `application_sessionInfo.txt`: numerical environment.

The `results/paper/historical_application/` directory contains the ten-row
original-manuscript comparison, default/tight E weights and rounding traces,
tight-E versions of the original tables, and session information. The reference
counts in the audit are used only for comparison after optimization.

All 13 preserved original R entry scripts were also executed. That execution
check is separate from the detailed manuscript comparison above, which covers
the two application tables. The original trailing Stage II plot still has its
previously documented undefined variable; the replacement plotting code is
identified in the Table/Figure index.
