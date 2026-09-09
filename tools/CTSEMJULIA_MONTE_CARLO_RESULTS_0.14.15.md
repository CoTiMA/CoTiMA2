# CoTiMA2 0.14.15 native ctsem/Julia Monte Carlo validation

## Purpose and design

This report records the first replicated parameter-recovery experiment for
CoTiMA2's experimental `backend = "ctsemjulia"`. It is an engineering
validation, not a general validation of hierarchical continuous-time models.

- CoTiMA2 0.14.15, repository commit `d3c4eb1`
- ctsem development `juliaFit` branch; Julia 1.12.7
- Ubuntu, AMD Ryzen Threadripper 3990X, 64 physical cores, 251 GiB RAM
- 8 concurrent workers and 4 Julia cores requested per fit
- 50 independently simulated data sets per scenario
- 16 studies, 35 persons per study, and 20 observations per person
- independent normal study effects on ctsem's unconstrained raw scale
- separate CINT, DRIFT, DIFFUSION, and T0VAR scenarios
- person-level CINT SD = 0.4 in the CINT and DRIFT scenarios

For nonlinear parameters, `summary(fit)$popsd.study` reports the induced SD
after transformation. The recovery target is therefore the integrated
transformed-scale population SD, not the raw Gaussian SD.

## Results

All 200 fits completed and converged.

| Scenario | True reported SD | Mean estimate | Bias | Relative bias | RMSE | 95% coverage | Wilson interval |
|---|---:|---:|---:|---:|---:|---:|---:|
| CINT | 0.3000 | 0.2955 | -0.0045 | -1.51% | 0.0477 | 0.98 | [0.895, 0.996] |
| DIFFUSION | 0.1620 | 0.1559 | -0.0061 | -3.79% | 0.0484 | 0.90 | [0.786, 0.957] |
| DRIFT | 0.4260 | 0.4139 | -0.0120 | -2.83% | 0.0957 | 0.92 | [0.812, 0.968] |
| T0VAR | 0.4654 | 0.4592 | -0.0063 | -1.35% | 0.1124 | 0.94 | [0.838, 0.979] |

The Monte Carlo standard errors of the four biases were 0.0068, 0.0069,
0.0136, and 0.0160. Every corresponding approximate 95% Monte Carlo interval
included zero. With only 50 replications, coverage remains imprecise; the
Wilson intervals should therefore accompany the point estimates.

Person-level CINT bias, RMSE, and coverage were 0.0052, 0.0178, and 0.98 in
the CINT scenario, and -0.0030, 0.0147, and 0.94 in the DRIFT scenario.

The mean absolute Laplace-check gap per person was `8.4e-16` (CINT),
`1.5e-5` (DIFFUSION), `1.4e-3` (DRIFT), and `1.3e-4` (T0VAR), all below the
suite's engineering warning threshold of 0.01 per person.

## Interpretation and limits

The results support the native Julia/Laplace path for the tested univariate,
single-family random-effect scenarios. They do not yet validate simultaneous
study variation in several parameter families, recovery of their
correlations, multivariate systems, boundary cases, or difficult nonlinear
models.

The next stage is `tools/ctsemjulia-joint-validation.R`. It simulates
independent study effects for all four families simultaneously. The current
ctsem API estimates a full raw-scale correlation matrix, so the test evaluates
both SD recovery and whether the six estimated correlations remain near zero.

## Provenance

The complete result object is stored outside the repository at
`/home/cdormann/SynologyDrive/R/CoTiMA2-Benchmark/ubuntu-mc-random50/monte-carlo-results.rds`.
The run started on 2026-09-09 with base seed 20260909 and
`--study-effects=random`.
