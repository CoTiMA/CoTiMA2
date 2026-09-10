# CoTiMA2 0.14.15 joint native ctsem/Julia validation

## Design

This experiment simultaneously estimated study-level random effects in CINT,
DRIFT, DIFFUSION, and T0VAR, a freely estimated 4 by 4 raw-scale correlation
matrix, and a person-level CINT random effect. Deterministic orthogonal study
effects made all six true study correlations exactly zero.

The quick design used 8 studies, 20 persons per study, and 12 observations per
person. The full design used 16 studies, 35 persons per study, and 20
observations per person. Both used five nodes in `ctLaplaceCheck()`.

## Recovery

| Design | Component | True reported SD | Estimate | Relative error | Covers truth |
|---|---|---:|---:|---:|---:|
| Quick | CINT | 0.300 | 0.253 | 15.7% | yes |
| Quick | DRIFT | 0.426 | 0.408 | 4.2% | yes |
| Quick | DIFFUSION | 0.162 | 0.164 | 1.2% | yes |
| Quick | T0VAR | 0.465 | 0.461 | 1.0% | yes |
| Full | CINT | 0.300 | 0.297 | 1.0% | yes |
| Full | DRIFT | 0.426 | 0.405 | 4.9% | yes |
| Full | DIFFUSION | 0.162 | 0.159 | 1.9% | yes |
| Full | T0VAR | 0.465 | 0.416 | 10.6% | yes |

The person-level CINT SD was 0.398 in the quick design and 0.371 in the full
design, compared with its true value of 0.400.

With only eight studies, two raw correlation estimates were sizeable (-0.304
and 0.499). With 16 studies, all six estimates were between -0.091 and 0.003.
This supports the expected interpretation that covariance recovery is weak in
the very small quick design but becomes substantially more stable with more
studies.

## Numerical diagnostics

| Design | Converged | Gradient norm | Fit time | Check time | Gap/person |
|---|---:|---:|---:|---:|---:|
| Quick | yes | 5.58e-9 | 91.7 s | 2,545.7 s | 0.03980 |
| Full | yes | 5.32e-9 | 281.6 s | 14,143.7 s | 0.03155 |

The optimizations were numerically clean and parameter recovery was good. The
Laplace-check gap nevertheless remained above the suite's engineering warning
threshold of 0.01 per person. Increasing the number of studies reduced but did
not remove it. The next diagnostic therefore compares nested random-effect
sets with `tools/ctsemjulia-dimensional-validation.R` to determine whether the
gap grows primarily with integration dimension or enters with a particular
parameter combination.

## Stored results

- Quick: `/home/cdormann/SynologyDrive/R/CoTiMA2-Benchmark/ubuntu-joint-quick/`
- Full: `/home/cdormann/SynologyDrive/R/CoTiMA2-Benchmark/ubuntu-joint-full/`

The full result retained its fitted ctsem object for later diagnostic work.
