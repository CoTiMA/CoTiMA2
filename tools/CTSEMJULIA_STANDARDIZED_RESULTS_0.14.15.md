# Covariance-stationary standardized validation (CoTiMA2 0.14.15)

## Model tested

This complementary validation assumes an idealized stationary latent process
with unit marginal variances. Process noise is derived rather than estimated
independently:

`Q = -(A P + P A')`.

The univariate model varies one DRIFT parameter between studies. The bivariate
model varies both diagonal DRIFT parameters, fixes the cross-effects at 0.08
and 0.04, and fixes the T0 correlation at zero. Both models fix measurement
error at zero. The data were generated exactly from their stationary
transition covariance.

## Full result and start-value robustness

The full design used 16 studies, 35 subjects per study, 20 observations per
subject, four requested Julia cores, and deterministic, mutually uncorrelated
study-effect sequences. Six optimizer configurations (with and without the
10-iteration careful-fit pass and with initial jitter from 0.01 to 0.20) all
converged to the same bivariate solution and log likelihood. The fastest was
the plain fit with jitter 0.05 (31.5 seconds, excluding Hessian and quadrature).

| Seed | Converged | Gradient norm | Fit seconds | Gap/person | True raw SDs | Estimated SDs | Estimated correlation |
|---|---:|---:|---:|---:|---|---|---:|
| 20260912 | yes | 7.72e-9 | 94.4 | 2.78e-4 | 0.200, 0.160 | 0.159, 0.169 | -0.152 |
| 20261010 | yes | 1.41e-9 | 105.1 | 3.47e-4 | 0.200, 0.160 | 0.214, 0.185 | -0.158 |

The corresponding raw population means were 0.326 and 0.163 in the first run,
and 0.311 and 0.139 in the alternate-seed run, for generating values 0.35 and
0.15. Both 95-percent intervals covered both generating study SDs and zero
correlation. The adaptive-quadrature comparisons show small Laplace error.

Numerical generator checks gave an exact unit diagonal for every target
covariance and a maximum Lyapunov/Cholesky reconstruction error of 4.44e-16.
Overall empirical variances were 0.970 and 0.981 for the first seed, and 0.999
and 1.012 for the second.

An exploratory formulation with a nonzero fixed T0 correlation allowed the
optimizer to enter regions where the implied Q was not positive definite.
That was a validation-model parameterization error, not a Laplace failure.
The final zero-correlation formulation keeps Q positive definite over the
complete parameter domain used here.

## Extended correlation validation

The ctsem development branch now accepts the required static,
parameter-dependent T0VAR expression.  The extended validation therefore
varies a bounded T0 correlation between studies while fixing both marginal
latent variances exactly at one and deriving DIFFUSION from the same stationary
covariance.  Six independent 16-study fits all converged.  For generating raw
study SDs `(0.20, 0.16, 0.18)`, the mean estimates were
`(0.188, 0.189, 0.183)` and the median estimates were
`(0.182, 0.200, 0.187)`.

Three of the six fits flagged one weak correlation in the unstructured
study-random-effect covariance: `rawcor_drift2__drift1.study` in two seeds and
`rawcor_rho_raw__drift1.study` in one.  This does not invalidate the marginal
SD recovery, but it supports adding independent or blockwise covariance
structures for parsimony and sensitivity analysis.  Static initial-state
expressions are supported; row-varying time-dependent predictors do not define
a single T0 covariance and remain outside this scope.

These assumptions apply to latent covariance stationarity. Standardizing
observed variables at each occasion does not by itself impose them in the
presence of measurement error, sampling variation, or nonstationary latent
dynamics.
