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

## Quick result

The quick design used 8 studies, 20 subjects per study, 12 observations per
subject, four requested Julia cores, and deterministic study quantiles.

| Scenario | Converged | Gradient norm | Fit seconds | Gap/person | True raw study SD | Estimated study SD |
|---|---:|---:|---:|---:|---|---|
| Univariate DRIFT | yes | 1.19e-9 | 50.1 | 5.27e-4 | 0.300 | 0.318 |
| Bivariate DRIFT | yes | 8.39e-9 | 61.1 | -1.69e-5 | 0.200, 0.160 | 0.153, 0.080 |

Both optimizations and both adaptive-quadrature checks passed numerically.
The small bivariate design underestimates especially the second study SD, so
the full design is required before drawing a recovery conclusion.

## Remaining API limitation

The native Julia backend currently rejects a parameter-dependent expression in
T0VAR. Hence a study-varying T0 correlation cannot yet be combined with both
marginal latent variances fixed exactly at one through the required Cholesky
identity. The validation fixes the T0 correlation rather than fitting a free
Cholesky element that would change a marginal variance.

These assumptions apply to latent covariance stationarity. Standardizing
observed variables at each occasion does not by itself impose them in the
presence of measurement error, sampling variation, or nonstationary latent
dynamics.
