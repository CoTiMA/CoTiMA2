# Proposed ctsem API for structured outer-level random-effect covariance

## Motivation

The native Julia/Laplace hierarchy currently estimates a full covariance among
all parameters marked `indvarying_<level>`.  In CoTiMA2's standardized
stationary validation, three study effects varied simultaneously: two bounded
DRIFT coordinates and one raw coordinate defining the off-diagonal T0
correlation.  A deliberately small eight-study diagnostic was rank-limited
(eight outer-level units for nine population parameters) and showed flat
correlation directions.  In the 16-study full validation, however, every
direction was identified and the full covariance fit converged with an accurate
Laplace check.  Structured covariance is therefore requested for parsimony,
sensitivity analysis, and smaller-study designs, not as a repair for a failure
of the native backend.

CoTiMA2 already offers independent, equicorrelated, full, and block covariance
structures in its separate Stan random-effects workflow.  These controls
cannot currently be translated to ctsem's native multilevel Julia model.

## Requested capability

Please expose a model- or `ctFit()`-level specification for the covariance of
random effects at each grouping level.  At minimum, the native Julia/Laplace
route needs:

1. a full covariance (the current behaviour);
2. an independent covariance, with every off-diagonal correlation fixed to
   zero; and
3. a block structure in which within-block correlations can be full or fixed
   to zero and between-block correlations are zero.

An equicorrelation block would be useful but is less urgent.

A possible interface would be conceptually similar to:

```r
model$randomEffectCovariance <- list(
  study = list(
    DRIFT = "full",
    T0VAR = "independent",
    between_blocks = 0
  ),
  id = "full"
)
```

The exact spelling is unimportant.  What matters is that the preparation,
parameter counting, raw-parameter naming, uncertainty calculation,
`ctLaplaceCheck()`, `ctIdentify()`, and `summary()` all use the same reduced set
of free correlation coordinates.

## Backward compatibility

Omitting the new specification should preserve the current full covariance.
Fixed-zero correlations should not appear as estimated parameters or be
reported as unidentified directions.  The fitted object should retain the
resolved structure for reproducibility.

## Why `poprank` is not a substitute

`poprank` imposes a reduced-rank regression representation.  It does not
represent a full-rank diagonal covariance with independently varying effects.
Using it here would change the marginal random-effects model rather than merely
remove unsupported or unnecessary correlations.

## Validation target

With the three study effects above and an independent study covariance, the
target validation is:

- all three marginal study SDs recovered without free correlation parameters;
- no structural or local unidentified correlation directions;
- the same stationary T0 covariance and derived DIFFUSION expressions;
- convergence with a small gradient norm; and
- a small Laplace-versus-quadrature gap per study/subject.
