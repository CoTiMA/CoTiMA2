# Proposed API: covariance-stationary standardized latent models

Status: design proposal for discussion; no public behaviour is changed by this
document.

## Purpose and terminology

CoTiMA2 commonly analyses variables standardized within measurement occasions.
That data preprocessing does not imply that the latent continuous-time process
has unit marginal variances or is covariance-stationary.  The proposed option
therefore uses the explicit name `stationaryStandardized`, and never infers the
constraint merely from standardized input data.

For study *s*, the constrained model defines a stable drift matrix `A_s` and a
target stationary latent covariance `P_s`, then derives the complete process
noise covariance

```text
Q_s = -(A_s P_s + P_s A_s')
```

Consequently, `DIFFUSION` is derived rather than independently estimated.  At
time zero, the model imposes `T0VAR_s = P_s`.

## Proposed user-facing argument

Add the following top-level argument to `ctmaFit()`:

```r
stationaryStandardized = c("none", "latent")
```

The default is `"none"`, preserving current behaviour exactly.  `"latent"`
means that the user explicitly assumes a covariance-stationary latent process
with unit marginal latent variances.  A logical compatibility shorthand may be
accepted temporarily (`FALSE` -> `"none"`, `TRUE` -> `"latent"`), but the
returned call and documentation should use the character values.

Do not call this option simply `standardized`: that would incorrectly suggest
that standardizing observed variables is sufficient.  Do not use `"auto"` in
the first public version because there is no reliable way to infer the latent
assumption from the data.

The option is a model-construction argument, not a `backendArgs` control.  The
same requested statistical model should ultimately be available to every
backend capable of representing it.  Initially it should be enabled only for
`backend = "ctsemjulia"`; other backends should stop with an informative
capability error rather than silently fit another model.

## Target covariance

For the first implementation, use

```r
stationaryCorrelation = NULL
```

with these meanings:

- `NULL`: identity target covariance, `P_s = I`.  All latent marginal
  variances equal one and all stationary correlations equal zero.
- A fixed positive-definite correlation matrix: common `P` across studies.
- A model or formula for study-varying correlations: eligible for the next
  validation stage now that the native Julia backend accepts and correctly
  differentiates static parameter-dependent expressions in `T0VAR`.

The diagonal is always fixed to one.  Covariance input is not accepted because
the option specifically represents standardized latent variables.

## Interaction with study fixed effects

`studyCoding` continues to determine how ordinary study indicators are coded;
it is not replaced by `stationaryStandardized`.

`studyVarying` identifies independently estimated study fixed effects.  Under
`stationaryStandardized = "latent"`:

- `"DRIFT"` is allowed.  Study-specific `A_s` values imply study-specific
  derived `Q_s` values.
- `"CINT"` and mean-related families retain their existing meanings.
- `"DIFFUSION"` is rejected because diffusion is derived from `A_s` and `P_s`.
- `"T0VAR"` is rejected in the initial identity-`P` implementation because
  `T0VAR` is fixed by the stationary target.

Thus a typical fixed-effects call would be conceptually equivalent to:

```r
ctmaFit(
  ctmaInitFit = initFit,
  backend = "ctsemjulia",
  stationaryStandardized = "latent",
  studyCoding = "weighted",
  studyVarying = c("CINT", "DRIFT")
)
```

This parameter-selective restriction is preferable to silently deleting
families from `studyVarying`, because the latter could conceal a substantively
different model.

## Interaction with hierarchical random effects

For native Julia/Laplace study random effects, retain the existing matrix-mask
API and apply the following rules:

- `reDRIFT`: allowed.  Random study drift effects propagate deterministically
  into study-specific diffusion matrices.
- `reDIFFUSION`: every entry must be `FALSE` (or the argument must be `NULL`).
  A `TRUE` entry is an error because there is no independent diffusion random
  effect under the constraint.
- `reT0VAR`: every entry must be `FALSE` (or `NULL`) for the initial identity
  target.  When a future study-varying stationary-correlation model is
  supported, only its off-diagonal correlation coordinates may vary; diagonal
  variances remain fixed at one.

Random effects on DRIFT may be correlated using the existing `covariance`,
`covarianceBlocks`, `reBlocks`, and `blockStructures` arguments.  The induced
dependence between DRIFT and derived DIFFUSION must not be represented as a
second random-effect block: it is deterministic, not an independently
estimated covariance.

Person-level random effects remain conceptually possible.  The first release
should allow them only if the native backend's parameter transformations and
Laplace integration have passed the same recovery checks as the study-only
case.  A person-level random effect in DRIFT likewise induces a person-specific
`Q`; it must not be combined with an independent person-level DIFFUSION effect.

## Required validation before fitting

The model builder should fail early unless all of the following hold:

1. `A` is parameterized so that stability is guaranteed over the admissible
   parameter domain, not only at the initial values.
2. `P` is symmetric positive definite with unit diagonal.
3. The derived `Q` is symmetric positive semidefinite throughout the
   admissible parameter domain.  Checking only the starting point is
   insufficient.
4. `T0VAR` represents exactly `P`, including the transformation used by ctsem.
5. No independently free DIFFUSION parameter remains.
6. Requested fixed and random study effects obey the compatibility rules
   above.
7. The selected backend supports parameter-dependent expressions in all
   affected matrices and returns correct gradients.  This is now established
   for static `T0VAR` expressions on ctsem `juliaFit` commit `2140ceb`, but not
   for expressions depending on observation rows, time-dependent predictors,
   or latent states.

For multivariate models, stable `A` does not by itself guarantee positive
semidefinite `Q` for a fixed `P`.  The initial implementation should therefore
support only parameterizations with a mathematical admissibility guarantee.
The validated bivariate identity-`P` prototype achieved this by bounding the
diagonal drift terms away from zero and fixing sufficiently small cross-effects.
A general free multivariate DRIFT requires a more principled constrained
parameterization before it becomes a public option.

## Returned-object metadata and summary

Every fit should retain:

```r
fit$model_constraints$stationaryStandardized
fit$model_constraints$stationaryCorrelation
fit$model_constraints$diffusion_derived
fit$model_constraints$t0var_fixed_to_stationary_target
```

`summary()` should print a short note such as:

> The latent process was constrained to be covariance-stationary with unit
> marginal variances. DIFFUSION was derived from DRIFT and the stationary
> covariance and has no independent fixed or random effects.

Study-specific summaries should report DRIFT and the corresponding derived
DIFFUSION, clearly labelling the latter as derived.  Standard errors or
intervals for derived diffusion quantities must propagate uncertainty from all
free coordinates rather than treating `Q` as fixed.

## Capability staging

Recommended implementation order:

1. Identity `P`, univariate DRIFT, native `ctsemjulia`, study-level random
   effects.
2. Identity `P`, the already validated restricted bivariate DRIFT design.
3. Fixed nonzero common `P`; the prerequisite parameter-dependent-`T0VAR`
   smoke test has passed and the stationary model itself remains to be tested.
4. Study-varying stationary correlations with fixed diagonal one.
5. Person-level random effects and more general multivariate DRIFT
   parameterizations.
6. Equivalent Stan and CoTiMA2-Julia representations, followed by controlled
   cross-backend comparison.

Each stage requires simulation recovery, multiple starts, exact-gradient
checks, a Laplace-versus-quadrature diagnostic where feasible, and verification
that derived `Q` remains admissible at every evaluated solution.

## Open decisions before implementation

1. Should `stationaryCorrelation` be exposed in the first public version, or
   kept internal until nonzero correlations are supported?
2. Should incompatible `studyVarying`/random-effect selections always be hard
   errors, as proposed, or should an explicit helper construct the compatible
   masks for users?
3. Which constrained multivariate DRIFT parameterization should replace the
   deliberately narrow validation bounds?
4. Should the feature first appear as an experimental helper that returns a
   modified ctsem model, before becoming a direct `ctmaFit()` argument?

The conservative recommendation is to begin with an experimental model-builder
helper and promote the argument to `ctmaFit()` only after the updated ctsem API
and the first two recovery stages have passed.

## Native Laplace smoke-test result

On 2026-09-11, `tools/ctsemjulia-t0var-laplace-smoke.R` was run against ctsem
3.12.0 from `juliaFit` commit `2140ceb`.  A bivariate model used two static
`PARS` coordinates in the composed expressions

```text
T0VAR[1,1] = log1p_exp(aa + bb)
T0VAR[2,1] = tanh(aa - bb)
T0VAR[2,2] = log1p_exp(aa)
```

and two subject-level random manifest means.  With
`backend = "julia"` and `intoverpop = "laplace"`, preparation succeeded, the
likelihood was finite, and optimization converged with gradient norm
`9.68e-08`.  The fit took 166 seconds after first-time Julia precompilation.
This closes the previously observed static-`T0VAR` capability gap for the
Laplace route.  It does not imply support for row-varying `tdpreds` or latent
state references in `T0VAR`; those are evaluated at the wrong conceptual time
for a single initial covariance and remain correctly rejected.

The subsequent covariance-stationary smoke test also passed technically.  It
used two bounded study-varying diagonal DRIFT coordinates and a third
study-varying raw coordinate defining `rho = 0.3 * tanh(rho_raw)`.  T0VAR had
unit diagonal and correlation `rho`; the complete DIFFUSION factor was derived
from the same three coordinates through the Lyapunov equation.  With eight
studies and 160 people, optimization converged with gradient norm `5.54e-09`
and a Laplace--quadrature gap of `1.73e-04` per person.  The deliberately small
design did not identify two correlations among the three study random effects,
so it establishes backend capability but not recovery.  A larger recovery run
or an explicitly independent study-effect covariance structure is required
before promoting this parameterization.
