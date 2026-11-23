# Internal Lohmann et al. replication adapter

CoTiMA2 0.14.59. This function is intentionally **not exported**.
Use `CoTiMA2:::ctmaLohmannFit()`. It does not modify ctmaInit/ctmaFit or their
settings, and is not a recommended production estimator.

## Purpose and provenance

The adapter reads the literal `OS_FE_char`, `OS_URE_char`, `OS_ERE_char` or
`OS_CRE_char` Stan program from the locally downloaded online supplement to
Lohmann, Luedtke and Robitzsch (2026), *A One-Stage Approach for Estimating a
Meta-Analytic Continuous-Time Random-Intercept Cross-Lagged Panel Model*.
It **never sources** their R compiler (which would load packages and compile all
four models). Original Stan code is not redistributed here. The returned object
records the compiler path, MD5 and complete selected code. The supplied simulation
and empirical compiler files were identical (MD5
`1032bdbb25255e14e7efba9a1f69497b`). Updated files are read literally, not silently
patched; inspect their code before treating them as the same simulation condition.

This is replication of the supplied executable specification, **not certification
that it implements every formula in the paper correctly**.

## Inputs and example

```r
# Matrices: complete correlations in X1,Y1,X2,Y2,... order.
# Studies may have different numbers of waves, but no missing correlations.
# Each study has one constant positive time interval between successive waves.
# V must use the same column-major lower-triangle ordering as R[lower.tri(R)].
prepared <- CoTiMA2:::ctmaLohmannFit(
  correlationMatrices = R_list,
  sampleSizes = N_vector,
  timeLag = lag_vector,
  samplingCovariances = V_list,
  supplementDirectory = "/path/to/downloaded/supplement",
  reStructure = "URE",
  fit = FALSE
)

# Compile only the chosen model and run Stan/NUTS.
# This can be computationally expensive; choose resources deliberately.
result <- CoTiMA2:::ctmaLohmannFit(
  stanData = prepared$stanData,
  supplementDirectory = "/path/to/downloaded/supplement",
  reStructure = "URE", chains = 4, iter = 2000, warmup = 1000,
  coresToUse = 4, seed = 20261106,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)
result$summary
rstan::check_hmc_diagnostics(result$fit)
```

Alternatively supply `primaryStudies` from ctmaPrep (or a complete fitted CoTiMA
object convertible with ctmaFitToPrep), plus `samplingCovariances`. Its `empcovs`
must already be correlation matrices in the stated order; covariance matrices
are deliberately not silently standardized. Its `deltas` must describe constant
spacing and match the matrix dimensions. No scaleTime is applied: choose one time
unit for the input intervals, which also determines the interpretation of priors.
Raw-data-only studies, binary-indicator likelihoods, moderators and general ctsem
constraints are not supported. `...` rejects unknown settings rather than
silently overriding ctmaFit arguments.

`stanData` also accepts the authors' prepared data directly, with dimensions
Nstudies, maxwaves, Nlatent, Nmanifest, max_cors, n_cors, Nmeasurements, TimeLag,
r_i and V_i. Only active observed blocks are validated; -99 padding is allowed.
For an exact empirical reproduction, use the authors' prepared V_i. Their
empirical script uses a pooled tssem1FEM correlation matrix to construct sampling
covariances. Computing V separately from each observed R is a different choice.
The adapter therefore requires explicit V rather than choosing an approximation.
Sample sizes are checked in the matrix interface but do not separately enter the
Stan likelihood once V is supplied. Pairwise sample sizes must be reflected in V.

## Arguments and results

- `reStructure`: case-insensitive FE, URE, ERE, CRE. FE has no study RE; URE has
  independent RE; ERE has one common RE correlation; CRE estimates their full
  correlation matrix. These are study-level effects, separate from person-level
  random intercepts represented by trait variances in all four programs.
- `supplementDirectory`: archive root or a specific modelcompiler.R. Multiple
  identical copies are accepted; differing copies require an explicit file.
- `fit=FALSE`: prepare/validate and return code/data without compiling or sampling.
- `chains=4`, `iter=2000`, `warmup=floor(iter/2)`, `coresToUse=1`,
  `seed=20261106`, `control=list()`, `init=0`: Stan NUTS controls. Iter includes
  warmup. Default zero initialization is reproducible but not an admissibility
  guarantee. RE variants retain the authors' centered/noncentered choices.
- `stanModel=NULL`: optionally reuse a compiled rstan model; its exact code must
  match the selected source. This avoids compilation for each MC replication.
- `activeDirectory=NULL`: optionally save the complete result as
  ctmaLohmannFit_DIMD_STRUCTURE.rds in an existing directory (overwrites that file).

Results contain fit (stanfit, or NULL), summary (raw rstan summary), stanData,
modelCode, source fingerprints, call, structure and caveats. This deliberately
uses a separate result class; it is not a ctsem fit and generic CoTiMA summary or
table methods must not be assumed to understand it. Population coordinates are
transformed/raw quantities, not automatically original-scale population means.
The parameter_order field records study-deviation ordering: A11raw, A21, A12,
A22raw, trait-X logit, trait-Y logit, trait correlation raw, within correlation raw.

## Important audit findings: preserved, not silently repaired

The split loops for the state loadings together produce the identity loading
map in wave-major order; they are not a variable-ordering error.

1. The later-wave residual covariance contains a standalone `-A12*A21`, whereas
   Phi - B Phi B' contains `-Phi12*A12*A21`. Consequently the supplied recurrence
   need not maintain the claimed stationary covariance. Again retained literally.
2. FE uses untransformed/unbounded Theta_Cors and Normal(0,100) priors; RE variants
   use tanh correlations and Normal(0,1) population-coordinate priors. Thus FE vs
   RE differences cannot be attributed solely to the presence of heterogeneity.
3. Negative diagonal drift transformations do not guarantee a stable full drift.
   Neither do the separate trait/state transforms guarantee a positive-semidefinite
   continuous-time diffusion. The likelihood checks V, not the physical validity
   of every model-implied process. Successful MCMC does not settle this issue.
4. URE is centered; CRE/ERE use standardized noncentered deviations. ERE maps a
   Normal(0,1) raw correlation into (-1/7,1); CRE uses LKJ(2). All RE SD priors are
   positive-half Cauchy(0,2.5). These choices are preserved.

Do not label simulations from this adapter as a validated implementation of the
intended mathematical model. A corrected mathematical version would be a separate,
explicit simulation condition after discussing these findings with the authors.
No new population fits or posterior-recovery claims accompany this first adapter.

## Univariate extension and runnable examples (0.14.59)

`n.latent=2` (default) selects literal original code. `n.latent=1` selects a
newly authored univariate extension; it needs no downloaded supplement. Its
study-specific correlations are tau + (1-tau)*exp(A*abs(t-s)), with A=-exp(raw),
tau=inv_logit(raw), total variance 1, and Q=-2*A*(1-tau). It has two study-level
random-effect coordinates (drift raw, trait logit). URE is centered; CRE uses
LKJ(2); ERE maps Normal(0,1) to (-1,1); RE SDs have half-Cauchy(0,2.5) priors.
FE retains the original broad Normal(0,100) coordinate-prior convention.
For two RE coordinates, ERE and CRE have the same possible correlation matrices
but different priors and therefore are not identical Bayesian models.
The univariate covariance is mathematically stationary; this does not transfer
a validation claim to the unchanged bivariate supplement code.

Two complete scripts are installed alongside this document:
`ctmaLohmannFit-univariate.R` and `ctmaLohmannFit-bivariate.R`.
Use `system.file("research", "ctmaLohmannFit-univariate.R", package="CoTiMA2")`
to find them. They simulate sample correlations with explicit Gaussian plug-in
sampling covariances. They run preparation only by default; set runSampling=TRUE
to run NUTS. Their five-study design illustrates syntax, not reliable RE precision.
The bivariate script needs a local supplement path and deliberately distinguishes
its mathematically stationary generating model from the preserved fitting-code
covariance discrepancy. These are specification examples, not recovery evidence.
