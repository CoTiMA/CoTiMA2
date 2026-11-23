# CoTiMA2 for Beginners: A Plain-English Guide

Reviewed against CoTiMA2 0.14.51, 2 November 2026.
Code examples are templates, not calculations run when this guide is opened.
The guide assumes that CoTiMA2 is installed and that you already have a study
list created by `ctmaPrep()`. Read the full User Guide for preparing that list.

## 1. What does CoTiMA2 do?

CoTiMA2 combines longitudinal evidence from primary studies that may use
different measurement intervals. Continuous-time effects describe change per
unit of time; their discrete-time counterparts depend on the chosen interval.

Three choices are separate:

- **Backend:** which numerical implementation performs the calculation?
- **Estimation:** maximum likelihood (ML), maximum a posteriori (MAP), or
  posterior sampling?
- **Model:** which parameters are common across studies, which differ between
  studies, and which vary between people?

Changing the backend does not intentionally change the scientific question.
However, supported models, integration methods and measurement approximations
can differ. Large differences between engines require investigation, especially
for binary outcomes, weakly identified models and random effects.

### A short glossary

- **DRIFT:** continuous-time autoregressive and cross-effects. Matrix row is
  the outcome; column is the predictor. These are not ordinary lag-1 regression
  coefficients.
- **DIFFUSION:** the process-noise parameterization. Report the transformed
  covariance (`DIFFUSIONcov`) when you mean a noise covariance; the model's
  DIFFUSION entries need not themselves be covariance entries.
- **T0MEANS / T0VAR:** initial-state means and the parameterization of initial
  covariance. Stable person effects can be represented separately, so a printed
  T0 covariance need not describe all observed baseline variation.
- **TI predictor:** a predictor constant within a person, such as a study
  indicator or a moderator.
- **Raw parameters:** internal estimation coordinates. They are transformed
  to obtain the quantities you interpret. A raw point estimate is not a set of
  uncertainty draws.
- **ML:** seeks parameters that maximize the likelihood. A returned answer is
  not necessarily a successful maximum: inspect convergence and identification.
- **MAP:** optimization with priors; not posterior sampling.
- **Bayesian sampling:** draws from a posterior distribution defined by the
  likelihood and priors. Draws require convergence and precision checks.

## 2. Choose a backend explicitly

| CoTiMA2 argument | Implementation | Main practical distinction |
| --- | --- | --- |
| `backend = "stan"` | ctsem's Stan route | Default for `ctmaInit()` and `ctmaFit()` when backend is omitted. |
| `backend = "julia"` | CoTiMA2's bundled Julia engine, connected through JuliaCall | Its own controls, optimizer and NUTS implementation. |
| `backend = "ctsemjulia"` | ctsem's native Julia engine, connected through JuliaConnectoR | Requires a compatible experimental ctsem installation; uses native optimizer and sampler controls. |

CoTiMA2's `backend = "ctsemjulia"` corresponds to `backend = "julia"` when
calling **ctsem directly**. The same word does not identify the same route in
the two packages. Installing experimental ctsem does not change CoTiMA2's
default backend.

Neither Julia route is guaranteed to be faster for every model. Compilation
can dominate the first small fit. Model support and compatible experimental
ctsem APIs can change between revisions.

### Setup: use the guide for the selected route

The two Julia routes have different setup helpers:

```r
# Bundled CoTiMA2 Julia engine: requires an installed Julia executable.
# Inspect first; consent permits dependency setup/precompilation.
CoTiMA2::ctmaJuliaSetup(dry_run = TRUE)
# CoTiMA2::ctmaJuliaSetup(agree = TRUE)

# Native ctsem Julia engine, after selecting a compatible ctsem installation:
# ctsem::ctJuliaSetup(threads = 2, agree = TRUE)
```

Do not run both setup blocks just because both use Julia. See
[Julia installation](Julia_Install.Rnw) and the
[Julia usage guide](CoTiMA2_Julia_Usage.Rnw) for installation instructions,
platform requirements and isolated R libraries. In particular, select the
experimental library **before** loading ctsem or CoTiMA2. Different ctsem
installations cannot be loaded simultaneously in one R session. A restart
must also avoid restoring a workspace that reloads the wrong namespace.

## 3. A first ML workflow

This example has two continuous variables and no person-level random effects.
Replace the output path and supply your own `studyList`.

```r
library(CoTiMA2)
activeDirectory <- "PATH/TO/AN/EXISTING/OUTPUT/FOLDER/"

# Fit each primary study separately.
initFit <- CoTiMA2::ctmaInit(
  primaryStudies = studyList,
  activeDirectory = activeDirectory,
  n.latent = 2,
  backend = "stan",
  optimize = TRUE,
  priors = FALSE,
  coresToUse = 2
)
saveRDS(initFit, file.path(activeDirectory, "initFit.rds"))

# Fit a pooled model. For native Julia, first complete its setup above.
pooledFit <- CoTiMA2::ctmaFit(
  ctmaInitFit = initFit,
  activeDirectory = activeDirectory,
  backend = "ctsemjulia",
  optimize = TRUE,
  priors = FALSE,
  coresToUse = 2
)
saveRDS(pooledFit, file.path(activeDirectory, "pooledFit.rds"))

summary(pooledFit)
CoTiMA2::ctmaFitDiagnostics(pooledFit)
```

Both `ctmaInit()` and `ctmaFit()` accept all three backend choices. Use the
same substantive specification when comparing them: variables, binary flags,
time scaling, intercept specification and restrictions must agree.

### What does `fit = FALSE` mean?

`fit = FALSE` generates the model and data without estimating the model. It is
useful for checking or editing the specification. It does **not** provide fitted
estimates. `fit = TRUE` fits using the requested backend; it does not mean
"fit Stan first, then Julia".

```r
setup <- CoTiMA2::ctmaFit(
  ctmaInitFit = initFit,
  backend = "stan",
  fit = FALSE
)
setup$ctModel
cat(setup$ctModelText)  # editable R specification defining m
# m <- eval(parse(text = setup$ctModelText))
```

The lower-level `ctmaJuliaFit()` can also use an unfitted setup or a compatible
saved fit. That is an advanced workflow; it is not necessary before using the
public backend argument.

### Resource configuration

Use `coresToUse` to state your simultaneous CPU budget. Avoid giving each of
several parallel chains the whole machine. Bundled Julia additionally exposes
`julia_threads`, `subject_threads` and `blas_threads`; those are not native
ctsem controls. A running Julia thread pool may require restarting its process
to change its size.

`ctmaFitConfig()` defaults to the **bundled Julia** route, so specify its
backend explicitly. Its detailed subject/chunk benchmark is specific to that
route. For Stan and native ctsem Julia it validates/resolves controls and
reports setup cost; it does not run the same tuning grid. Native configuration
currently supports ML only. Keep configuration and fitting specifications
identical, and do not transfer a recommended argument list between backends.

## 4. ML, MAP and sampling

For native ctsem Julia and Stan, distinguish:

- `optimize = TRUE, priors = FALSE`: ML.
- `optimize = TRUE, priors = TRUE`: MAP, where supported.
- `optimize = FALSE, priors = TRUE`: posterior sampling with proper priors.

`optimize = FALSE, priors = FALSE` does not select MAP. Bundled Julia supports
ML and Bayesian sampling, but does not support MAP optimization through this
interface.

Example native posterior controls:

```r
bayesFit <- CoTiMA2::ctmaFit(
  ctmaInitFit = initFit,
  backend = "ctsemjulia",
  optimize = FALSE,
  priors = TRUE,
  coresToUse = 2,
  backendArgs = list(
    sampleControl = list(
      sampler = "nuts",
      chains = 4,
      warmup = 1000,
      draws = 1000,
      seed = 1234
    )
  )
)
```

Native ctsem also exposes other sampler choices through `sampleControl`,
including SAEM, learned-map sampling and SMC. Their supported targets and
placement settings differ. Do not force `placement = "saem"` on a marginal
model without sampled random effects. The native `sampler = "map"` option
means learned-map **sampling**, not MAP optimization. Bundled Julia has a
different Bayesian interface and does not accept native `sampleControl`.

ML normally requests Hessian-based uncertainty. That uncertainty can be
unreliable when the Hessian needs substantial repair or parameters are not
identified. The absence of an error is not evidence that intervals are valid.

## 5. Differences between studies versus differences between people

### Study-level fixed effects and random effects

Standard fixed-effects CoTiMA does **not** automatically give every study a
separate drift coefficient. With the current default
`studyVarying = c("CINT", "DIFFUSION", "T0VAR")`, drift is shared across studies,
while eligible parameters in the listed families can have fixed study effects.
Add `"DRIFT"` only when study-specific drift is intended. Fixed study effects
are freely estimated contrasts; hierarchical study random effects instead
model selected study differences as arising from a distribution.

`studyCoding = "dummy"` chooses a reference-study coding.
`studyCoding = "weighted"` uses weighted effect coding. Coding is not itself
a random-effects model, and it does not guarantee better convergence. Under
nonlinear parameter transforms, distinguish a raw-scale center from an average
of transformed study parameters.

### Person-level random intercepts

A random intercept represents stable between-person differences, separated
from within-person dynamics. This is different from study-level random effects.

- `indVarying = "MANIFEST"` or `"CINT"` requests a person-level random-intercept
  structure with a shared population covariance across studies in the pooled
  model.
- `randomIntercepts = "MANIFEST"` or `"CINT"` uses CoTiMA's phantom-process
  formulation and can allow study-specific intercept covariances through study
  effects. It is not interchangeable with the preceding argument.
- Bundled Julia translates supported `indVarying` intercept specifications
  into its compatible phantom formulation; native ctsem Julia can use its own
  random-effect integration.

MANIFEST adds the stable effect to the measurement mean; CINT introduces it
through the dynamic intercept. Equivalence requires an appropriate mapping and
compatible model constraints; the numerical intercept parameters are not
interchangeable.

Phantom states are constant processes with no diffusion. They can enlarge model
matrices. Their names, order and dimension depend on the constructed model:
inspect `setup$ctModel$pars` and latent names rather than universally assuming
that every random-intercept model has exactly twice as many states. Fixing all
variance of a phantom random intercept to zero removes that random intercept.
See the [RI guide](CoTiMA2_RI_Random_Effects.Rnw) for matrix specifications.

## 6. Two distinct routes to hierarchical study random effects

### Native ctsem Julia route

`ctmaFit(backend = "ctsemjulia")` accepts study-level selections in
`backendArgs`, for example `study_random_effects`, `study_parameters`, or matrix
masks `reDRIFT`, `reDIFFUSION`, `reT0VAR`. It uses native ctsem integration and
optimizer/sampler controls. This is separate from `ctmaRePrep()`.

Begin with a few scientifically motivated study-varying parameters. Inspect
study membership and the generated specification before fitting. A mask must
match the actual generated latent matrix, and may select only eligible free
parameters. Native Laplace ML integrates study effects; it must not be confused
with the legacy joint-mode optimizer described below. The accuracy of the
integration and optimizer convergence require separate checks.

### Separate `ctmaRePrep()` / `ctmaReFit()` route

This route modifies the generated Stan model to add hierarchical study effects.
`ctmaReFit()` supports Stan and an initial bundled-Julia Bayesian implementation;
it does not accept `backend = "ctsemjulia"`.

The following example assumes an ordinary two-latent-state setup **without
additional phantom states**:

```r
setupRE <- CoTiMA2::ctmaFit(
  ctmaInitFit = initFit, backend = "stan", fit = FALSE,
  indVarying = FALSE, randomIntercepts = FALSE
)

reDRIFT <- matrix(FALSE, nrow = 2, ncol = 2)
reDRIFT[1, 1] <- TRUE
reDIFFUSION <- matrix(FALSE, nrow = 2, ncol = 2)
reT0VAR <- matrix(FALSE, nrow = 2, ncol = 2)

prepared <- CoTiMA2::ctmaRePrep(
  setupRE,
  reDRIFT = reDRIFT,
  reDIFFUSION = reDIFFUSION,
  reT0VAR = reT0VAR,
  covariance = "independent",
  parameterization = "auto"
)
prepared$specification

reFit <- CoTiMA2::ctmaReFit(
  prepared, backend = "stan", estimation = "bayes",
  chains = 4, iter = 2000, cores = 2
)
reSummary <- summary(reFit)
```

For symmetric covariance parameterizations, select eligible lower-triangle
entries. `"independent"`, `"equicorrelated"`, `"full"` and `"block"` describe
study-effect covariance structures. Rich covariance models need enough studies
and informative data; there is no universal safe minimum study count.

Centered and noncentered formulations describe the same statistical model when
transformed consistently, but their numerical behavior differs. `"auto"` is a
starting rule, not proof that the chosen parameterization is optimal.

The separate route's `estimation = "ml"` (including `ctmaReMLFit()`) maximizes
a **joint/conditional** objective over study effects; it is not marginal ML.
MAP also targets a joint mode with priors. Optimizer uncertainty defaults to
`"none"`; additional experimental uncertainty methods exist, but point
estimates alone provide no intervals. Bayesian fitting is the main inference
route for this separate implementation.

### Which hierarchical summary should I report?

For a `ctmaReFit()` result:

- `population_parameters_raw_centre` transforms the population's raw-scale
  center. Under nonlinear transforms, it is not the transformed-scale mean.
- `population_parameters_study_average` averages transformed parameters over
  the **observed studies with equal weight within each draw**. It is not
  automatically the expected parameter in a new study or a participant-weighted
  average.
- `study_parameters` contains study-specific summaries.
- `study_sd_raw` and `study_sd_transformed` describe heterogeneity, not the
  standard error of a population mean. Read the summary's definition of the
  transformed quantity before choosing what to report.

These field names belong to the separate RE implementation and should not be
assumed to exist on a native ctsem Julia fit. See the
[study RE guide](CoTiMA2_Random_Effects.Rnw) for interpretation, priors and
sensitivity checks. `ctmaReMoments()` provides a post-hoc heterogeneity summary
of an existing fixed-effects fit; it does not refit a hierarchical model or
provide its partial pooling. `ctmaReRedHet()` and `ctmaRedHet()` compare
heterogeneity reductions under their respective supported inputs; a negative
reduction means an increase, not an impossible result.

## 7. Binary outcomes, time scaling and stationarity

Use `binaries = c(1, 0)` when the first variable is binary and the second
continuous. Supply this consistently in primary and pooled fits. Binary data
should use the intended 0/1 coding. The binary likelihood accounts for the
outcome type, but does not guarantee identification of latent scale, diffusion
or intercepts, especially with rare outcomes and few observations per person.

`cint = "auto"` and `manifestMeans = "auto"` are alternative ways to request
appropriate free binary location parameters. Do not free redundant locations
without considering identification. Random-intercept specifications can change
which parameters provide location. Gaussian and binary measurement calculations
can differ between backends; a Stan comparison is not automatically a gold
standard for binary data.

`scaleTime` multiplies the long-data time values: `1/12` divides them by 12.
Drift and diffusion rates change with the time unit, so compare results on the
same time scale. Use the explicitly labelled original-time-scale results when
available. Neither `stationary = FALSE` nor standardized observations imply
stationary latent dynamics. Observed unit variances do not, by themselves,
determine latent diffusion in a model with random intercepts or measurement
error. Do not impose stationary covariance constraints without specifying the
intended within-person covariance.

## 8. Inspect and save results

A short `print()` overview does not mean the native fit lost its contents.

```r
# Primary-study fits remain a list of complete native fit objects.
firstStudyFit <- initFit$studyFitList[[1]]

# A pooled fit stores its single native fit directly in studyFitList.
nativePooledFit <- pooledFit$studyFitList

# Recover a study list, retaining wide data and additionally returning long data.
recovered <- CoTiMA2::ctmaFitToPrep(
  initFit, reUseEmprawData = TRUE, emprawLong = TRUE
)
recovered$emprawList[[1]]      # wide data
recovered$emprawLongList[[1]]  # additional long data in original time units
```

Per-study fits now retain `$empraw`; older fits may retain empirical data only
in the parent object's `$emprawList`. Do not regenerate pseudo-data accidentally
when the purpose is an exact comparison of two engines.

For ML/MAP, inspect convergence, gradients, curvature and identification.
An objective-change stopping rule alone can stop on a nearly flat surface.
A repaired covariance matrix and tiny printed standard errors do not establish
precision. For sampling, inspect divergences, chain mixing, R-hat near 1
(commonly below 1.01), bulk/tail effective sample sizes and Monte Carlo error.
Thresholds are diagnostic aids, not guarantees. A short pilot demonstrates
computational feasibility rather than posterior validity.

`ctmaFitDiagnostics()` reads stored diagnostics without refitting. Summary
methods and table requests can calculate missing summaries or transformations;
consult their messages. `ctmaTables()` supports study/parameter/statistic
sorting through `sortBy`; a requested sorting statistic must be selected in the
table. Version 0.14.51 includes the sorting changes and Stan intercept-metadata repair. See the [tables guide](CoTiMA2_Tables.Rnw).

## 9. A sensible order of work

1. Verify installation and the actually loaded package versions.
2. Inspect the generated model and data with `fit = FALSE`.
3. Fit a small version using an explicit backend and CPU budget.
4. Resolve convergence and identification problems before interpreting effects.
5. Introduce person-level and study-level variation deliberately, one change at
   a time. Do not select a complex RE covariance structure merely because it is
   available.
6. Check sensitivity to starting values, priors and justified parameterizations.
7. Compare backends at matched parameters/model definitions when their answers
   disagree; a faster run is not necessarily a better estimate.
8. Save complete fits and record versions, time units, data provenance and
   uncertainty methods.

For more detail, read the versioned **CoTiMA2 User Guide**, the
[Julia fitting diagnostics guide](CoTiMA2_Julia_Fit-Diagnostics.Rnw), and the
specialist guides linked above.
