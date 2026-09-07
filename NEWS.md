# CoTiMA2 0.14.14

- Added `ctmaStudyCodingBenchmark()` for a controlled, paired comparison of
  dummy and weighted study coding with identical Stan settings, seeds, model
  arguments, and data. It reports runtime per post-warmup iteration,
  divergences, tree-depth saturation, leapfrog effort, acceptance, R-hat, and
  effective sample size while clearly treating short runs as computational
  diagnostics rather than substantive posterior analyses.
- Documented that simulations with identical true effects across studies are
  useful for isolating coding-related computational geometry but do not test
  behavior under genuine between-study heterogeneity.

# CoTiMA2 0.14.13

- Fixed Stan Bayesian fits started through `ctmaFit()`, `ctmaAllInvFit()`, or
  `ctmaInit()` so that an unspecified warmup is set to half of `iter` instead
  of being passed to ctsem as zero warmup. An explicit `backendArgs$warmup`
  override is now supported and validated.
- Added a short study-coding diagnostic under `tools/` to compare dummy and
  weighted design geometry and, optionally, bounded Julia optimization pilots.

# CoTiMA2 0.14.12

- Hierarchical summaries now distinguish the transformed raw-scale population
  centre from the equal-study average of transformed study parameters. The
  latter is calculated within every posterior draw and is normally the
  meta-analytic natural-scale estimand.
- Added equal-study posterior summaries for jointly transformed T0VAR and T0
  covariance elements. Calling `summary()` prints a concise interpretation of
  the different population, study-specific, and heterogeneity estimands.

# CoTiMA2 0.14.11

- Added backend-independent, parameter-selective primary-study contrasts via
  `studyCoding` and `studyVarying`. The defaults use dummy coding and allow
  `CINT`, `DIFFUSION`, and `T0VAR` to differ between studies while keeping
  DRIFT invariant unless explicitly requested.
- Deprecated `WEC=TRUE`; it now maps to weighted study coding for `CINT`,
  `DIFFUSION`, and `T0VAR`. Weighted DRIFT heterogeneity analyses must include
  `"DRIFT"` explicitly in `studyVarying`.
- Updated `ctmaRedHet()` validation and the user and Julia-diagnostics guides
  for the new study-coding interface.

# CoTiMA2 0.14.10

- Added `ctmaTransformJuliaFit()` to calculate constrained ctsem and
  study-specific point transformations from a completed Julia fit without
  optimization, a Hessian, covariance estimation, or uncertainty draws.
- Added the standalone `CoTiMA_Julia_Fit-Diagnostics.Rnw` guide, including a
  decision tree for parallel-backend, optimizer, study-design, local-parameter,
  transformed-point, and Hessian diagnostics.
- `ctAddJuliaSamples()` and `ctmaTransformJuliaFit()` now also accept diagnostic
  preparation objects that store the matching ctsem setup under `$fit`.

# CoTiMA2 0.14.9

- Added a guarded Julia ML subset warm start. With
  `backendArgs$warm_start = "auto"`, CoTiMA2 samples complete persons within
  studies, enlarges the requested fraction for base-parameter count,
  repeated-measurement, study-coverage, and binary-category requirements, and
  optimizes base parameters while TI effects remain zero.
- A subset candidate is accepted only when its full-data objective is finite
  and no worse than the ordinary start. Explicit `inits` disable the pilot.
  Detailed selection, optimization, validation, and timing diagnostics are
  retained in `fit$warm_start` and `fit$warm_start_timing`.

# CoTiMA2 0.14.8

- `ctmaFit(backend="julia", indVarying="CINT")` and its `"MANIFEST"`
  counterpart now use a Julia-compatible augmented phantom-process model
  instead of ctsem's unsupported state-dependent `nindvarying` branch.
- The restricted model freely estimates one pooled person-level
  random-intercept covariance matrix while fixing all TI effects on phantom
  variances and T0--phantom covariances. Phantom means and their TI effects
  remain available. The returned argument list records
  `randomInterceptStructure = "pooled"`, providing an explicit extension point
  for later study-random-effects or Laplace estimation.

# CoTiMA2 0.14.7

- Removed captured scalar and array aliases from threaded likelihood and
  gradient reductions. Each worker now updates only its explicitly indexed
  partial accumulator. Subject-level diagnostics established that the
  numerical subject kernel itself was exactly reproducible and localized the
  race to this reduction layer.

# CoTiMA2 0.14.6

- Subject-parallel likelihood and gradient loops now use static thread
  assignment. A new subject-level consistency diagnostic localizes any
  remaining serial-versus-parallel discrepancy.
- ctmaFitConfig() no longer ranks thread counts when both concurrent
  components have fallen back to their deterministic serial implementations.
  Candidate timings are repeated and ranked by their medians.
- Short carefulFit pilots now consider raw-gradient reduction in addition to
  objective improvement per second, making their recommendation more relevant
  to eventual convergence.

# CoTiMA2 0.14.5

- Added an independent repeatability and serial-agreement gate for concurrent
  numeric likelihood evaluation. Likelihood and gradient parallelism are now
  enabled separately only when their own complete-model checks pass.
- Replaced spawned likelihood reductions with bounded, indexed
  `Threads.@threads` partial sums. If either parallel component is unstable,
  the optimizer automatically uses its deterministic serial counterpart, so
  the objective and gradient always describe the same function.

# CoTiMA2 0.14.4

- Reintroduced bounded parallel subject gradients using indexed
  `Threads.@threads` reductions instead of spawned task closures. Before ML or
  Bayesian estimation, Julia now compares the complete serial and parallel
  gradients at the common starting vector. Parallel AD is enabled only when
  the maximum relative coordinate difference is at most `1e-10`; otherwise it
  automatically retains the exact serial gradient while keeping numeric
  likelihood evaluations parallel.

# CoTiMA2 0.14.3

- Julia now uses subject-task parallelism only for ordinary numeric likelihood
  evaluations. ForwardDiff likelihoods, gradients, and Hessians use the
  serial-safe path because concurrent dual-number evaluation was shown to
  disagree with the exact serial derivative on TI-predictor models.
- Setup-only `ctmaJuliaFit(estimation="none")` calls no longer calculate
  Hessian uncertainty unless `uncertainty=TRUE` is explicitly requested.

# CoTiMA2 0.14.2

- Fixed concurrent subject reductions so every spawned Julia task captures
  its own immutable subject range and owns its likelihood, gradient, or
  Hessian accumulator. This restores agreement between serial and concurrent
  differentiation for TI-predictor models.

# CoTiMA2 0.14.1

- Julia ML optimization now retains the best finite objective evaluated in
  every optimization stage. A stalled `ObjectiveIncreased` or `NoXChange`
  endpoint can no longer replace a better accepted point or be propagated to
  the next stage.
- Joint and staged Julia optimization automatically retry stalled stages with
  initial steps reduced by factors of ten, while reporting the total optimizer
  work across attempts.

# CoTiMA2 0.14.0

- Added `ctmaFitConfig()`, a time-budgeted configuration pilot accepting the
  model-building arguments of `ctmaFit()`. Julia pilots rank combinations of
  subject concurrency and ForwardDiff chunk size, then compare short joint and
  staged (`carefulFit`) optimizer pilots. Stan controls are validated and
  retained with setup timing reported.
- Julia thread-pool capacity now defaults to `julia_threads = "auto"`, resolved
  as the CPUs available to the R process minus one. The new runtime-adjustable
  `subject_threads` control limits concurrent subject-processing tasks without
  requiring an R restart.
- Replaced thread-ID-indexed static subject reductions with independent
  task-local reductions for likelihoods, gradients, study-RE gradients, and
  Hessians.

# CoTiMA2 0.13.8

- Fixed the bundled Julia environment by resolving `ForwardDiff` and `Optim`
  as direct dependencies and added an actionable setup error if required Julia
  packages are unavailable.

# CoTiMA2 0.13.7

- Added `Julia_Install.Rnw`, a beginner-oriented installation guide with
  separate copy-and-paste instructions for Windows, macOS, and Ubuntu/Linux,
  including R upgrades, JuliaCall rebuilding, verification, and troubleshooting.
- Julia backend requests now print the R command that opens the Julia
  installation and troubleshooting guide.
- Simplified the internal ctsem setup notice, suppressed duplicate version
  notices from nested CoTiMA calls, and silenced Julia's project-activation
  output.
- Declared `ForwardDiff` and `Optim` as direct dependencies of the bundled
  Julia project so that `ctmaJuliaSetup()` installs everything loaded by
  `CtsemJulia.jl`.

# CoTiMA2 0.13.6

- Added `CoTiMA_RI_Random_EFfects.Rnw` as the fourth maintained vignette,
  documenting random-intercept CoTiMA with study-level random effects.

# CoTiMA2 0.13.5

- Prevented extreme binary-logit predictions from rounding to exact zero or
  one and producing a singular Bernoulli innovation variance.
- Julia's staged optimizer now treats positive-definiteness, singularity, and
  threaded task failures in gradient callbacks as rejectable line-search trial
  points instead of aborting the complete fit.
- Renamed the Julia and random-effects guide sources to
  `CoTiMA_Julia_Usage.Rnw` and `CoTiMA_Random_Effects.Rnw`.

# CoTiMA2 0.13.4

- Consolidated the combined hierarchical and moment-based random-effects guide
  into the canonical package vignette directory.
- Added the general CoTiMA 1.0.4 guide and its figure assets to the vignette
  sources, yielding three maintained Rnw guides: general use, Julia, and
  random effects.
- Added the small random-effects recovery simulation under
  `tools/random-effects-validation/` and adapted it to load CoTiMA2 directly.
- The local build helper now compiles every Rnw source and installs its PDF
  under `inst/doc` before building the source tarball.

# CoTiMA2 0.13.3

- Consolidated the Julia guide into one canonical package vignette, removed
  duplicate/generated vignette build files, and added a reproducible local
  build-and-install verification script.
- Updated the guide's integration-status section to describe the Julia backend
  as an existing package feature rather than future work.

# CoTiMA2 0.13.2

- Expanded the Julia usage vignette with binary-variable syntax, CINT versus
  MANIFESTMEANS parameterizations, their identification rationale, validation
  rules and exemptions, and details of the Julia binary likelihood.

# CoTiMA2 0.13.1

- Added Julia estimation for mixed continuous and binary manifest-variable
  models. The Julia likelihood now includes binary-logit probabilities,
  Bernoulli log-likelihood contributions, Bernoulli observation variance in
  the extended Kalman update, and ctsem-compatible measurement Jacobians.

# CoTiMA2 0.13.0

- Added model-specification checks for binary manifest variables. Every
  binary variable must have either a freely estimated connected CINT or a
  freely estimated MANIFESTMEANS parameter; specifications with neither now
  stop before fitting. CINT and manifest-mean vectors are also checked against
  their latent and manifest dimensions.

# CoTiMA2 0.12.9

- Julia maximum-likelihood fits now estimate Hessian-based uncertainty and
  transformed standard errors and confidence intervals by default, including
  direct `ctmaJuliaFit()` calls. Stan continues to estimate uncertainty by
  default through `optimcontrol$estonly = FALSE` and positive
  `optimcontrol$finishsamples`.

# CoTiMA2 0.12.8

- Restored the Julia setup and executable-discovery source in source-package
  builds so automatic Julia discovery works on macOS, Linux, and Windows.

# CoTiMA2 0.12.7

- Increased the default Julia ML optimization ceiling from 100 to 5,000
  iterations for direct `ctmaJuliaFit()` calls and higher-level workflows.

# CoTiMA2 0.12.6

- Julia fits now print the effective Julia and BLAS thread counts.
- `ctmaFit(backend = "julia")` now automatically translates `coresToUse`
  into a Julia resource configuration when users have not explicitly supplied
  Julia worker/thread controls. ML uses the resolved cores for subject
  parallelism (or divides them across multistarts), Bayesian estimation divides
  them across concurrent chains, and BLAS remains single-threaded to avoid
  nested oversubscription.
- Audited the other Julia fitting entry points. `ctmaInit()`,
  `ctmaAllInvFit()`, `ctmaEqual()`, and `ctmaPower()` now identify themselves
  correctly in their automatic-allocation notes; `ctmaOptimizeFit()` now also
  derives its Julia allocation from `coresToUse`; and `ctmaReFit()` reports its
  effective main-process Julia and BLAS thread counts.

# CoTiMA2 0.12.5

- Suppressed the misleading "No model is fitted" banner when `ctmaFit()`
  internally prepares a `fit = FALSE` ctsem model for immediate Julia
  estimation. The banner remains visible when users explicitly request
  `fit = FALSE`.

# CoTiMA2 0.12.4

- `ctmaGenData(ctmaExtract = FALSE)` now prints a follow-up
  `ctmaExtract()` command containing the actual `activeDirectory`. The new
  optional `resultName` argument inserts the user's result-object name; when
  it is omitted, the note clearly identifies the placeholder that must be
  replaced.

# CoTiMA2 0.12.3

- Improved conflict handling in `ctmaLoadExamples()`. Its default error now
  explains how to remove or rename conflicting objects or deliberately use
  `overwrite = TRUE`. Deliberately overwritten values are additionally saved
  in `ctmaSecuredObjects`, including distinct backups across repeated loads.

# CoTiMA2 0.12.2

- Fixed false modification reports in `ctmaUnloadExamples()` for complex fit
  objects containing environments, external references, or lazily evaluated
  components. Loaded values are now retained by reference and compared with
  `identical()`, relying on R's copy-on-modify behavior. Added an integration
  test that loads every packaged example, changes one object, and verifies
  that only this object remains.

# CoTiMA2 0.12.1

- Replaced all internal dependence on the exported `CoTiMAStanctArgs` object
  with a private factory that returns a fresh Stan-default list for every
  call. `ctmaBackendDefaults("stan")` is now the canonical public accessor.
- Retained `CoTiMAStanctArgs` temporarily as a deprecated compatibility
  snapshot and retained the legacy function arguments with deprecation
  notices. Removed its redundant user-guide example-data copy.

# CoTiMA2 0.12.0

- Migrated `ctmaRedHet(dt=)` to a backend-neutral raw-draw map. Discrete-time
  residual effects are now calculated draw by draw as the difference between
  matrix exponentials of the complete study DRIFT and its moderator-predicted
  DRIFT, before estimating residual heterogeneity.
- Packaged user-guide objects are no longer attached automatically. Use
  `ctmaLoadExamples()` to load selected examples explicitly and
  `ctmaUnloadExamples()` to remove only examples that remain exactly
  unchanged. `ctmaPrep()` now constrains all name-based lookups to its caller.

# CoTiMA2 0.11.1

- Migrated the continuous-time `ctmaRedHet()` path to the corrected shared
  two-stage meta-analysis calculations used by `ctmaBiG()`.
- Heterogeneity reduction is now based on WEC-reconstructed residual study
  effects around moderator-predicted study means. The aligned result table
  reports unconditional and moderated residual variance and proportional
  reduction, while retaining Q, I-squared, and descriptive dispersion as
  secondary diagnostics.
- Added backend-neutral extraction, explicit handling of zero denominators and
  negative reductions, S3 print/summary methods, and validation tests.

# CoTiMA2 0.11.0

- Added ctmaReRedHet() to compare residual study-level random-effect
  variances before and after including moderators. It reports
  parameter-specific pseudo-R-squared reductions, independent cross-posterior
  Monte Carlo intervals, threshold probabilities, and an overall reduction in
  the sum of marginal variances.
- Added strict comparability checks for studies, random-effect
  specifications, transformations, covariance structures, and priors.
- Documented the moderator heterogeneity-reduction workflow and its
  interpretation limits in the Julia/CoTiMA methodology vignette.

# CoTiMA2 0.10.2

- Corrected `ctmaPub()` author matching across common `Surname, Initials` and
  `Initials; Surname` formats and stopped treating labels such as `Study1` as
  authors.
- Corrected recency aggregation for functions other than `sum`, and corrected
  exclusion of individual publication years.
- Added input validation and returned author-matching diagnostics.

# CoTiMA2 0.10.1

- Added backend-independent primary-study extraction to `ctmaBiG()`.
  Transformed results from Stan and Julia `ctmaInit()` fits now enter the same
  validated meta-analysis implementation, including mixed-backend objects.
- `ctmaBiG()` now records the detected backend for every included study and
  gives an actionable error when a Julia fit lacks transformed draws.

# CoTiMA2 0.10.0

- Corrected `ctmaBiG()` common- and random-effects calculations: two-sided
  p-values, non-negative DerSimonian--Laird heterogeneity, parameter-specific
  study counts, corrected H-squared/I-squared intervals, and prediction
  intervals that include uncertainty in the pooled mean.
- Corrected T0 covariance standard-error extraction and discrete-time z-curve
  input selection in `ctmaBiG()`.

# CoTiMA2 0.9.99

- Added `backend` and `backendArgs` to `ctmaEqual()`, with Stan retained as
  the default and Julia ML/Bayesian fitting routed through the shared backend
  infrastructure.
- `ctmaEqual()` now accepts both conventional `CoTiMAFit` inputs and the
  compact `ctmaJuliaFit` objects returned by Julia-based `ctmaFit()` calls.
- Equality-constrained study-parameter summaries now use unique display row
  names when several matrix cells share one parameter.
- Fixed `ctmaEqual(NULL)` so that it reports the intended input error instead
  of dereferencing the missing object first.

# CoTiMA2 0.9.98

- Added `backend` and `backendArgs` to `ctmaPower()`. Its internal invariant
  model now uses either Stan (default) or Julia consistently, including
  resource allocation, finish-sample mapping, and stored backend metadata.
- Fixed one-variable long-data filtering in `ctmaPower()`.

# CoTiMA2 0.9.97

- Added a central backend-aware core resolver. Stan/R worker counts are
  conservatively capped to avoid exhausting R connections on many-core hosts;
  Julia keeps its thread-based allocation.
- Added the public `backend` and `backendArgs` interface and Julia ML/Bayesian
  route to `ctmaAllInvFit()`; Stan remains the default backend.

# CoTiMA2 0.9.96

- `ctmaInit(finishsamples=...)` now maps to Julia's `transformed_draws`
  unless `backendArgs$transformed_draws` is supplied explicitly.

# CoTiMA2 0.9.95

- Resolves negative `coresToUse` values before allocating Julia resources;
  for example, `coresToUse = -1` uses every detected core except one. ML,
  multistart ML, and Bayesian chains receive an automatic allocation unless
  the user supplies an explicit Julia thread control in `backendArgs`.
- Labels the elapsed time reported on exit from each principal workflow as its
  total wall-clock execution time, distinguishing it from nested fit,
  transformation, or per-study timings.

# CoTiMA2 0.9.94

- Maps `indVarying = "MANIFEST"` and `indVarying = "CINT"` to CoTiMA's
  Julia-compatible latent pseudo-process random-intercept formulation within
  `ctmaInit()` only. When not overridden in `backendArgs`, `coresToUse` now
  supplies `julia_threads` for this workflow. Julia ML single-study fits use
  5,000 optimizer iterations; automatically translated pseudo-process fits use
  a 0.05 raw-gradient tolerance unless the user supplies another value.
- Fixes pseudo-process construction for one or any other number of latent
  variables and retains a ctsem transformation object so Julia fits can feed
  CoTiMA's random-intercept covariance summaries.

# CoTiMA2 0.9.93

- Adds the experimental Julia ML and Bayesian backend route to `ctmaInit()`,
  retaining ctsem-transformed single-study summaries for downstream CoTiMA2
  workflows.

# CoTiMA2 0.9.92

- Reports wall-clock execution time when each principal workflow exits,
  including exits caused by an error.

# CoTiMA2 0.9.91

- Establishes the experimental CoTiMA 2.0 development line under the separate
  R package name `CoTiMA2`, allowing the original `CoTiMA` package to remain
  installed independently.
- Adds the optional Julia ML and Bayesian backends, study-level random-effects
  workflows, backend argument infrastructure, and the associated validation
  and documentation developed before this preview snapshot.
- Principal modeling workflows print an explicit experimental-version note at
  entry so saved scripts and console logs identify which implementation ran.
