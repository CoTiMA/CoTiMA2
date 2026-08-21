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
