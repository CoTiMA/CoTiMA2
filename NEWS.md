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
