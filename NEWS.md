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
