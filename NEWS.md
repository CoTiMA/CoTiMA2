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
