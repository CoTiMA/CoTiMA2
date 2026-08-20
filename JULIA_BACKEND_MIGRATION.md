# CoTiMA backend migration

This file records the incremental migration from the legacy
`CoTiMAStanctArgs` structure to a common backend interface. Every stage should
be a narrow commit for which the package tests pass and the source package
builds. Stan remains the default backend throughout the migration.

| Stage | Scope | Status |
|---|---|---|
| 1 | Add backend defaults, recursive validation, merge helpers, and tests without changing fitting behavior | Complete |
| 2 | Replace fragile legacy Stan-list merges internally, one function at a time | Complete for `ctmaInit()`, `ctmaFit()`, `ctmaAllInvFit()`, `ctmaPower()`, and `ctmaOptimizeFit()` |
| 3 | Add `backend` and `backendArgs` to `ctmaFit()` while retaining all old calls | Complete; Stan remains default and Julia requests identify the staged boundary |
| 4 | Route validated Julia ML/Bayesian estimation through `ctmaFit()` | Complete |
| 5 | Repair `ctmaOptimizeFit()` workflow issues and make it backend-aware | Workflow correctness repaired; backend routing pending |
| 6 | Migrate additional CoTiMA fitting functions individually | Pending |
| 7 | Deprecate `CoTiMAStanctArgs` with a compatibility translator | Pending |
| 8 | Remove the duplicate data object only after the transition period | Pending |

## Argument ownership

- Top-level CoTiMA workflow/model controls: `scaleTime`, `randomScaleTime`,
  `scaleTI`, `randomScaleTI`, `reFits`, and similar settings that change the
  model or prepared data.
- Backend controls: compiler, optimizer, sampler, threading, and diagnostic
  settings supplied through `backendArgs`.
- Outer repeated fits (`reFits`) rebuild and refit CoTiMA models. Inner Julia
  `multistart` attempts optimize one fixed prepared likelihood from multiple
  starting values; the two concepts remain distinct.

## Compatibility rules

1. `backend = "stan"` is always the default.
2. Existing calls using `CoTiMAStanctArgs` remain supported during migration.
3. If legacy and new backend controls are supplied together, the eventual
   public interface must report a clear conflict instead of guessing.
4. Unknown backend control names are errors; nested lists merge recursively.
5. Julia options are exposed only after their workflow has been validated.
