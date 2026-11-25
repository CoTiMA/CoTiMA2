# CoTiMA2 0.14.70

- Added ctmaFitToModel with original, init, fit and prep conversion targets.
- Added ctmaInit(model=) and reusable-model support in ctmaFit(model=).
- Original long-data times are reconstructed and the requested time factor is applied once.
- Explicit arguments override inherited settings, including NULL and FALSE.
- Clarified WEC-centred baseline drift matrix terminology.

# CoTiMA2 0.14.69.9001

- Fix ctmaRePrep compatibility with ctsem 4 additive T0 covariance and its Julia-default setup routing. Preserve legacy Stan templates.

# CoTiMA2 0.14.69

- Enable native Julia Bayesian/MAP estimation for determined diffusions. Priors are defined on the stationary raw coordinates, not transported from unrestricted diffusion priors.
- Map study-level DRIFT random effects to the shared auxiliary coordinates while retaining saturated stationary covariance contrasts. Reject independent random effects on derived diffusion/initial covariance and person-level process coordinates.
- Preserve ML behavior and the fixed stationary-variance option. Record the prior-coordinate interpretation in fit metadata.

# CoTiMA2 0.14.64

* Support named location specifications and selective random intercepts in ctmaInit and ctmaFit. Warn about binary indicators with all connected location terms fixed at zero.

# CoTiMA2 0.14.62

* Distribute primary-study fits across local and SSH workers with `ctmaInit(backend = "ctsemjulia", compute = ctmaCompute(tasks = "studies"))`.

# CoTiMA2 0.14.60 (2026-11-13)

* Adds experimental `compute` plans to `ctmaFit()`, `ctmaInit()` and `ctmaJuliaFit()` for bundled Julia maximum likelihood: distribute one shared likelihood/gradient calculation or independent optimization starts.
* Adds `ctmaHost()` and `ctmaCompute()` with per-host worker/thread budgets, multiple SSH hosts, and an external coordinator with zero workers on the submitting laptop.
* Adds SSH connectivity checks, isolated Git source snapshots, checkpoint recovery, and persistent background jobs with `ctmaJobStatus()`, `ctmaCollect()` and `ctmaResume()`. Passwords are not stored; unattended multi-host jobs require coordinator-side authentication.
* Documents examples, setup, limitations and measured Mac/Ryzen scaling in the versioned User Guide. Distributed NUTS chains are outside scope; ctmaInit study scheduling remains sequential. Integrated subject/study random effects and priors are not supported by this distributed ML path.

# CoTiMA2 0.14.59

* Includes compact initialization notes and explicit fitted/original summary time units.
* Adds internal Lohmann FE/URE/ERE/CRE replication and a separately labelled univariate extension, with two preparation/sampling example scripts and standalone documentation.

## Development after 0.14.58: explicit summary time scales

- Conventional and pooled Julia summaries now expose fitted/original-time aliases and label the time convention. Native pooled DRIFT and DIFFUSIONcov estimates, SD/SEs and intervals are converted using the stored scaleTime; stored legacy original-scale tables are retained. Raw parameters, likelihood, z-values, ESS and R-hat are not multiplied. Native discrete-time coefficients remain at the requested fitted-time interval; the corresponding original interval is recorded. No refitting is required.

## Development after 0.14.58: compact initialization diagnostics

- Initialization identification notes now print a study-indexed table and one explanation per flag. Full original messages and parameter names are retained; use `print(fit$summary$identification_notes, details=TRUE)` or `as.character(...)`. During a single ctmaInit call, only the first two-wave warning prints its full explanation; subsequent warnings retain the study label and reference that explanation. No fits, masks or exclusions are changed.

## CoTiMA2 0.14.58 (2026-11-05)

### Study-order reproduction

- `ctmaFit(primaryStudyList=...)` now selects and orders initialization records by original study ID, keeping fits, caches, moderators and metadata aligned. Invalid or ambiguous IDs are rejected before fitting. This corrects the shuffled-optimizer raw-start handoff; raw starts still require the same model and parameter layout.

### Removed legacy controls

- Removed `experimental` and `activateRPB` from all CoTiMA2 function interfaces, internal calls and help. Behavior follows their former `FALSE` defaults; remove these obsolete arguments from existing scripts. Saved older fits can still be read. Pushbullet notifications and their package dependency were removed.

- Targeted installation tests and updated versioned Rnw guides/PDFs accompany this release. The RI reproduction does not validate projected Hessian uncertainty.

## CoTiMA2 0.14.57

- Fixed partial drift invariance: explicit invariantDrift requests enable study variation in the remaining drift effects when studyVarying is omitted; contradictory explicit settings are rejected.
- ctmaEqual rejects incompatible transformations and uses fresh starting values for its changed parameter layout.
- Corrected Inv23 and equality refits completed successfully; updated Rnw guides and PDFs.

## Development fixes after 0.14.55

- Restore native pooled common-drift plots from transformed matrices, retain time scaling, and reject conditional drift plots explicitly.
- Clarify User Guide prerequisites/backend scope and repair wrapped comments, obsolete arguments, a missing comma and LCS save code.

# CoTiMA2 0.14.55 (2026-11-04)

- Fix the 10000-byte R-name limit in determined-diffusion matrix setup with realistic parameter labels. Direct covariance expressions replace redundant Cholesky expansion without changing parameter coordinates or covariance constraints.
- Add exact matrix-setup regression tests for default labels and up to 45 study contexts, plus Lyapunov-identity and native fixed-point/gradient checks.

# CoTiMA2 0.14.54 (2026-11-04)

- Study-design partitions, identification-aware uncertainty reporting and native study-specific determined diffusion are included in this release. See the detailed entries below.

## CoTiMA2 0.14.53

- Fixes explicit indVaryingT0=FALSE being ignored in ctmaInit() and ctmaFit(). Only individual T0MEANS random effects are disabled; manifest/CINT random effects and within-person T0VAR are retained. NULL/TRUE preserve the prior automatic behavior.
- Adds setup-only regression tests for both public constructors and updates help/Rnw/PDF documentation. This fixes the documented stationary manifest-RI determined-diffusion specification without relaxing its unsupported-model guards.

## CoTiMA2 0.14.52

- Adds opt-in randomDeterminedDiffusions=TRUE to the gambling pooled ML search. Stationarity, model, scales and study order remain fixed; both parametrizations are prepared before fitting. Variant choices, R seeds, failures and convergence diagnostics are recorded; raw starting vectors are discarded.
- Updates argument help and versioned Rnw/PDF guides. Determined diffusion retains the original drift transforms and native z random-intercept safeguards. Unsupported covariance constructions are rejected.

## Development after 0.14.51

- Adds explicit stationary=FALSE and determinedDiffusions=FALSE arguments to ctmaInit() and ctmaFit(); existing default models are unchanged.
- The opt-in determined option estimates the stationary within-person covariance for supported one/two-process Gaussian ML models. It derives valid diffusion and preserves the original drift transforms and native random-intercept correlation safeguards. Unsupported constraints, priors, binary data and study-effect mappings are rejected explicitly.
- Provides an exact explicit-initial-covariance compatibility path for zero-CINT determined models on published ctsem 3.11.1, whose stationary switch is disabled. The native ctsemjulia route forwards stationary directly.
- Documents the restrictions in help and new English Rnw/PDF guides. Generated model text preserves auxiliary parameter declarations and covariance construction.

- Report unassessed weak-direction counts as NA and expose curvature_assessed. Legacy summary zeroes without recorded curvature provenance are treated as unknown; retained objects are not modified.

## CoTiMA2 0.14.51

- Add planOnly=TRUE for reviewing candidate nominations and case plans without fitting, likelihood evaluation or checkpoint writes.

- Add automatic reference and exclusion-candidate nominations with stored reasons; no refitting is performed during candidate analysis and exclusions remain explicit.

- Keep `ctmaOptimizeFit(gambling = TRUE)` as the unchanged legacy search; add an opt-in `gambling = FALSE` pooled ML workflow for fixed-model reference checks and explicitly nominated study exclusions.
- Add bounded pilots, per-case diagnostics, owned worker cleanup and fingerprint-checked checkpoints. Subset likelihoods never compete with the full-data fit; short pilots are not final inference.
- Document supported Stan/native ctsemjulia models and systematic examples in the User Guide, Julia Usage and Fit-Diagnostics Rnw files.

## CoTiMA2 0.14.50

- Preserve the native pooled Julia fit in `studyFitList` and delegate standard fit methods.
- Retain empirical wide data and study metadata in each primary-study fit.
- Repair `ctmaFitToPrep()` metadata/data recovery and explicit cached-data reuse.
- Add `emprawLong = FALSE`; opt in to additional long data in `emprawLongList` without replacing wide data or changing time units.

# CoTiMA2 0.14.49

- Match CPU control names exactly so the native optimcontrol$coreset logical flag is not mistaken for an optimcontrol$cores count. This fixes a pre-fit failure discovered in the isolated installation test.

# CoTiMA2 0.14.48

- Resource resolution honors available CPU restrictions and the requested core budget, including bundled Julia worker/BLAS products. Invalid or excessive explicit allocations fail before starting Julia.
- Owned subprocesses inherit single-thread BLAS/OpenMP limits, with the caller's environment restored on exit. This does not reconfigure an R BLAS library already loaded in the caller.
- Bundled Julia setup and fitting use a content-keyed writable user project instead of the installed package tree. Dry runs do not copy the project; dependency setup requires explicit consent.

# CoTiMA2 0.14.47

- Native ctsemjulia supports explicit quadrature refinement for optimized random-effects fits.
- Fit diagnostics distinguish checked and Laplace log likelihoods and expose source-reported optimizer resolution metadata.
- Quadrature retains integrated random-effects diagnostic safeguards; sampling and models without random effects reject this optimization-only request early.

# CoTiMA2 0.14.46

* Automatic Stan study-RE initialization follows the fitting seed (including a
  fit_args$seed override), separates chain streams and restores caller RNG state.
* ctmaOptimizeInit validates required inputs before retry execution and uses
  owned Unix fork workers without replacing the caller's foreach backend.
* ctmaGenData composes cluster cleanup with option restoration on exit.
* Existing output files now cause a clear noninteractive error unless
  silentOverwrite=TRUE; interactive rename choices recheck existing targets.

# CoTiMA2 0.14.45

* Reports cached Laplace approximation-accuracy status separately from optimizer
  convergence in ctmaFitDiagnostics(), native summary footers and table notes.
* Preserves unresolved or missing accuracy verdicts: an operational status of
  "checked" is not interpreted as an accuracy pass. No check or fit is rerun.
* Reports stored integration gap, maximum standardized shift and correction
  application when available, without altering parameter estimates or intervals.

# CoTiMA2 0.14.44

* Corrects mixed per-manifest mean specifications such as `c("mmV1", 0)`;
  explicit vectors no longer fail a scalar `if` condition or lose fixed entries.
* Adds regression coverage for vectors, numeric fixed values, scalar zero and auto.
* Documents the distinction between a binary measurement offset and a dynamic
  CINT, including why an initial mean fixed at zero does not encode a nonzero
  binary baseline prevalence. This is model guidance, not automatic refitting.

# CoTiMA2 0.14.43

* Reports stored bundled NUTS Rhat, bulk ESS and divergences in ctmaFitDiagnostics without refitting.

* Rejects unsupported bundled Julia MAP requests instead of silently fitting ML.
  Use `ctsemjulia` or `stan` with `optimize=TRUE, priors=TRUE` for MAP.
* Bayesian bundled Julia uncertainty uses retained posterior draws without
  requiring a likelihood Hessian at the initialization point.
* Filters removed inert Stan optimizer defaults only when the installed ctsem
  rejects them; explicit unsupported non-default requests produce an error.
* Routes pooled Stan fitting through the compatibility wrapper.
* Supports ctsem 4's renamed matrix registry in transformed study summaries.

# CoTiMA2 0.14.42

* Implements ctsem 4 z covariance conversion in the bundled Julia engine; previous off-diagonal code-2 fits need refitting. Unknown future codes fail before engine startup.
* Accepts native sampleControl$iter without changing its semantics.
* Separates posterior convergence from placement status, handles native pooled storage, and reports stored transformed diagnostics and SMC's non-applicable Rhat correctly.
* Corrects table convergence notes and reuses retained native matrix summaries.
* Clarifies the explicit Hessian default and native runtime selection.

# CoTiMA2 0.14.41

* Table notes now distinguish posterior convergence from optimizer placement, flag stored transformed Rhat, and use native resultsSummary caches.

* `ctmaTables()` now accepts `fontSize = 10` (points) for titles, headers,
  cells and notes in Word, XML, TeX and PDF exports.

## Regular Stan table compatibility (2026-10-29)

* Reads joint population DRIFT arrays from regular optimized and sampled Stan fits.
* Reconstructs missing optimized-fit caches from retained raw draws through the
  fitted Stan program, with Kalman filtering disabled and no model refitting.
* Correctly labels native ML, MAP and Bayesian fits from their direct args fields.
* Retains native-summary precision for four-decimal tables, avoids partial cache
  name matching, and preserves a wrapper's CT representation during DT refresh.
* Validated against regular ctsem 3.11.1 and an existing 2000-draw CoTiMA benchmark.

# CoTiMA2 0.14.40

## Table presentation options (2026-10-29)

* layout selects portrait (long, default) or landscape (wide) pages consistently
  across Word, XML, TeX and PDF exports. Oversized tables continue in column panels.
* sortBy groups ctmaInit estimates by study (default) or by parameter.
* rowLayout keeps one parameter per row (default), combines DRIFT entries, or
  combines all selected parameters with grouped parameter/statistic headers.
* Numeric long-format values and safe staged output publication are preserved.

# CoTiMA2 0.14.39

## Saved-fit table summaries (2026-10-29)

* ctmaTables defaults to four display decimals.
* Missing statistics and DT uncertainty can be regenerated from saved native fits
  with summary(timeinterval=lag), with one consent for all primary studies.
* recomputeSummary="yes" permits batch use; no models are refitted and raw point
  estimates are never treated as uncertainty draws.

# CoTiMA2 0.14.38

## Table selectors and primary-study summaries (2026-10-28)

* Accepts bare, case-insensitive matrix families such as parameters=c(drift, diffusion).
  Quoted selectors and character-vector variables remain supported.
* Presents every primary study from a ctmaInit aggregate summary with original study
  labels, retained uncertainty and explicit aggregate fit-statistic notes; no refitting.
* Preserves unavailable statistics as unavailable and aligns study/parameter text
  columns in Word and TeX exports.

# CoTiMA2 0.14.37

## APA tables (2026-10-28)

* Adds ctmaTables() with checked summary agreement, named/matrix parameter selection,
  source-aware statistics, interactive selection and APA-style Word/XML/TeX/PDF exports.
* Adds joint-draw DT effects at selected fitted-time lags and coordinate-aligned model comparisons.
* Exports are staged; existing outputs require explicit overwrite. Unsupported summaries
  and unavailable statistics fail explicitly rather than fabricating transformed values.
* Documents optional Word/PDF dependencies and the limits of legacy summary provenance.

# CoTiMA2 0.14.36

## Power analysis refactoring (2026-10-28)

* Separates validation, invariant fitting/loading, matrix extraction, implied moments,
  regression effects, power probabilities, integer sample-size search and result assembly.
* Verifies unrestricted Gaussian regression eligibility for direct covariance calculations.
  Supports explicit restricted Gaussian regression syntax through a lavaan fallback, with
  console and stored summary notes identifying the calculation route.
* Corrects negative-effect p-values, two-wave and higher-dimensional calculations, lag
  bounds, reduced hypotheses, root boundaries and tied optimal-lag selection.
* Fast default power includes both rejection tails. useMBESS=TRUE retains the legacy
  probability convention. Missing post-hoc powers stay NA instead of becoming zero.
* Restores original time units from invariant-fit scaling metadata. Reduced-fit caches
  now validate covariance/model metadata and are reused across requested power levels.

# CoTiMA2 0.14.35

## Plotting corrections (2026-10-28)

* Correct forest pooled confidence intervals and funnel confidence limits.
* Correct time-range indexing, axis positions, power-only plots and fractional-lag selection.
* Handle multiple fits, legacy example objects and per-fit moderator/lag mappings.
* Export uniquely named PNG files in the requested directory, including headless sessions.

# CoTiMA2 0.14.34

- Separates optimizer, identification, transform-saturation and posterior
  diagnostics in Julia fit summaries. Adds ctmaFitDiagnostics() for stored
  fits without refitting or starting Julia. Unknown diagnostics remain unknown.
- New native fits cache population matrices for a labelled heuristic flag
  when binary linear-predictor SD exceeds 10; no parameter is fixed or changed.
  Optimizer convergence is explicitly not a certificate of identification,
  reliable standard errors or likelihood accuracy.

# CoTiMA2 0.14.33

- Updates the ctsemjulia adapter for upstream revision 57a9758d: supports
  automatic sampling plus explicit nuts, saem, map and smc kernels, particles
  and coreset controls. NULL defaults delegate tuning and process allocation
  to ctsem, avoiding implicit NUTS controls that would invalidate SMC.
- Exposes backendArgs$intoverpop independently of the sampler in ctmaInit()
  and ctmaFit(). Rejects contradictory targets and SAEM placement without
  random effects early. Reports the native resolved route; joint sampling
  now has ctsemjulia_mode='joint', independent of its sampler.

# CoTiMA2 0.14.32

- Restricts `ctmaBiG()` SE validation to the DRIFT parameters used in its
  analysis. Zero SEs for fixed DIFFUSION/T0VAR cells no longer stop the fit.
  Invalid drift SEs are reported with study, matrix and parameter name;
  excluded studies are ignored. Empty estimation-problem lists no longer
  trigger a warning about failed primary-study fits.

# CoTiMA2 0.14.31

- Makes generated `ctModelText` easier to edit: the default model name is
  `m`, matrices use `matrix(..., nrow = ..., ncol = ...)`, and numeric values
  have no integer `L` suffix. The parameter table uses `data.frame()` and
  preserves post-construction constraints.

# CoTiMA2 0.14.30

- Fixes a post-fit `ctmaInit()` crash when a DIFFUSION or T0VAR interval
  diagnostic occurs without a preceding DRIFT diagnostic. Problems are now
  recorded independently with study and matrix-family labels; diagnostics
  cannot index an empty list or append a later study's problem to an earlier
  study. Missing parameter labels receive matrix-index fallbacks.

# CoTiMA2 0.14.29

- Adds Bayesian sampling to `backend = "ctsemjulia"` for `ctmaFit()` and
  `ctmaInit()` with compatible ctsem 4.x development versions. Native sampler
  settings are supplied through `backendArgs$sampleControl`; both `"nuts"`
  and `"saem"` are supported together with their common and sampler-specific
  controls. Random-effect NUTS fits retain the Laplace-marginal target, while
  SAEM uses `placement = "saem"` by default and samples the joint posterior.

# CoTiMA2 0.14.28

- Supports both generations of the experimental ctsem Julia/Laplace API:
  earlier versions exporting a separate `ctLaplaceCorrect()` function and
  ctsem 4.0.0 or newer versions integrating correction and refinement into
  `ctLaplaceCheck()`. The `ctsemjulia` availability check now reports the
  detected interface in `laplace_correction` instead of incorrectly rejecting
  the integrated API because `ctLaplaceCorrect()` is no longer exported.

- Suppresses ctsem's Stan-specific binary-approximation comparison while
  generating setup data for `backend = "julia"`. The Julia likelihood remains
  active, and all genuine setup warnings continue to be shown.

# CoTiMA2 0.14.27

- Preserves the complete fitted DRIFT and DIFFUSION matrices in `ctmaInit()`
  summaries when users fix selected matrix elements numerically. Fixed elements
  are no longer removed before assigning the full matrix-element names, and
  DIFFUSION results are no longer filtered using the free-parameter positions
  of DRIFT.

# CoTiMA2 0.14.26

- Fixes custom diffusion specifications in `ctmaInit()`. Free labels and fixed
  numerical values are now transferred correctly instead of referring to the
  nonexistent `diffParamsParams` object. Custom drift specifications now also
  preserve fixed nonzero values, and both specifications receive explicit
  dimension checks.
- Suppresses ctsem's internal “Maximum likelihood estimation requested” message
  during setup-only generation. The CoTiMA2 Julia backend now reports ML or
  Bayesian NUTS explicitly according to the estimator actually requested.
- Makes Bayesian NUTS robust to invalid intermediate covariance proposals in
  the CoTiMA2 Julia backend. Non-positive-definite, singular, and related
  numerical trial points are now rejected normally instead of terminating the
  complete chain with a `PosDefException` or nested `TaskFailedException`.
- Clarifies that the OTL printed by `summary()` is calculated from the mean
  estimated drift matrix. `ctmaOTL(..., CI = TRUE)` now repeats the numerical
  search across transformed drift draws and reports median OTLs and uncertainty
  intervals for ML/MAP or Bayesian fits. Lists and three-dimensional arrays of
  drift matrices are also supported.
- Preserves the requested binary manifest types when `ctmaInit()` rebuilds a
  model for process-based `randomIntercepts = "CINT"` or `"MANIFEST"`.
  Previously, the replacement model reverted to ctsem's all-continuous
  default even when `binaries` had been supplied.
- Corrects missing-row detection in `ctmaShapeRawData()`. Observed values or
  timestamps that sum to zero are no longer mistaken for entirely missing
  rows, including person-standardized data. The `withinTimeA` T0 check now
  tests actual completeness rather than numeric sums.
- Uses exact time-point suffixes in `ctmaShapeRawData()`, so matching `T1` no
  longer unintentionally includes `T10`, `T11`, and later waves. Long-to-wide
  wave counting now also uses the distinct numeric time suffixes.
- Restores the Stan path with the current ctsem development version. The
  removed historical `vb = FALSE` control is now discarded before calling
  `ctsem::ctFit()`, while a genuine `vb = TRUE` request produces an explicit
  CoTiMA2 error recommending ML/MAP optimization or Bayesian sampling.

# CoTiMA2 0.14.23

- `ctmaInit(..., fit = FALSE)` and `ctmaFit(..., fit = FALSE)` now return an
  executable `ctModelText` representation of the generated ctsem model. The
  returned summary explains that `eval(parse(text = result$ctModelText))` can
  be used to reconstruct the model, including CoTiMA2's final parameter table
  and post-construction constraints.

# CoTiMA2 0.14.22

- Adds `backend = "ctsemjulia"` support to `ctmaFitConfig()`. The function
  validates the matching unfitted model and native ctsem Julia/Laplace
  controls, retains them in `$recommended`, and reports setup time. It does
  not run CoTiMA2-specific subject-thread or ForwardDiff-chunk benchmarks,
  which are unavailable through the native ctsem interface.

# CoTiMA2 0.14.21

- Prefixes warnings generated while fitting or summarizing primary-study
  models in `ctmaInit()` with the study's position and, when available, its
  study number and source. Deferred warnings therefore remain attributable to
  the correct primary study.

# CoTiMA2 0.14.20

- Revises the Julia usage guide around the current high-level workflow:
  CoTiMA2's `backend = "julia"` can now fit both primary-study models in
  `ctmaInit()` and the pooled model in `ctmaFit()`.
- Adds a comparison of the `"stan"`, `"julia"`, and `"ctsemjulia"`
  backends, documents `ctmaInit()` thread allocation and uncertainty behavior,
  and adds workflow-specific `ctmaBackendDefaults()` examples.
- Retains `ctmaJuliaFit()` as the documented low-level validation and
  diagnostic interface, while removing obsolete `julia_file` placeholders
  from ordinary installed-package examples.
- Updates the workflow diagram and end-to-end template to use direct backend
  selection and `ctmaFitConfig()`.

# CoTiMA2 0.14.19

- Clarifies in every help page exposing backend controls that
  `ctmaBackendDefaults()` without arguments returns the Stan list rather than
  defaults for every backend.
- Makes the effective `ctsemjulia` ML default visible as
  `optimcontrol$uncertainty = "hessian"`, matching ctsem and providing standard
  errors and intervals for ordinary fits. Development diagnostics can still
  skip post-fit uncertainty explicitly with `estonly = TRUE`.

# CoTiMA2 0.14.18

- Normalizes optimizer diagnostics returned by the current ctsem 3.12.0
  `juliaFit` API under `$ctsemjulia_diagnostics`, while retaining the native
  `$optim` and `$estimate` components unchanged. The normalized view records
  convergence, gradient norm, iterations, and the number of converged inner
  Laplace modes.
- Marks native `subject_loglik` output as unsuitable for personwise LOO or
  WAIC when study-level random effects are integrated, or when the backend
  returns a constant placeholder vector. The original values remain available
  for reproducibility.
- Documents a reproducible native-Laplace reference-fit and inner-tolerance
  sensitivity workflow.

# CoTiMA2 0.14.17

- `ctmaFit(backend = "ctsemjulia")` now accepts pooled models containing
  subject-level random effects without requiring an additional study-level
  random effect. Thus, for example, `indVarying = "MANIFEST"` can use native
  Laplace integration while ordinary study-indicator effects remain fixed.
  Such subject-only models no longer construct an unused outer study-grouping
  level or print the corresponding study-membership message.

# CoTiMA2 0.14.16

- `ctmaInit(backend = "ctsemjulia")` now supports native ctsem Julia
  maximum-likelihood fits of individual primary studies. Models without
  person-level random effects use direct ML; ctsem's `intoverpop = "auto"`
  selects Laplace integration when the primary-study model contains such
  effects. The initial implementation requires `doPar = 1` to avoid nested
  R/Julia parallelism.
  Native `popsd` and raw-correlation summaries are converted into CoTiMA2's
  random-effect SD, covariance, and correlation output. In `ctmaInit()`,
  `indVarying` requests person-level variation and activates Laplace
  integration; `randomIntercepts` instead adds explicit latent intercept
  processes.

- Every CoTiMA function exposing a `backend` argument now detects an explicit
  `backend = "ctsemjulia"` request. If the compatible native ctsem Julia API is
  absent, the error provides an at-your-own-risk recipe that installs Charles
  Driver's `juliaFit` branch and `JuliaConnectoR` into an isolated R library,
  without replacing the user's production ctsem installation.

- `ctmaFit(backend = "ctsemjulia")` now uses native ctsem Julia ML directly
  for models without random effects and retains the Laplace path when subject-
  or study-level random effects are requested. The fitted object records the
  selected route in `$ctsemjulia_mode`.

# CoTiMA2 0.14.14

- Renamed package documentation sources from `CoTiMA_...` to `CoTiMA2_...`
  and synchronized the User Guide filename and title page with package version
  0.14.14. The tarball build script now stops when the User Guide name or
  embedded version does not match `DESCRIPTION`.
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
# CoTiMA2 0.14.49

- Match CPU control names exactly so the native optimcontrol$coreset logical flag is not mistaken for an optimcontrol$cores count. This fixes a pre-fit failure discovered in the isolated installation test.

# CoTiMA2 0.14.48

- Resource resolution honors available CPU restrictions and the requested core budget, including bundled Julia worker/BLAS products. Invalid or excessive explicit allocations fail before starting Julia.
- Owned subprocesses inherit single-thread BLAS/OpenMP limits, with the caller's environment restored on exit. This does not reconfigure an R BLAS library already loaded in the caller.
- Bundled Julia setup and fitting use a content-keyed writable user project instead of the installed package tree. Dry runs do not copy the project; dependency setup requires explicit consent.

# CoTiMA2 0.14.47

- Native ctsemjulia supports explicit quadrature refinement for optimized random-effects fits.
- Fit diagnostics distinguish checked and Laplace log likelihoods and expose source-reported optimizer resolution metadata.
- Quadrature retains integrated random-effects diagnostic safeguards; sampling and models without random effects reject this optimization-only request early.

# CoTiMA2 0.14.46

* Automatic Stan study-RE initialization follows the fitting seed (including a
  fit_args$seed override), separates chain streams and restores caller RNG state.
* ctmaOptimizeInit validates required inputs before retry execution and uses
  owned Unix fork workers without replacing the caller's foreach backend.
* ctmaGenData composes cluster cleanup with option restoration on exit.
* Existing output files now cause a clear noninteractive error unless
  silentOverwrite=TRUE; interactive rename choices recheck existing targets.

# CoTiMA2 0.14.45

* Reports cached Laplace approximation-accuracy status separately from optimizer
  convergence in ctmaFitDiagnostics(), native summary footers and table notes.
* Preserves unresolved or missing accuracy verdicts: an operational status of
  "checked" is not interpreted as an accuracy pass. No check or fit is rerun.
* Reports stored integration gap, maximum standardized shift and correction
  application when available, without altering parameter estimates or intervals.

# CoTiMA2 0.14.44

* Corrects mixed per-manifest mean specifications such as `c("mmV1", 0)`;
  explicit vectors no longer fail a scalar `if` condition or lose fixed entries.
* Adds regression coverage for vectors, numeric fixed values, scalar zero and auto.
* Documents the distinction between a binary measurement offset and a dynamic
  CINT, including why an initial mean fixed at zero does not encode a nonzero
  binary baseline prevalence. This is model guidance, not automatic refitting.

# CoTiMA2 0.14.43

* Reports stored bundled NUTS Rhat, bulk ESS and divergences in ctmaFitDiagnostics without refitting.

* Rejects unsupported bundled Julia MAP requests instead of silently fitting ML.
  Use `ctsemjulia` or `stan` with `optimize=TRUE, priors=TRUE` for MAP.
* Bayesian bundled Julia uncertainty uses retained posterior draws without
  requiring a likelihood Hessian at the initialization point.
* Filters removed inert Stan optimizer defaults only when the installed ctsem
  rejects them; explicit unsupported non-default requests produce an error.
* Routes pooled Stan fitting through the compatibility wrapper.
* Supports ctsem 4's renamed matrix registry in transformed study summaries.

# CoTiMA2 0.14.42

* Implements ctsem 4 z covariance conversion in the bundled Julia engine; previous off-diagonal code-2 fits need refitting. Unknown future codes fail before engine startup.
* Accepts native sampleControl$iter without changing its semantics.
* Separates posterior convergence from placement status, handles native pooled storage, and reports stored transformed diagnostics and SMC's non-applicable Rhat correctly.
* Corrects table convergence notes and reuses retained native matrix summaries.
* Clarifies the explicit Hessian default and native runtime selection.

# CoTiMA2 0.14.41

* `ctmaTables()` now accepts `fontSize = 10` (points) for titles, headers,
  cells and notes in Word, XML, TeX and PDF exports.

## Regular Stan table compatibility (2026-10-29)

* Reads joint population DRIFT arrays from regular optimized and sampled Stan fits.
* Reconstructs missing optimized-fit caches from retained raw draws through the
  fitted Stan program, with Kalman filtering disabled and no model refitting.
* Correctly labels native ML, MAP and Bayesian fits from their direct args fields.
* Retains native-summary precision for four-decimal tables, avoids partial cache
  name matching, and preserves a wrapper's CT representation during DT refresh.
* Validated against regular ctsem 3.11.1 and an existing 2000-draw CoTiMA benchmark.

# CoTiMA2 0.14.40

## Table presentation options (2026-10-29)

* layout selects portrait (long, default) or landscape (wide) pages consistently
  across Word, XML, TeX and PDF exports. Oversized tables continue in column panels.
* sortBy groups ctmaInit estimates by study (default) or by parameter.
* rowLayout keeps one parameter per row (default), combines DRIFT entries, or
  combines all selected parameters with grouped parameter/statistic headers.
* Numeric long-format values and safe staged output publication are preserved.

# CoTiMA2 0.14.39

## Saved-fit table summaries (2026-10-29)

* ctmaTables defaults to four display decimals.
* Missing statistics and DT uncertainty can be regenerated from saved native fits
  with summary(timeinterval=lag), with one consent for all primary studies.
* recomputeSummary="yes" permits batch use; no models are refitted and raw point
  estimates are never treated as uncertainty draws.

# CoTiMA2 0.14.38

## Table selectors and primary-study summaries (2026-10-28)

* Accepts bare, case-insensitive matrix families such as parameters=c(drift, diffusion).
  Quoted selectors and character-vector variables remain supported.
* Presents every primary study from a ctmaInit aggregate summary with original study
  labels, retained uncertainty and explicit aggregate fit-statistic notes; no refitting.
* Preserves unavailable statistics as unavailable and aligns study/parameter text
  columns in Word and TeX exports.

# CoTiMA2 0.14.37

## APA tables (2026-10-28)

* Adds ctmaTables() with checked summary agreement, named/matrix parameter selection,
  source-aware statistics, interactive selection and APA-style Word/XML/TeX/PDF exports.
* Adds joint-draw DT effects at selected fitted-time lags and coordinate-aligned model comparisons.
* Exports are staged; existing outputs require explicit overwrite. Unsupported summaries
  and unavailable statistics fail explicitly rather than fabricating transformed values.
* Documents optional Word/PDF dependencies and the limits of legacy summary provenance.

# CoTiMA2 0.14.36

## Power analysis refactoring (2026-10-28)

* Separates validation, invariant fitting/loading, matrix extraction, implied moments,
  regression effects, power probabilities, integer sample-size search and result assembly.
* Verifies unrestricted Gaussian regression eligibility for direct covariance calculations.
  Supports explicit restricted Gaussian regression syntax through a lavaan fallback, with
  console and stored summary notes identifying the calculation route.
* Corrects negative-effect p-values, two-wave and higher-dimensional calculations, lag
  bounds, reduced hypotheses, root boundaries and tied optimal-lag selection.
* Fast default power includes both rejection tails. useMBESS=TRUE retains the legacy
  probability convention. Missing post-hoc powers stay NA instead of becoming zero.
* Restores original time units from invariant-fit scaling metadata. Reduced-fit caches
  now validate covariance/model metadata and are reused across requested power levels.

# CoTiMA2 0.14.35

# CoTiMA2 0.14.35

## Plotting corrections (2026-10-28)

* Correct forest pooled confidence intervals and funnel confidence limits.
* Correct time-range indexing, axis positions, power-only plots and fractional-lag selection.
* Handle multiple fits, legacy example objects and per-fit moderator/lag mappings.
* Export uniquely named PNG files in the requested directory, including headless sessions.

# CoTiMA2 0.14.34

- Separates optimizer, identification, transform-saturation and posterior
  diagnostics in Julia fit summaries. Adds ctmaFitDiagnostics() for stored
  fits without refitting or starting Julia. Unknown diagnostics remain unknown.
- New native fits cache population matrices for a labelled heuristic flag
  when binary linear-predictor SD exceeds 10; no parameter is fixed or changed.
  Optimizer convergence is explicitly not a certificate of identification,
  reliable standard errors or likelihood accuracy.

# CoTiMA2 0.14.33

- Updates the ctsemjulia adapter for upstream revision 57a9758d: supports
  automatic sampling plus explicit nuts, saem, map and smc kernels, particles
  and coreset controls. NULL defaults delegate tuning and process allocation
  to ctsem, avoiding implicit NUTS controls that would invalidate SMC.
- Exposes backendArgs$intoverpop independently of the sampler in ctmaInit()
  and ctmaFit(). Rejects contradictory targets and SAEM placement without
  random effects early. Reports the native resolved route; joint sampling
  now has ctsemjulia_mode='joint', independent of its sampler.

# CoTiMA2 0.14.32

- Restricts `ctmaBiG()` SE validation to the DRIFT parameters used in its
  analysis. Zero SEs for fixed DIFFUSION/T0VAR cells no longer stop the fit.
  Invalid drift SEs are reported with study, matrix and parameter name;
  excluded studies are ignored. Empty estimation-problem lists no longer
  trigger a warning about failed primary-study fits.

# CoTiMA2 0.14.31

- Makes generated `ctModelText` easier to edit: the default model name is
  `m`, matrices use `matrix(..., nrow = ..., ncol = ...)`, and numeric values
  have no integer `L` suffix. The parameter table uses `data.frame()` and
  preserves post-construction constraints.

# CoTiMA2 0.14.30

- Fixes a post-fit `ctmaInit()` crash when a DIFFUSION or T0VAR interval
  diagnostic occurs without a preceding DRIFT diagnostic. Problems are now
  recorded independently with study and matrix-family labels; diagnostics
  cannot index an empty list or append a later study's problem to an earlier
  study. Missing parameter labels receive matrix-index fallbacks.

# CoTiMA2 0.14.29

- Adds Bayesian sampling to `backend = "ctsemjulia"` for `ctmaFit()` and
  `ctmaInit()` with compatible ctsem 4.x development versions. Native sampler
  settings are supplied through `backendArgs$sampleControl`; both `"nuts"`
  and `"saem"` are supported together with their common and sampler-specific
  controls. Random-effect NUTS fits retain the Laplace-marginal target, while
  SAEM uses `placement = "saem"` by default and samples the joint posterior.

# CoTiMA2 0.14.28

- Supports both generations of the experimental ctsem Julia/Laplace API:
  earlier versions exporting a separate `ctLaplaceCorrect()` function and
  ctsem 4.0.0 or newer versions integrating correction and refinement into
  `ctLaplaceCheck()`. The `ctsemjulia` availability check now reports the
  detected interface in `laplace_correction` instead of incorrectly rejecting
  the integrated API because `ctLaplaceCorrect()` is no longer exported.

- Suppresses ctsem's Stan-specific binary-approximation comparison while
  generating setup data for `backend = "julia"`. The Julia likelihood remains
  active, and all genuine setup warnings continue to be shown.

# CoTiMA2 0.14.27

- Preserves the complete fitted DRIFT and DIFFUSION matrices in `ctmaInit()`
  summaries when users fix selected matrix elements numerically. Fixed elements
  are no longer removed before assigning the full matrix-element names, and
  DIFFUSION results are no longer filtered using the free-parameter positions
  of DRIFT.

# CoTiMA2 0.14.26

- Fixes custom diffusion specifications in `ctmaInit()`. Free labels and fixed
  numerical values are now transferred correctly instead of referring to the
  nonexistent `diffParamsParams` object. Custom drift specifications now also
  preserve fixed nonzero values, and both specifications receive explicit
  dimension checks.
- Suppresses ctsem's internal “Maximum likelihood estimation requested” message
  during setup-only generation. The CoTiMA2 Julia backend now reports ML or
  Bayesian NUTS explicitly according to the estimator actually requested.
- Makes Bayesian NUTS robust to invalid intermediate covariance proposals in
  the CoTiMA2 Julia backend. Non-positive-definite, singular, and related
  numerical trial points are now rejected normally instead of terminating the
  complete chain with a `PosDefException` or nested `TaskFailedException`.
- Clarifies that the OTL printed by `summary()` is calculated from the mean
  estimated drift matrix. `ctmaOTL(..., CI = TRUE)` now repeats the numerical
  search across transformed drift draws and reports median OTLs and uncertainty
  intervals for ML/MAP or Bayesian fits. Lists and three-dimensional arrays of
  drift matrices are also supported.
- Preserves the requested binary manifest types when `ctmaInit()` rebuilds a
  model for process-based `randomIntercepts = "CINT"` or `"MANIFEST"`.
  Previously, the replacement model reverted to ctsem's all-continuous
  default even when `binaries` had been supplied.
- Corrects missing-row detection in `ctmaShapeRawData()`. Observed values or
  timestamps that sum to zero are no longer mistaken for entirely missing
  rows, including person-standardized data. The `withinTimeA` T0 check now
  tests actual completeness rather than numeric sums.
- Uses exact time-point suffixes in `ctmaShapeRawData()`, so matching `T1` no
  longer unintentionally includes `T10`, `T11`, and later waves. Long-to-wide
  wave counting now also uses the distinct numeric time suffixes.
- Restores the Stan path with the current ctsem development version. The
  removed historical `vb = FALSE` control is now discarded before calling
  `ctsem::ctFit()`, while a genuine `vb = TRUE` request produces an explicit
  CoTiMA2 error recommending ML/MAP optimization or Bayesian sampling.

# CoTiMA2 0.14.23

- `ctmaInit(..., fit = FALSE)` and `ctmaFit(..., fit = FALSE)` now return an
  executable `ctModelText` representation of the generated ctsem model. The
  returned summary explains that `eval(parse(text = result$ctModelText))` can
  be used to reconstruct the model, including CoTiMA2's final parameter table
  and post-construction constraints.

# CoTiMA2 0.14.22

- Adds `backend = "ctsemjulia"` support to `ctmaFitConfig()`. The function
  validates the matching unfitted model and native ctsem Julia/Laplace
  controls, retains them in `$recommended`, and reports setup time. It does
  not run CoTiMA2-specific subject-thread or ForwardDiff-chunk benchmarks,
  which are unavailable through the native ctsem interface.

# CoTiMA2 0.14.21

- Prefixes warnings generated while fitting or summarizing primary-study
  models in `ctmaInit()` with the study's position and, when available, its
  study number and source. Deferred warnings therefore remain attributable to
  the correct primary study.

# CoTiMA2 0.14.20

- Revises the Julia usage guide around the current high-level workflow:
  CoTiMA2's `backend = "julia"` can now fit both primary-study models in
  `ctmaInit()` and the pooled model in `ctmaFit()`.
- Adds a comparison of the `"stan"`, `"julia"`, and `"ctsemjulia"`
  backends, documents `ctmaInit()` thread allocation and uncertainty behavior,
  and adds workflow-specific `ctmaBackendDefaults()` examples.
- Retains `ctmaJuliaFit()` as the documented low-level validation and
  diagnostic interface, while removing obsolete `julia_file` placeholders
  from ordinary installed-package examples.
- Updates the workflow diagram and end-to-end template to use direct backend
  selection and `ctmaFitConfig()`.

# CoTiMA2 0.14.19

- Clarifies in every help page exposing backend controls that
  `ctmaBackendDefaults()` without arguments returns the Stan list rather than
  defaults for every backend.
- Makes the effective `ctsemjulia` ML default visible as
  `optimcontrol$uncertainty = "hessian"`, matching ctsem and providing standard
  errors and intervals for ordinary fits. Development diagnostics can still
  skip post-fit uncertainty explicitly with `estonly = TRUE`.

# CoTiMA2 0.14.18

- Normalizes optimizer diagnostics returned by the current ctsem 3.12.0
  `juliaFit` API under `$ctsemjulia_diagnostics`, while retaining the native
  `$optim` and `$estimate` components unchanged. The normalized view records
  convergence, gradient norm, iterations, and the number of converged inner
  Laplace modes.
- Marks native `subject_loglik` output as unsuitable for personwise LOO or
  WAIC when study-level random effects are integrated, or when the backend
  returns a constant placeholder vector. The original values remain available
  for reproducibility.
- Documents a reproducible native-Laplace reference-fit and inner-tolerance
  sensitivity workflow.

# CoTiMA2 0.14.17

- `ctmaFit(backend = "ctsemjulia")` now accepts pooled models containing
  subject-level random effects without requiring an additional study-level
  random effect. Thus, for example, `indVarying = "MANIFEST"` can use native
  Laplace integration while ordinary study-indicator effects remain fixed.
  Such subject-only models no longer construct an unused outer study-grouping
  level or print the corresponding study-membership message.

# CoTiMA2 0.14.16

- `ctmaInit(backend = "ctsemjulia")` now supports native ctsem Julia
  maximum-likelihood fits of individual primary studies. Models without
  person-level random effects use direct ML; ctsem's `intoverpop = "auto"`
  selects Laplace integration when the primary-study model contains such
  effects. The initial implementation requires `doPar = 1` to avoid nested
  R/Julia parallelism.
  Native `popsd` and raw-correlation summaries are converted into CoTiMA2's
  random-effect SD, covariance, and correlation output. In `ctmaInit()`,
  `indVarying` requests person-level variation and activates Laplace
  integration; `randomIntercepts` instead adds explicit latent intercept
  processes.

- Every CoTiMA function exposing a `backend` argument now detects an explicit
  `backend = "ctsemjulia"` request. If the compatible native ctsem Julia API is
  absent, the error provides an at-your-own-risk recipe that installs Charles
  Driver's `juliaFit` branch and `JuliaConnectoR` into an isolated R library,
  without replacing the user's production ctsem installation.

- `ctmaFit(backend = "ctsemjulia")` now uses native ctsem Julia ML directly
  for models without random effects and retains the Laplace path when subject-
  or study-level random effects are requested. The fitted object records the
  selected route in `$ctsemjulia_mode`.

# CoTiMA2 0.14.14

- Renamed package documentation sources from `CoTiMA_...` to `CoTiMA2_...`
  and synchronized the User Guide filename and title page with package version
  0.14.14. The tarball build script now stops when the User Guide name or
  embedded version does not match `DESCRIPTION`.
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

### Primary-study identification reporting (unreleased)
- Preserve explicit native ctsem no-width flags as NA in parameter-matrix and combined primary-study uncertainty tables, including stored-fit summary refresh without refitting. Keep point estimates unchanged.
- Add study-labelled observed two-wave random-intercept design warnings and combined identification, boundary and unstable-interval notes. Wave count is not a universal identification test.
- Interval-problem labels now use original study identifiers rather than list positions.

### Observed study-design partitions (unreleased)
- ctmaInit returns studyDesign and ctmaPrep-compatible studyLists$twoWaves / studyLists$threeOrMoreWaves (NULL for empty groups), preserving study numbers, order and metadata. No studies are excluded.
- Study-labelled two-wave warnings distinguish separately estimated random-intercept covariance from the shared covariance requested by indVarying in pooled ctmaFit models. Explicitly pooled phantom covariance is also exempt from the separate-study warning.

### Native pooled study-specific determined diffusion (unreleased)
- Support shared drift with fully study-specific stationary within-person covariance in ctmaFit(backend="ctsemjulia", stationary=TRUE, determinedDiffusions=TRUE). Saturated dummy/weighted study contrasts act on the covariance chart; derived diffusion varies by study and the free parameter count is retained.
- Reject general process moderators, partially varying covariance coordinates and study-varying drift rather than silently changing invariance constraints. Existing ML/Gaussian/random-effect safeguards remain. Other backend support for common-process determined models is unchanged.

* Added the non-exported research adapter ctmaLohmannFit for literal FE/URE/ERE/CRE
  replication from a locally supplied Lohmann et al. supplement. Separate internal
  help and inst/research documentation explain input ordering, explicit sampling
  covariances, source fingerprints and preserved mathematical/code limitations.
  No existing ctmaFit behavior changes and no posterior validation is claimed.

## Development after 0.14.59: User Guide selection compatibility

* Retain the historical empty NA-study closing row produced by ctmaPrep.
  ctmaFit now recognizes that row at its primaryStudyList selection boundary,
  only when the recorded count and study payloads confirm it is empty. Real
  missing IDs, populated closing rows and misalignment remain errors. No change
  to study counts, contrast counts or the existing legacy consumers is intended.

* Undefined transformed random-effect correlation draws no longer abort ctmaInit
  during quantile extraction. Affected quantiles remain NA with a study-labelled
  warning; finite-coordinate quantiles and model estimation are unchanged.
