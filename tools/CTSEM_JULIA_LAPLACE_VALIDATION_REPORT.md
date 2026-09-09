# Validation report: ctsem native Julia/Laplace multilevel backend

Dear Charlie,

Thank you for sending the multilevel Laplace vignette. I have now reproduced
the example and tested the development API through an experimental internal
CoTiMA2 adapter. Overall, the results are very encouraging. The native
Julia/Laplace approach appears to be a strong candidate for hierarchical
random-effects estimation in CoTiMA2, particularly as an alternative to the
currently very expensive Stan/NUTS route.

## Test environment

- macOS on Apple silicon (`aarch64-apple-darwin23`)
- R 4.6.1
- Julia 1.12.7 installed through `juliaup`
- `ctsem` 3.12.0, branch `juliaFit`
- tested commit: `92c2d17f8648afdabc23a144e3a914de993e70af`
- `JuliaConnectoR` 1.1.6

The development version of `ctsem` and `JuliaConnectoR` were installed in a
separate R library. The existing production installation of `ctsem` was not
replaced.

## 1. Reproduction of the vignette example

I reproduced the 30-subject, five-study example from
`laplace-multilevel.html`, with subject-level random effects on `T0MEANS` and
`CINT` and a study-level random effect on `CINT`.

- optimization elapsed time: approximately 31 seconds
- `ctLaplaceCheck()` elapsed time: approximately 10 seconds
- Laplace log marginal: -287.32628
- quadrature log marginal: -287.32628
- gap: numerically zero
- largest reported first-order corrections: approximately `1e-12` standard
  errors

This confirms the expected exactness for the identity-transformed,
linear-Gaussian random effects in this model.

The Julia session had four threads available. The chunk tuner sometimes chose
one or two chunks for these small problems, which was sensible and should not
be confused with a failure to start the requested Julia thread pool.

## 2. CoTiMA2 adapter checks with linear random effects

I implemented an internal, currently unexported adapter that:

1. checks for the required development API;
2. adds a strict subject-within-study hierarchy to a CoTiMA/ctsem model;
3. adds the appropriate `indvarying_<level>` and `sdscale_<level>` columns;
4. marks selected subject- and study-level random effects;
5. removes redundant fixed study-dummy effects for parameters selected as
   study random effects; and
6. calls `ctFit(..., backend = "julia", intoverpop = "laplace")`.

The first genuine CoTiMA2 test used three simulated studies with identical true
effects, 30 subjects per study, and subject- and study-level CINT effects.

- optimization: 42.3 seconds
- `ctLaplaceCheck()`: 33.0 seconds
- marginal-likelihood gap: `-2.27e-13`
- several population SDs collapsed to zero
- five information-matrix directions were weakly identified

The collapsed scales and weak directions are consistent with attempting to
estimate heterogeneity from only three studies whose true effects are equal.
They do not indicate a Laplace approximation problem.

A second test used nine groups with deliberately imposed, oppositely directed
shifts in two outcomes.

- 270 subjects in total
- optimization: 81.0 seconds
- `ctLaplaceCheck()`: 102.2 seconds
- convergence after 299 iterations
- gradient norm: `7.67e-6`
- marginal-likelihood gap: numerically zero
- estimated study-level CINT SDs: 0.135 and 0.147
- estimated study-level CINT correlation: -0.727, in the expected direction

No study-level scale collapsed in this heterogeneous condition.

## 3. Nonlinear study-level DRIFT validation

I next simulated a univariate OU process with subject-level variation in CINT
and study-level variation in DRIFT. This is a genuinely approximate case,
because the DRIFT random effect enters through ctsem's nonlinear transform.

### Initial nonlinear test

- 12 studies
- 30 subjects per study
- 10 observations per subject
- optimization: 44.8 seconds
- `ctLaplaceCheck()`: 10.8 seconds
- converged: yes
- gradient norm: `5.75e-9`
- estimated study-level DRIFT SD: 0.144
- Laplace/quadrature gap: 0.2255 log units in total, or 0.000626 per subject
- largest relevant correction: 0.042 standard errors

A subject-level CINT scale collapsed to zero in this relatively short design.
The resulting extremely large `delta_se` for that parameter belonged to a
dropped, unidentified boundary direction and was not interpretable as an
approximation error. Reading `delta_se` together with `dropped_directions`, as
the vignette recommends, was essential here.

### Stronger parameter-recovery test

For a cleaner recovery test, I simulated DRIFT on the exact raw scale used by
ctsem. The raw study effects had a realized SD of exactly 0.30 and were then
mapped through the model's DRIFT transform.

- 16 studies
- 35 subjects per study
- 20 observations per subject
- true raw study-level DRIFT SD: 0.300
- optimization: 61.7 seconds
- `ctLaplaceCheck()`: 17.5 seconds
- converged: yes
- gradient norm: `9.88e-7`
- estimated raw study-level DRIFT SD: 0.376
- 95% interval: 0.228 to 0.574
- relative point-recovery error: 25.3%
- Laplace/quadrature gap: 1.525 log units in total, or 0.00272 per subject

The true study-level SD was inside the reported interval. With only 16 studies,
the point deviation does not by itself imply systematic bias.

The largest approximation corrections were:

- measurement-error parameter: 0.201 SE
- diffusion parameter: -0.158 SE
- study-level DRIFT SD: 0.042 SE
- population DRIFT: 0.020 SE
- subject-level CINT SD: -0.031 SE

Thus, all corrections were well below one standard error. In particular, the
random-effects quantity of primary interest was almost unaffected by the
Laplace approximation.

`ctLaplaceCorrect()` also completed successfully:

- correction elapsed time: 32.7 seconds
- DRIFT SD before correction: 0.376
- DRIFT SD after correction: 0.378
- subject-level CINT SD before/after correction: approximately 0.392/0.391

The correction therefore made no substantive difference in this test.

## 4. Issue: study-only random effects fail the early Laplace gate

**Update, 9 September 2026:** Charlie resolved this issue in commit
`32e8d86`. I subsequently repeated a study-level-only fit with zero
subject-level random effects and one study-level random effect. The model
passed setup, produced a finite log likelihood of -58.952, and converged. The
CoTiMA2 adapter now detects the corrected outer-level API and no longer
requires an artificial subject-level random effect.

The development version originally tested rejected a model that had a random
effect only at an outer grouping level and no subject-level random effect. For
example, the following intended structure failed before commit `32e8d86`:

```r
m <- ctModel(
  type = "ct",
  id = c("id", "study"),
  manifestNames = "Y1",
  latentNames = "eta1",
  LAMBDA = matrix(1),
  CINT = matrix("cint_eta1"),
  MANIFESTMEANS = matrix(0)
)

m$pars$indvarying <- FALSE
m$pars$indvarying_study <- FALSE
m$pars$indvarying_study[m$pars$matrix == "CINT"] <- TRUE

fit <- ctFit(
  dat, m,
  backend = "julia",
  intoverpop = "laplace",
  cores = 4
)
```

The error is:

```text
intoverpop='laplace' was requested but no free parameters are marked
indvarying, so there is nothing to integrate over.
```

The earlier check in `ctFit.R` examined only
`ctm$pars$indvarying`:

```r
if (!any(ctm$pars$indvarying[is.na(ctm$pars$value)])) ...
```

It appears that this gate should consider every hierarchy-level varying column,
including columns matching `^indvarying_`. The downstream representation may
already support such a model; I have not bypassed the check to assume that.

Before the upstream correction, the CoTiMA2 adapter refused a study-only
specification rather than silently introducing person heterogeneity merely to
pass the gate. That temporary restriction has now been removed.

This issue is now resolved and retained here to document what was tested.

## 5. Juliaup shim and `JULIA_BINDIR`

With:

```text
which julia
/Users/cdormann/.juliaup/bin/julia
```

the initial `JuliaConnectoR` startup treated the juliaup shim directory as the
runtime directory and attempted to load:

```text
/Users/cdormann/.juliaup/bin/../lib/julia/sys.dylib
```

That file does not exist. Pointing `JULIA_BINDIR` to the real versioned Julia
binary directory solved the problem:

```r
Sys.setenv(
  JULIA_BINDIR = paste0(
    "/Users/cdormann/.julia/juliaup/",
    "julia-1.12.7+0.aarch64.apple.darwin14/",
    "Julia-1.12.app/Contents/Resources/julia/bin"
  )
)
```

This behavior may originate in `JuliaConnectoR` rather than ctsem. However,
ctsem's setup helper may be able to make juliaup installations more robust by
resolving the shim through a command such as:

```sh
julia -e 'println(Sys.BINDIR)'
```

before passing the directory to the connector.

## 6. Overall assessment and next integration step

Across the exact CINT tests and the nonlinear DRIFT tests:

- optimization was stable and fast;
- exact cases agreed with quadrature to numerical precision;
- nonlinear DRIFT corrections remained small relative to uncertainty;
- `ctLaplaceCheck()` and `ctLaplaceCorrect()` behaved as documented; and
- the hierarchical representation is much more promising computationally than
  the present CoTiMA2 Stan/NUTS study-random-effects route.

The CoTiMA2 adapter remains internal and feature-gated. We have not exposed a
public `backend = "ctsemjulia"`, because that should wait until the API is
stable and the outer-level-only behavior is clarified.

I would be grateful for your comments on:

1. whether study-only random effects should work without any subject-level
   random effect;
2. whether the API used here is sufficiently stable for an experimental
   CoTiMA2 backend;
3. whether you recommend additional validation cases, especially for multiple
   nonlinear random effects or deeper nesting; and
4. whether resolving a juliaup shim belongs in ctsem or should be reported
   directly to `JuliaConnectoR`.

Best wishes,

Christian
