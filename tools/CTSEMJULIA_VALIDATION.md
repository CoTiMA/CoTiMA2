# CoTiMA2 native ctsem/Julia validation suite

This suite tests whether `backend = "ctsemjulia"` recovers known
hierarchical parameters and whether the Laplace approximation agrees with
adaptive quadrature closely enough for further development.

It is deliberately separate from ordinary package tests. A native Julia fit
requires Julia, the Julia packages used by ctsem, and a ctsem version that
provides the experimental multilevel API. CRAN checks must not launch such a
fit.

## What is tested

The simulations use a one-variable continuous-time Ornstein--Uhlenbeck model.
Between-study standard deviations are generated on ctsem's unconstrained raw
scale, which is also the scale of the Gaussian random effects estimated by
ctsem.

- `cint`: subject- and study-level CINT heterogeneity;
- `drift`: subject-level CINT and study-level DRIFT heterogeneity;
- `diffusion`: study-level DIFFUSION heterogeneity;
- `t0var`: study-level T0VAR heterogeneity.

The quick mode runs the first two cases. Full mode runs all four with more
studies, people, and observations. Each fit is followed by
`ctsem::ctLaplaceCheck(nodes = 5)`.

## Running it on macOS or Ubuntu

Open Terminal, change to the CoTiMA2 project directory, and execute:

```sh
Rscript tools/ctsemjulia-validation-suite.R --mode=quick --cores=4
```

For the full validation:

```sh
Rscript tools/ctsemjulia-validation-suite.R --mode=full --cores=4
```

To repeat only selected scenarios, provide a comma-separated list:

```sh
Rscript tools/ctsemjulia-validation-suite.R \
  --mode=full --scenarios=diffusion,t0var --cores=4
```

If the required ctsem development version is installed in a separate R
library, specify it explicitly:

```sh
Rscript tools/ctsemjulia-validation-suite.R \
  --mode=full \
  --cores=4 \
  --development-library=/path/to/the/R/library
```

If Julia cannot be found automatically, add its binary directory:

```sh
Rscript tools/ctsemjulia-validation-suite.R \
  --mode=quick \
  --julia-bindir=/path/to/julia/bin
```

Options use the portable `--name=value` form. The same commands work on macOS
and Ubuntu; only paths may differ.

## Results

By default, results are written to
`~/CoTiMA_Benchmark/ctsemjulia-validation/`:

- `REPORT.md`: readable automatic report;
- `summary.csv`: compact comparison table;
- `validation-results.rds`: provenance, metrics, and diagnostics;
- one `.rds` file per scenario.

Use `--output=/another/directory` to change the location and
`--retain-fits=true` to retain the potentially large fitted objects.

`PASS` and `WARN` are engineering screening decisions. They are not evidence
of unbiased estimation or nominal confidence-interval coverage. Promotion of
the backend from experimental status requires replicated simulations across
seeds, model dimensions, amounts of heterogeneity, and difficult boundary
conditions.

A fixed-effects Stan or CoTiMA2-Julia fit is intentionally not reported as a
reference estimate for the hierarchical Laplace cases: it optimizes a
different likelihood and targets a different estimand. Cross-backend timing
and numerical comparisons become valid only when the backends express the
same random-effects model.

## Resumable Monte Carlo validation

The single-run checks above can be replicated with independent simulation
seeds by the portable controller:

```sh
Rscript tools/ctsemjulia-monte-carlo.R \
  --mode=full \
  --replications=10 \
  --workers=4 \
  --cores-per-fit=4
```

Every replication runs in its own R/Julia process, while its four scenarios
run sequentially in that process to amortize Julia startup. Do not request
more than the available physical CPU cores; leaving at least one core free is
enforced. On a laptop, start with one worker. On a many-core workstation,
increase the number of workers only after a small timing pilot.

Results are checkpointed under `replications/rep-NNNN/` immediately after
each scenario. Repeating the same command reuses only completed files whose
scenario, seed, and design dimensions agree; failed or interrupted scenarios
are run again. `--replication-start=` permits a separate, non-overlapping
batch. `--aggregate-only=true` rebuilds the reports without fitting.

The controller creates `replication-summary.csv`, `monte-carlo-summary.csv`,
`MONTE_CARLO_REPORT.md`, `monte-carlo-results.rds`, per-replication logs, and
a run manifest. Aggregate diagnostics include convergence, bias, relative
bias, RMSE, absolute error, Laplace/quadrature gaps, timings, and 95-percent
interval coverage whenever ctsem supplies finite limits for the random-effect
standard deviation.
