# Julia performance notes

This document records practical findings about the experimental Julia backend
so that recommendations from development and testing remain available outside
individual Codex tasks.

## Thread controls

- `julia_threads` is the size of Julia's thread pool. It is capacity, not a
  guarantee that every thread is active during every fitting stage.
- `subject_threads` limits the number of concurrent subject-processing tasks.
  For subject-wise likelihood and gradient work, this value can be the actual
  concurrency bottleneck even when `julia_threads` is larger.
- `blas_threads` controls BLAS-internal parallelism. A value of `1L` normally
  avoids oversubscription when Julia already parallelizes work at subject level.

For example, a report of

```text
9 Julia pool thread(s), 1 concurrent subject task(s), and 1 BLAS thread(s)
```

means that Julia has nine pool threads available but CoTiMA processes only one
subject task at a time. Low total CPU utilization is therefore expected.

## macOS Activity Monitor

JuliaCall embeds Julia in the R process; `julia_threads = 9L` does not create
nine separate Julia processes. In Activity Monitor, inspect the R, RStudio, or
`rsession` process. Roughly 100% CPU represents one fully occupied core, while
900% would represent about nine fully occupied cores. Utilization can be much
lower when the current fitting stage is serial or `subject_threads` is small.

The configured Julia thread count can be checked in the same R session with:

```r
JuliaCall::julia_eval("Threads.nthreads()")
```

## Consistency safeguards

Parallel likelihoods and parallel AD gradients are enabled independently only
after exact serial-consistency checks. If either check fails, CoTiMA2 uses the
deterministic serial-safe path for that operation. Increasing thread settings
then does not necessarily increase speed. Hessians use the serial-safe path.

## Current observation to investigate

For a large random-intercept model with 157 parameters, `ctmaFitConfig()`
recommended approximately:

```r
julia_threads   = 9L
subject_threads = 1L
blas_threads    = 1L
```

The resulting fit showed low CPU utilization and required about ten hours to
approach 5,000 iterations. We need to inspect the complete configuration
benchmark results and determine whether the short benchmark systematically
favours `subject_threads = 1L` for large RI models even when larger values would
reduce full-fit time.

Do not assume that the fastest short likelihood benchmark also yields the
fastest or most stable full optimization. Candidate values such as 1, 2, 4,
and 8 subject tasks should be compared while checking serial/parallel agreement.
