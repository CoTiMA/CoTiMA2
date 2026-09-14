# Reproducible report: Laplace inner-mode tolerance blocks outer optimisation

Dear Charlie,

While validating a pooled CoTiMA2 model with person-level random effects only,
we found that the native Julia/Laplace outer optimiser can stop at zero
iterations because none of the subject-level modes satisfies the current inner
convergence criterion. The model has 571 subjects, one random-effects level,
and four person-level random effects. No outer study grouping is present.

The installed and current `juliaFit` revision is
`2140ceb97f0ec4919d6f129b0edda83c2298007e` (`ctsem` 3.12.0).

## Results

1. With the unmodified defaults (`inner_maxiter = 50`,
   `inner_tol = 1e-10`), all 571 inner solves reached 50 iterations. The largest
   remaining inner gradient was approximately `1.35e-5`. Zero objective steps
   were accepted, the outer optimiser reported zero iterations, and the fit was
   marked as stalled.
2. Disabling ctsem's data-derived starting values produced the same result:
   0/571 inner modes converged and the outer optimiser stopped at zero
   iterations. Thus the data-start heuristic was not the cause.
3. In an isolated copy of the Julia engine, increasing
   `inner_maxiter` to 200 while using `inner_tol = 1e-6` improved inner
   convergence to 281/571, but the remaining modes stalled near `1.36e-5` and
   the outer optimiser still could not start.
4. With `inner_maxiter = 200` and `inner_tol = 2e-5`, all 571 inner modes
   converged. The outer optimiser then completed its requested 20 diagnostic
   iterations, improving the log posterior from approximately -5378 to -3535.
   It was not expected to converge within this deliberately short outer run,
   but it was no longer stalled.

These comparisons used `optimcontrol$estonly = TRUE`, so no Hessian or
posterior approximation was computed. They also confirmed that subject-only
models do not require an outer study-membership level.

## API requests

Could `inner_maxiter` and `inner_tol` be exposed through `optimcontrol`, rather
than remaining fixed Julia-constructor defaults? A scale-aware convergence
criterion or a somewhat more permissive default may also be appropriate, but
the `2e-5` value above should be treated as diagnostic evidence rather than a
proposed universal default until likelihood and gradient accuracy have been
checked.

It would also be helpful if `.ctFitJuliaBackendImpl()` skipped uncertainty and
exact-Hessian computation whenever the outer optimiser is stalled or the
inner modes have not converged. In the original full call, the failed
zero-iteration fit proceeded to a lengthy exact Hessian calculation, where all
parameter columns were `NaN` because the inner solve failed at the Hessian
steps.

Best,

Christian
