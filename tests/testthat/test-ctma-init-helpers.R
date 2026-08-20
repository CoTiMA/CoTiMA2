test_that("Stan control overrides are merged by name and recursively", {
  defaults <- list(
    iter = 1000,
    verbose = 0,
    optimcontrol = list(is = FALSE, finishsamples = 1000)
  )
  overrides <- list(
    verbose = 2,
    iter = 500,
    optimcontrol = list(is = NULL)
  )

  result <- CoTiMA:::ctmaMergeStanctArgs(defaults, overrides)

  expect_identical(result$iter, 500)
  expect_identical(result$verbose, 2)
  expect_true("is" %in% names(result$optimcontrol))
  expect_null(result$optimcontrol$is)
  expect_identical(result$optimcontrol$finishsamples, 1000)
})

test_that("default importance-sampling controls are NULL", {
  expect_null(CoTiMA::CoTiMAStanctArgs$optimcontrol$is)
  expect_null(CoTiMA::CoTiMAStanctArgs$optimcontrol$isitersize)
  expect_null(CoTiMA::CoTiMAStanctArgs$optimcontrol$isESS)
})

test_that("Stan control overrides reject unnamed and unknown values", {
  defaults <- list(iter = 1000)

  expect_error(
    CoTiMA:::ctmaMergeStanctArgs(defaults, list(500)),
    "fully named list"
  )
  expect_error(
    CoTiMA:::ctmaMergeStanctArgs(defaults, list(unknown = 1)),
    "Unknown CoTiMAStanctArgs: unknown"
  )
})

test_that("ctmaInit resolves legacy Stan controls recursively by name", {
  result <- CoTiMA:::.ctmaResolveInitStanArgs(
    overrides = list(
      verbose = 2,
      iter = 750,
      optimcontrol = list(finishsamples = 250)
    )
  )

  expect_identical(result$verbose, 2)
  expect_identical(result$iter, 750)
  expect_identical(result$optimcontrol$finishsamples, 250)
  expect_identical(
    result$optimcontrol$carefulfit,
    CoTiMAStanctArgs$optimcontrol$carefulfit
  )
})

test_that("explicit ctmaInit arguments override the legacy Stan list", {
  result <- CoTiMA:::.ctmaResolveInitStanArgs(
    overrides = list(
      scaleTI = FALSE,
      iter = 400,
      verbose = 1,
      optimcontrol = list(finishsamples = 100)
    ),
    scaleTI = TRUE,
    iter = 900,
    verbose = 3,
    finishsamples = 500
  )

  expect_true(result$scaleTI)
  expect_identical(result$iter, 900)
  expect_identical(result$verbose, 3)
  expect_identical(result$optimcontrol$finishsamples, 500)
})

test_that("ctmaInit backend merge detects misordered and unknown controls", {
  result <- CoTiMA:::.ctmaResolveInitStanArgs(list(
    verbose = 4,
    iter = 321
  ))
  expect_identical(result$verbose, 4)
  expect_identical(result$iter, 321)

  expect_error(
    CoTiMA:::.ctmaResolveInitStanArgs(list(not_a_control = TRUE)),
    "Unknown CoTiMAStanctArgs: not_a_control"
  )
  expect_error(
    CoTiMA:::.ctmaResolveInitStanArgs(list(
      optimcontrol = list(not_a_control = TRUE)
    )),
    "Unknown CoTiMAStanctArgs\\$optimcontrol: not_a_control"
  )
})

test_that("ctmaFit resolves partial legacy Stan controls recursively", {
  result <- CoTiMA:::.ctmaResolveFitStanArgs(
    overrides = list(
      verbose = 2,
      control = list(adapt_delta = 0.95),
      optimcontrol = list(finishsamples = 300)
    )
  )

  expect_identical(result$verbose, 2)
  expect_identical(result$control$adapt_delta, 0.95)
  expect_identical(
    result$control$max_treedepth,
    CoTiMAStanctArgs$control$max_treedepth
  )
  expect_identical(result$optimcontrol$finishsamples, 300)
  expect_identical(
    result$optimcontrol$carefulfit,
    CoTiMAStanctArgs$optimcontrol$carefulfit
  )
})

test_that("ctmaFit direct controls have defined precedence", {
  inherited <- CoTiMA:::.ctmaResolveFitStanArgs(
    overrides = list(verbose = 2, iter = 400, priors = TRUE),
    optimize = TRUE,
    priors = FALSE,
    iter = 900,
    verbose_was_supplied = FALSE
  )
  expect_identical(inherited$verbose, 2)
  expect_identical(inherited$iter, 900)
  expect_false(inherited$priors)

  explicit <- CoTiMA:::.ctmaResolveFitStanArgs(
    overrides = list(verbose = 2),
    verbose = 0,
    verbose_was_supplied = TRUE
  )
  expect_identical(explicit$verbose, 0)
})

test_that("ctmaFit rejects unknown legacy Stan controls", {
  expect_error(
    CoTiMA:::.ctmaResolveFitStanArgs(list(not_a_control = TRUE)),
    "Unknown CoTiMAStanctArgs: not_a_control"
  )
  expect_error(
    CoTiMA:::.ctmaResolveFitStanArgs(list(
      control = list(not_a_control = TRUE)
    )),
    "Unknown CoTiMAStanctArgs\\$control: not_a_control"
  )
})

test_that("ctmaFit applies backend validation before fitting", {
  expect_error(
    ctmaFit(
      ctmaInitFit = NULL,
      CoTiMAStanctArgs = list(control = list(adapt_delta = 0.9))
    ),
    "A fitted CoTiMA object has to be supplied"
  )
  expect_error(
    ctmaFit(
      ctmaInitFit = NULL,
      CoTiMAStanctArgs = list(not_a_control = TRUE)
    ),
    "Unknown CoTiMAStanctArgs: not_a_control"
  )
})

test_that("ctmaFit keeps Stan as its default backend", {
  expect_identical(eval(formals(ctmaFit)$backend), c("stan", "julia"))

  resolved <- CoTiMA:::.ctmaResolveFitBackend()
  expect_identical(resolved$backend, "stan")
  expect_null(resolved$stan_args)
})

test_that("ctmaFit accepts either new or legacy Stan controls", {
  modern <- CoTiMA:::.ctmaResolveFitBackend(
    backendArgs = list(control = list(adapt_delta = 0.9))
  )
  legacy <- CoTiMA:::.ctmaResolveFitBackend(
    legacy_stan_args = list(control = list(adapt_delta = 0.8))
  )
  expect_identical(modern$stan_args$control$adapt_delta, 0.9)
  expect_identical(legacy$stan_args$control$adapt_delta, 0.8)

  expect_error(
    CoTiMA:::.ctmaResolveFitBackend(
      backendArgs = list(iter = 500),
      legacy_stan_args = list(iter = 1000)
    ),
    "either 'CoTiMAStanctArgs' or 'backendArgs'"
  )
  expect_error(
    CoTiMA:::.ctmaResolveFitBackend(
      backendArgs = list(not_a_control = TRUE)
    ),
    "Unknown backendArgs: not_a_control"
  )
})

test_that("ctmaFit reports the staged Julia routing boundary", {
  expect_error(
    CoTiMA:::.ctmaResolveFitBackend("julia"),
    "not yet routed through ctmaFit"
  )
  expect_error(
    CoTiMA:::.ctmaResolveFitBackend(
      "julia", legacy_stan_args = list(iter = 1000)
    ),
    "cannot be used with backend = 'julia'"
  )
})

test_that("ctmaFit backendArgs work through the public entry point", {
  expect_error(
    ctmaFit(
      ctmaInitFit = NULL,
      backend = "stan",
      backendArgs = list(control = list(adapt_delta = 0.9))
    ),
    "A fitted CoTiMA object has to be supplied"
  )
  expect_error(
    ctmaFit(
      ctmaInitFit = NULL,
      CoTiMAStanctArgs = list(iter = 1000),
      backendArgs = list(iter = 500)
    ),
    "either 'CoTiMAStanctArgs' or 'backendArgs'"
  )
  expect_error(
    ctmaFit(ctmaInitFit = NULL, backend = "julia"),
    "ctmaJuliaFit"
  )
})

test_that("ctmaAllInvFit resolves partial nested Stan controls", {
  result <- CoTiMA:::.ctmaResolveAllInvStanArgs(
    overrides = list(
      verbose = 2,
      iter = 400,
      control = list(adapt_delta = 0.95),
      optimcontrol = list(finishsamples = 250)
    ),
    scaleTime = 0.5,
    iter = 900,
    finishsamples = 500
  )

  expect_identical(result$verbose, 2)
  expect_identical(result$scaleTime, 0.5)
  expect_identical(result$iter, 900)
  expect_identical(result$control$adapt_delta, 0.95)
  expect_identical(
    result$control$max_treedepth,
    CoTiMAStanctArgs$control$max_treedepth
  )
  expect_identical(result$optimcontrol$finishsamples, 500)
  expect_identical(
    result$optimcontrol$carefulfit,
    CoTiMAStanctArgs$optimcontrol$carefulfit
  )
})

test_that("ctmaAllInvFit explicit verbose overrides the legacy list", {
  inherited <- CoTiMA:::.ctmaResolveAllInvStanArgs(
    overrides = list(verbose = 2),
    verbose_was_supplied = FALSE
  )
  explicit <- CoTiMA:::.ctmaResolveAllInvStanArgs(
    overrides = list(verbose = 2),
    verbose = 0,
    verbose_was_supplied = TRUE
  )
  expect_identical(inherited$verbose, 2)
  expect_identical(explicit$verbose, 0)
})

test_that("ctmaAllInvFit validates controls before requiring a fit", {
  expect_error(
    ctmaAllInvFit(
      ctmaInitFit = NULL,
      CoTiMAStanctArgs = list(control = list(adapt_delta = 0.9))
    ),
    "A fitted CoTiMA object has to be supplied"
  )
  expect_error(
    ctmaAllInvFit(
      ctmaInitFit = NULL,
      CoTiMAStanctArgs = list(not_a_control = TRUE)
    ),
    "Unknown CoTiMAStanctArgs: not_a_control"
  )
})

test_that("ctmaPower resolves partial nested Stan controls", {
  result <- CoTiMA:::.ctmaResolvePowerStanArgs(
    overrides = list(
      scaleTime = 0.5,
      verbose = 2,
      control = list(adapt_delta = 0.95),
      optimcontrol = list(finishsamples = 250)
    ),
    iter = 900,
    finishsamples = 500
  )

  expect_identical(result$scaleTime, 0.5)
  expect_identical(result$verbose, 2)
  expect_identical(result$iter, 900)
  expect_identical(result$control$adapt_delta, 0.95)
  expect_identical(
    result$control$max_treedepth,
    CoTiMAStanctArgs$control$max_treedepth
  )
  expect_identical(result$optimcontrol$finishsamples, 500)
  expect_identical(
    result$optimcontrol$carefulfit,
    CoTiMAStanctArgs$optimcontrol$carefulfit
  )
})

test_that("ctmaPower direct scale and verbose controls take precedence", {
  result <- CoTiMA:::.ctmaResolvePowerStanArgs(
    overrides = list(scaleTime = 0.5, verbose = 2),
    scaleTime = 2,
    verbose = 0,
    verbose_was_supplied = TRUE
  )
  expect_identical(result$scaleTime, 2)
  expect_identical(result$verbose, 0)
})

test_that("ctmaPower validates controls before requiring a fit", {
  expect_error(
    ctmaPower(
      ctmaInitFit = NULL,
      statisticalPower = 0.8,
      CoTiMAStanctArgs = list(control = list(adapt_delta = 0.9))
    ),
    "A fitted CoTiMA object.*has to be supplied"
  )
  expect_error(
    ctmaPower(
      ctmaInitFit = NULL,
      statisticalPower = 0.8,
      CoTiMAStanctArgs = list(not_a_control = TRUE)
    ),
    "Unknown CoTiMAStanctArgs: not_a_control"
  )
})

test_that("ctmaOptimizeFit merges defaults, inherited controls, and overrides", {
  result <- CoTiMA:::.ctmaResolveOptimizeStanArgs(
    inherited = list(
      verbose = 2,
      iter = 700,
      control = list(adapt_delta = 0.9),
      optimcontrol = list(finishsamples = 300)
    ),
    overrides = list(
      iter = 900,
      control = list(max_treedepth = 12L)
    ),
    finishsamples = 500
  )

  expect_identical(result$verbose, 2)
  expect_identical(result$iter, 900)
  expect_identical(result$control$adapt_delta, 0.9)
  expect_identical(result$control$max_treedepth, 12L)
  expect_identical(result$optimcontrol$finishsamples, 500)
  expect_identical(
    result$optimcontrol$carefulfit,
    CoTiMAStanctArgs$optimcontrol$carefulfit
  )
})

test_that("ctmaOptimizeFit rejects unknown stored and supplied controls", {
  expect_error(
    CoTiMA:::.ctmaResolveOptimizeStanArgs(
      inherited = list(not_a_control = TRUE)
    ),
    "Unknown stored CoTiMAStanctArgs: not_a_control"
  )
  expect_error(
    CoTiMA:::.ctmaResolveOptimizeStanArgs(
      overrides = list(not_a_control = TRUE)
    ),
    "Unknown CoTiMAStanctArgs: not_a_control"
  )
})

test_that("ctmaOptimizeFit recognizes only finite fit criteria", {
  valid <- list(summary = list(minus2ll = 123.4))
  nonfinite <- list(summary = list(minus2ll = Inf))
  missing <- list(summary = list())
  failed <- simpleError("failed")

  expect_identical(CoTiMA:::.ctmaOptimizeFitValue(valid), 123.4)
  expect_true(is.na(CoTiMA:::.ctmaOptimizeFitValue(nonfinite)))
  expect_true(is.na(CoTiMA:::.ctmaOptimizeFitValue(missing)))
  expect_true(is.na(CoTiMA:::.ctmaOptimizeFitValue(failed)))
})

test_that("ctmaOptimizeFit validates backend controls before input objects", {
  expect_error(
    ctmaOptimizeFit(
      CoTiMAStanctArgs = list(control = list(adapt_delta = 0.9))
    ),
    "ctmaFitFit or ctmaIniFit"
  )
  expect_error(
    ctmaOptimizeFit(CoTiMAStanctArgs = list(not_a_control = TRUE)),
    "Unknown CoTiMAStanctArgs: not_a_control"
  )
})

test_that("parallel fit selection uses the greatest finite log-posterior", {
  attempt <- function(value) {
    list(
      ok = TRUE,
      fit = list(stanfit = list(optimfit = list(value = value))),
      warn_hessinv = FALSE,
      warnings = character(),
      error = NA_character_
    )
  }

  result <- CoTiMA:::ctmaSelectBestParallelFit(list(
    attempt(-100),
    attempt(-10),
    attempt(NA_real_)
  ))

  expect_identical(result$logposterior, -10)
  expect_identical(result$fit$stanfit$optimfit$value, -10)
})

test_that("parallel fit selection ignores errors and reports total failure", {
  failed <- list(
    ok = FALSE,
    fit = NULL,
    warn_hessinv = FALSE,
    warnings = character(),
    error = "fit failed"
  )
  successful <- list(
    ok = TRUE,
    fit = list(stanfit = list(optimfit = list(value = -5))),
    warn_hessinv = FALSE,
    warnings = character(),
    error = NA_character_
  )

  expect_identical(
    CoTiMA:::ctmaSelectBestParallelFit(list(failed, successful))$logposterior,
    -5
  )
  expect_error(
    CoTiMA:::ctmaSelectBestParallelFit(list(failed)),
    "All parallel fitting attempts failed"
  )
})
