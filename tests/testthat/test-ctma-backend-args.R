test_that("Stan is the default backend and preserves legacy defaults", {
  expect_identical(ctmaBackendDefaults(), CoTiMAStanctArgs)
  expect_identical(
    ctmaBackendDefaults("stan", "Bayes", "ctmaReFit"),
    CoTiMAStanctArgs
  )
})

test_that("Julia defaults are workflow and estimation specific", {
  ml <- ctmaBackendDefaults("julia", "ML", "ctmaJuliaFit")
  bayes <- ctmaBackendDefaults("julia", "Bayes", "ctmaJuliaFit")
  re_bayes <- ctmaBackendDefaults("julia", "Bayes", "ctmaReFit")

  expect_identical(ml$carefulFit, TRUE)
  expect_identical(ml$multistart, 1L)
  expect_false("chains" %in% names(ml))
  expect_identical(bayes$metric, "diagonal")
  expect_identical(bayes$chains, 4L)
  expect_false("maxiters" %in% names(bayes))
  expect_identical(re_bayes$metric, "dense")
  expect_identical(re_bayes$progress_every, 25L)

  expect_error(
    ctmaBackendDefaults("julia", "ML", "ctmaReFit"),
    "currently support only"
  )
})

test_that("backend overrides merge recursively and by name", {
  defaults <- list(
    iter = 1000L,
    verbose = 0L,
    control = list(adapt_delta = 0.8, max_treedepth = 10L)
  )
  result <- CoTiMA:::.ctmaMergeBackendArgs(defaults, list(
    control = list(max_treedepth = 12L),
    iter = 500L
  ))

  expect_identical(result$iter, 500L)
  expect_identical(result$verbose, 0L)
  expect_identical(result$control$adapt_delta, 0.8)
  expect_identical(result$control$max_treedepth, 12L)
})

test_that("empty top-level and nested backend lists are valid", {
  defaults <- list(control = list(), iter = 1000L)
  expect_identical(
    CoTiMA:::.ctmaMergeBackendArgs(defaults, list()),
    defaults
  )
  expect_identical(
    CoTiMA:::.ctmaMergeBackendArgs(defaults, list(control = list())),
    defaults
  )
})

test_that("backend overrides reject unnamed and unknown controls", {
  defaults <- list(iter = 1000L, control = list(adapt_delta = 0.8))

  expect_error(
    CoTiMA:::.ctmaMergeBackendArgs(defaults, list(500L)),
    "fully named list"
  )
  expect_error(
    CoTiMA:::.ctmaMergeBackendArgs(defaults, list(unknown = 1)),
    "Unknown backendArgs: unknown"
  )
  expect_error(
    CoTiMA:::.ctmaMergeBackendArgs(
      defaults, list(control = list(unknown = 1))
    ),
    "Unknown backendArgs\\$control: unknown"
  )
})

test_that("backend defaults are isolated between calls", {
  first <- ctmaBackendDefaults("julia", "ML", "ctmaJuliaFit")
  first$maxiters <- 999L
  second <- ctmaBackendDefaults("julia", "ML", "ctmaJuliaFit")
  expect_identical(second$maxiters, 100L)

  stan_first <- ctmaBackendDefaults()
  stan_first$control$adapt_delta <- 0.99
  expect_identical(ctmaBackendDefaults()$control$adapt_delta, 0.8)
})
