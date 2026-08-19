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
