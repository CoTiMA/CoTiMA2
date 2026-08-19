test_that("Julia executable discovery accepts explicit paths", {
  executable <- file.path(R.home("bin"), "R")
  expect_identical(
    CoTiMA:::.ctmaJuliaFindExecutable(julia = executable),
    normalizePath(executable)
  )
  expect_error(
    CoTiMA:::.ctmaJuliaFindExecutable(julia = tempfile()),
    "Julia was not found"
  )
})

test_that("ctmaJuliaSetup validates controls without executing Julia", {
  executable <- file.path(R.home("bin"), "R")

  expect_error(
    ctmaJuliaSetup(julia = executable, precompile_tasks = 0, dry_run = TRUE),
    "positive integer"
  )
  expect_error(
    ctmaJuliaSetup(julia = executable, project = tempfile(), dry_run = TRUE),
    "project was not found"
  )
})

test_that("ctmaJuliaSetup constructs a portable dry-run plan", {
  executable <- file.path(R.home("bin"), "R")
  result <- ctmaJuliaSetup(
    julia = executable,
    dry_run = TRUE,
    verbose = FALSE,
    precompile_tasks = 1
  )

  expect_s3_class(result, "ctmaJuliaSetup")
  expect_identical(result$julia, normalizePath(executable))
  expect_true(dir.exists(result$project))
  expect_length(result$commands, 3L)
  expect_match(result$commands[[1L]], "Suppressor")
  expect_match(result$commands[[1L]], "RCall")
  expect_match(result$commands[[1L]], "JULIA_PKG_PRECOMPILE_AUTO=0")
  expect_false(grepl("Pkg.precompile", result$commands[[1L]], fixed = TRUE))
  expect_match(result$commands[[2L]], "allow_autoprecomp=false", fixed = TRUE)
  expect_match(result$commands[[3L]], "JULIA_NUM_PRECOMPILE_TASKS=1")
})
