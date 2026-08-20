test_that("Julia executable discovery accepts explicit paths", {
  executable <- file.path(R.home("bin"), "R")
  expect_identical(
    CoTiMA2:::.ctmaJuliaFindExecutable(julia = executable),
    normalizePath(executable)
  )
  expect_error(
    CoTiMA2:::.ctmaJuliaFindExecutable(julia = tempfile()),
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
  expect_error(
    ctmaJuliaSetup(
      julia = executable, bootstrap_juliacall = "sometimes", dry_run = TRUE
    ),
    '"auto", TRUE, or FALSE',
    fixed = TRUE
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

test_that("Linux Julia subprocesses discard inherited native library paths", {
  command <- list(
    label = "Test Julia",
    executable = "/opt/julia/bin/julia",
    args = c("--startup-file=no", "--version"),
    env = c("R_HOME=/usr/lib/R", "JULIA_PKG_PRECOMPILE_AUTO=0")
  )
  invocation <- CoTiMA2:::.ctmaJuliaCommandInvocation(command, linux = TRUE)

  expect_match(invocation$executable, "env$")
  expect_identical(invocation$label, "Test Julia")
  expect_identical(
    invocation$args[seq_len(4L)],
    c("-u", "LD_LIBRARY_PATH", "-u", "LD_PRELOAD")
  )
  expect_true("R_HOME=/usr/lib/R" %in% invocation$args)
  expect_true("'/opt/julia/bin/julia'" %in% invocation$args)
  expect_length(invocation$env, 0L)

  unchanged <- CoTiMA2:::.ctmaJuliaCommandInvocation(command, linux = FALSE)
  expect_identical(unchanged, command)
})

test_that("JuliaCall bootstrap can be skipped or forced", {
  executable <- file.path(R.home("bin"), "R")
  skipped <- ctmaJuliaSetup(
    julia = executable,
    bootstrap_juliacall = FALSE,
    dry_run = TRUE,
    verbose = FALSE
  )
  forced <- ctmaJuliaSetup(
    julia = executable,
    bootstrap_juliacall = TRUE,
    dry_run = TRUE,
    verbose = FALSE
  )

  expect_length(skipped$commands, 2L)
  expect_length(forced$commands, 3L)
  expect_match(forced$commands[[1L]], 'Pkg.build("RCall")', fixed = TRUE)
})

test_that("fitting interfaces default to automatic Julia discovery", {
  expect_null(eval(formals(ctmaJuliaFit)$julia_home))
  expect_null(eval(formals(ctmaJuliaSetup)$julia_home))
  expect_identical(eval(formals(ctmaJuliaSetup)$bootstrap_juliacall), "auto")
})
