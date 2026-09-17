test_that("automatic Julia capacity reserves one available CPU", {
  expect_identical(
    CoTiMA2:::.ctmaResolveJuliaThreads("auto", available = 64L),
    63L
  )
  expect_identical(
    CoTiMA2:::.ctmaResolveJuliaThreads("auto", available = 1L),
    1L
  )
  expect_identical(CoTiMA2:::.ctmaResolveJuliaThreads(8L), 8L)
  expect_error(CoTiMA2:::.ctmaResolveJuliaThreads(0L), "positive integer")
})

test_that("subject concurrency is bounded by Julia pool capacity", {
  expect_identical(CoTiMA2:::.ctmaResolveSubjectThreads("auto", 16L), 16L)
  expect_identical(CoTiMA2:::.ctmaResolveSubjectThreads(8L, 16L), 8L)
  expect_identical(CoTiMA2:::.ctmaResolveSubjectThreads(32L, 16L), 16L)
  expect_error(CoTiMA2:::.ctmaResolveSubjectThreads(0L, 16L), "positive integer")
})

test_that("ctmaFitConfig exposes all backend routes and tuning controls", {
  expect_true(is.function(ctmaFitConfig))
  arguments <- names(formals(ctmaFitConfig))
  expect_true(all(c(
    "ctmaInitFit", "...", "backend", "backendArgs", "time_limit",
    "candidate_subject_threads", "candidate_chunksizes",
    "candidate_carefulFit", "pilot_iterations", "benchmark_repetitions"
  ) %in% arguments))
  expect_identical(
    eval(formals(ctmaFitConfig)$backend),
    c("julia", "stan", "ctsemjulia")
  )
})

test_that("ctsemjulia configuration validates and retains native controls", {
  checked <- FALSE
  local_mocked_bindings(
    .ctmaCheckCtsemJuliaRequest = function(backend, workflow, supported = FALSE) {
      checked <<- identical(backend, "ctsemjulia") &&
        identical(workflow, "ctmaFitConfig") && isTRUE(supported)
      invisible(NULL)
    },
    ctmaFit = function(...) list(prepared = TRUE),
    .ctmaResolveFitBackend = function(...) list(
      backend_args = list(
        inits = NULL,
        optimcontrol = list(uncertainty = "hessian", maxiter = 250L)
      )
    ),
    .package = "CoTiMA2"
  )

  config <- ctmaFitConfig(
    backend = "ctsemjulia",
    backendArgs = list(optimcontrol = list(maxiter = 250L)),
    verbose = FALSE
  )

  expect_true(checked)
  expect_s3_class(config, "ctmaFitConfig")
  expect_identical(config$backend, "ctsemjulia")
  expect_identical(config$setup, list(prepared = TRUE))
  expect_identical(config$recommended, config$resolved_backend_args)
  expect_identical(config$recommended$optimcontrol$uncertainty, "hessian")
  expect_identical(config$recommended$optimcontrol$maxiter, 250L)
  expect_equal(nrow(config$benchmarks), 0L)
  expect_match(config$note, "native ctsem Julia/Laplace backend", fixed = TRUE)
  expect_match(config$note, "coresToUse", fixed = TRUE)
})

test_that("Julia configuration uses repeat timing and convergence-aware pilots", {
  body_text <- paste(deparse(body(ctmaFitConfig)), collapse = "\n")
  expect_match(body_text, "stats::median", fixed = TRUE)
  expect_match(body_text, "gradient_reduction_factor", fixed = TRUE)
  expect_match(body_text, "convergence_adjusted_score", fixed = TRUE)
  expect_match(body_text, "deterministic serial path only", fixed = TRUE)
})
