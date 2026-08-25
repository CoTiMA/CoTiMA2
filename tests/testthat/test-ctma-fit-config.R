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

test_that("ctmaFitConfig exposes both backend routes and tuning controls", {
  expect_true(is.function(ctmaFitConfig))
  arguments <- names(formals(ctmaFitConfig))
  expect_true(all(c(
    "ctmaInitFit", "...", "backend", "backendArgs", "time_limit",
    "candidate_subject_threads", "candidate_chunksizes",
    "candidate_carefulFit", "pilot_iterations", "benchmark_repetitions"
  ) %in% arguments))
  expect_identical(eval(formals(ctmaFitConfig)$backend), c("julia", "stan"))
})

test_that("Julia configuration uses repeat timing and convergence-aware pilots", {
  body_text <- paste(deparse(body(ctmaFitConfig)), collapse = "\n")
  expect_match(body_text, "stats::median", fixed = TRUE)
  expect_match(body_text, "gradient_reduction_factor", fixed = TRUE)
  expect_match(body_text, "convergence_adjusted_score", fixed = TRUE)
  expect_match(body_text, "deterministic serial path only", fixed = TRUE)
})
