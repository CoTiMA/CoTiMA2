test_that("Julia backend resources are installed with CoTiMA", {
  julia_file <- system.file("julia", "CtsemJulia.jl", package = "CoTiMA2")
  project_file <- system.file("julia", "Project.toml", package = "CoTiMA2")

  expect_true(nzchar(julia_file))
  expect_true(file.exists(julia_file))
  expect_true(nzchar(project_file))
  expect_true(file.exists(project_file))
})

test_that("Julia user-facing functions and controls are available", {
  expect_true(is.function(ctmaJuliaFit))
  expect_true(is.function(ctAddJuliaSamples))
  expect_identical(ctsmAddJuliaSamples, ctAddJuliaSamples)

  fit_arguments <- names(formals(ctmaJuliaFit))
  expect_true(all(c(
    "estimation", "inits", "carefulFit", "multistart",
    "parallel_starts", "threads_per_start", "parallel_chains",
    "threads_per_chain", "warm_start", "warm_start_fraction",
    "warm_start_subjects_per_parameter", "warm_start_minimum_per_study",
    "warm_start_maximum_subjects", "warm_start_binary_minimum",
    "warm_start_maxiters"
  ) %in% fit_arguments))
})

test_that("the packaged Julia project declares the backend dependencies", {
  project_file <- system.file("julia", "Project.toml", package = "CoTiMA2")
  backend_file <- system.file("julia", "CtsemJulia.jl", package = "CoTiMA2")
  project <- readLines(project_file, warn = FALSE)
  backend <- readLines(backend_file, warn = FALSE)

  expect_true(any(grepl("AdvancedHMC", project, fixed = TRUE)))
  expect_true(any(grepl("ForwardDiff", project, fixed = TRUE)))
  expect_true(any(grepl("Optim", project, fixed = TRUE)))
  expect_true(any(grepl("using ForwardDiff", backend, fixed = TRUE)))
  expect_true(any(grepl("using Optim", backend, fixed = TRUE)))
  expect_true(any(grepl("validate_parallel_subjects", backend, fixed = TRUE)))
  expect_true(any(grepl("ctsem_optimize_warm_start", backend, fixed = TRUE)))
  expect_true(any(grepl("Threads.@threads :static", backend, fixed = TRUE)))
})

test_that("model-derived Julia starts respect ctsem matrix transforms", {
  testthat::local_mocked_bindings(
    .ctmaJuliaRawParameterCount = function(fit) 7L,
    .package = "CoTiMA2"
  )
  fit <- list(setup = list(popsetup = data.frame(
    param = 1:7,
    transform = c(1L, 0L, 1L, 1L, 3L, 1L, 0L),
    row.names = c(
      "DRIFT", "DRIFT1", "DIFFUSION", "T0VAR",
      "T0VAR1", "MANIFESTVAR", "CINT"
    )
  )))

  expect_equal(
    CoTiMA2:::.ctmaJuliaDefaultRawStart(fit),
    c(0.9, 0, -1.2, -0.75, 0, -1.2, 0)
  )
})
