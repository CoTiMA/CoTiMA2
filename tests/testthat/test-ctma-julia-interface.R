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
  expect_true(is.function(ctmaTransformJuliaFit))
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

test_that("point transformation validates its public inputs", {
  expect_error(
    ctmaTransformJuliaFit(list(), setup = list()),
    "must inherit from 'ctmaJuliaFit'",
    fixed = TRUE
  )
  empty_fit <- structure(list(), class = "ctmaJuliaFit")
  expect_error(
    ctmaTransformJuliaFit(empty_fit, setup = list(), cores = 0L),
    "'cores' must be a positive integer",
    fixed = TRUE
  )
  expect_error(
    ctmaTransformJuliaFit(empty_fit, setup = list(), digits = 0L),
    "'digits' must be a positive integer",
    fixed = TRUE
  )
})

test_that("point transformation does not require uncertainty information", {
  captured_rows <- NULL
  testthat::local_mocked_bindings(
    .ctmaJuliaReadObject = identity,
    .ctmaJuliaIsUnfittedCoTiMA = function(object) FALSE,
    .ctmaJuliaExtractFit = function(object) list(stanfit = list(rawest = 1:2)),
    .ctmaJuliaAttachConstraintModel = identity,
    .ctmaJuliaTransformSummary = function(
        fit, estimates, covariance, loglik, draws, cores, seed = NULL,
        raw_draws = NULL, digits = 4L) {
      captured_rows <<- raw_draws
      expect_equal(covariance, matrix(0, 2, 2))
      expect_identical(draws, 2L)
      list(
        parameters = "transformed point",
        summary = list(
          popmeans = "population point",
          tipreds = "TI point",
          parmatrices = "matrix point"
        ),
        ctsem_fit = "constraint bridge"
      )
    },
    .ctmaJuliaCoTiMAWEC = function(object, fit, raw_draws) {
      list(study_specific_parameters = "study point")
    },
    .package = "CoTiMA2"
  )
  fit <- structure(
    list(
      coefficients = list(raw = c(0.25, -0.5)),
      optimization = list(loglik = -10),
      information = NULL
    ),
    class = "ctmaJuliaFit"
  )

  result <- ctmaTransformJuliaFit(fit, setup = list())

  expect_equal(captured_rows, rbind(c(0.25, -0.5), c(0.25, -0.5)))
  expect_identical(result$coefficients$transformed, "transformed point")
  expect_identical(result$coefficients$matrices, "matrix point")
  expect_identical(result$coefficients$study, "study point")
  expect_null(result$transformed_draws)
  expect_null(result$ctsem_summary)
  expect_false(result$point_transformation$uncertainty)
  expect_false(result$point_transformation$draws)
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
  expect_true(any(grepl("_nuts_safe_evaluate", backend, fixed = TRUE)))
  expect_true(any(grepl("_nuts_rejectable_error", backend, fixed = TRUE)))
  expect_true(any(grepl("CompositeException", backend, fixed = TRUE)))
  expect_true(any(grepl("return (-Inf, zeros", backend, fixed = TRUE)))
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
test_that("Julia setup suppresses only misleading backend messages", {
  expect_true(.ctmaJuliaSuppressSetupMessage(
    "Maximum likelihood estimation requested"
  ))
  expect_true(.ctmaJuliaSuppressSetupMessage(paste(
    "Binary indicators use a linearised (moment-matched Gaussian)",
    "measurement update on the stan backend, which makes DRIFT unreliable."
  )))
  expect_false(.ctmaJuliaSuppressSetupMessage(
    "Binary likelihood validation failed."
  ))
  expect_false(.ctmaJuliaSuppressSetupMessage(
    "A genuine ctsem setup warning"
  ))
})

test_that("bundled Julia asks ctsem for the Stan layout without estimation", {
  object <- list(data = data.frame(id = 1L), ctModel = list())
  setup <- function(ct_fit) {
    testthat::local_mocked_bindings(ctFit = ct_fit, .package = "ctsem")
    suppressMessages(.ctmaJuliaGenerateSetup(object))
  }

  modern <- setup(function(fit, backend = "julia", ...) {
    list(fit = fit, backend = backend)
  })
  legacy <- setup(function(fit, ...) list(fit = fit, extras = list(...)))

  expect_identical(modern, list(fit = FALSE, backend = "stan"))
  expect_false(legacy$fit)
  expect_false("backend" %in% names(legacy$extras))
})
