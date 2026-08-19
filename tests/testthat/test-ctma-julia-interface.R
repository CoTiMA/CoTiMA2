test_that("Julia backend resources are installed with CoTiMA", {
  julia_file <- system.file("julia", "CtsemJulia.jl", package = "CoTiMA")
  project_file <- system.file("julia", "Project.toml", package = "CoTiMA")

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
    "threads_per_chain"
  ) %in% fit_arguments))
})

test_that("the packaged Julia project declares the backend dependencies", {
  project_file <- system.file("julia", "Project.toml", package = "CoTiMA")
  backend_file <- system.file("julia", "CtsemJulia.jl", package = "CoTiMA")
  project <- readLines(project_file, warn = FALSE)
  backend <- readLines(backend_file, warn = FALSE)

  expect_true(any(grepl("AdvancedHMC", project, fixed = TRUE)))
  expect_true(any(grepl("using ForwardDiff", backend, fixed = TRUE)))
  expect_true(any(grepl("using Optim", backend, fixed = TRUE)))
})
