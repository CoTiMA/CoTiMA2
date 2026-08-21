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
  package_ml <- ctmaBackendDefaults("julia", "ML", "ctmaFit")

  expect_identical(ml$carefulFit, TRUE)
  expect_identical(ml$multistart, 1L)
  expect_false("chains" %in% names(ml))
  expect_identical(bayes$metric, "diagonal")
  expect_identical(bayes$chains, 4L)
  expect_false("maxiters" %in% names(bayes))
  expect_identical(re_bayes$metric, "dense")
  expect_identical(re_bayes$progress_every, 25L)
  expect_true(package_ml$uncertainty)
  expect_true(package_ml$transformed)

  expect_error(
    ctmaBackendDefaults("julia", "ML", "ctmaReFit"),
    "currently support only"
  )
})

test_that("ctmaAllInvFit exposes Julia backend defaults", {
  defaults <- ctmaBackendDefaults("julia", "ML", "ctmaAllInvFit")
  expect_true(defaults$uncertainty)
  expect_true(defaults$transformed)
  expect_true(all(c("backend", "backendArgs") %in%
                    names(formals(ctmaAllInvFit))))
})

test_that("central core resolution is backend aware", {
  stan <- CoTiMA2:::.ctmaResolveCores(
    -1L, backend = "stan", max_stan_workers = 8L, announce = FALSE
  )
  expect_gte(stan$resolved, 1L)
  expect_lte(stan$resolved, 8L)

  tasks <- CoTiMA2:::.ctmaResolveCores(
    12L, backend = "stan", task_count = 3L,
    max_stan_workers = 32L, announce = FALSE
  )
  expect_identical(tasks$resolved, min(3L, tasks$available))

  julia <- CoTiMA2:::.ctmaResolveCores(
    -1L, backend = "julia", max_stan_workers = 1L, announce = FALSE
  )
  expect_identical(julia$resolved, max(1L, julia$available - 1L))
})

test_that("backend overrides merge recursively and by name", {
  defaults <- list(
    iter = 1000L,
    verbose = 0L,
    control = list(adapt_delta = 0.8, max_treedepth = 10L)
  )
  result <- CoTiMA2:::.ctmaMergeBackendArgs(defaults, list(
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
    CoTiMA2:::.ctmaMergeBackendArgs(defaults, list()),
    defaults
  )
  expect_identical(
    CoTiMA2:::.ctmaMergeBackendArgs(defaults, list(control = list())),
    defaults
  )
})

test_that("backend overrides reject unnamed and unknown controls", {
  defaults <- list(iter = 1000L, control = list(adapt_delta = 0.8))

  expect_error(
    CoTiMA2:::.ctmaMergeBackendArgs(defaults, list(500L)),
    "fully named list"
  )
  expect_error(
    CoTiMA2:::.ctmaMergeBackendArgs(defaults, list(unknown = 1)),
    "Unknown backendArgs: unknown"
  )
  expect_error(
    CoTiMA2:::.ctmaMergeBackendArgs(
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
test_that("CoTiMA2 preview identity and major-call notices are explicit", {
  expect_identical(as.character(utils::packageVersion("CoTiMA2")), "0.11.0")
  expect_message(
    CoTiMA2:::.ctmaVersionNote(),
    "Note: This is the experimental CoTiMA 2.0 package \\(version 0.11.0\\)."
  )

  major_functions <- c(
    "ctmaPrep", "ctmaInit", "ctmaFit", "ctmaEqual", "ctmaOptimizeFit",
    "ctmaPower", "ctmaBiG", "ctmaPub", "ctmaRedHet", "ctmaSV",
    "ctmaJuliaFit", "ctmaJuliaSetup", "ctAddJuliaSamples", "ctmaRePrep",
    "ctmaReFit", "ctmaReMLFit", "ctmaReMoments", "ctmaReRedHet"
  )
  for (function_name in major_functions) {
    function_body <- deparse(body(getExportedValue("CoTiMA2", function_name)))
    expect_true(
      any(grepl(".ctmaVersionNote()", function_body, fixed = TRUE)),
      info = function_name
    )
    expect_true(
      any(grepl(".ctmaRegisterExecutionTimer", function_body, fixed = TRUE)),
      info = paste(function_name, "execution timer")
    )
  }
})

test_that("execution timer reports seconds for short calls", {
  timer <- CoTiMA2:::.ctmaRegisterExecutionTimer("exampleFunction")
  expect_message(
    elapsed <- timer(),
    paste0(
      "exampleFunction total wall-clock execution time: ",
      "[0-9]+\\.[0-9]{3} seconds\\."
    )
  )
  expect_true(is.numeric(elapsed))
  expect_gte(elapsed, 0)
})
