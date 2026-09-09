test_that("native ctsem Julia/Laplace capability detection is structured", {
  status <- .ctmaCtsemJuliaLaplaceStatus()

  expect_s3_class(status, "ctmaCtsemJuliaLaplaceStatus")
  expect_type(status$available, "logical")
  expect_length(status$available, 1L)
  expect_named(
    status$api,
    c(
      "ctLaplaceCheck", "ctLaplaceCorrect", "ctJuliaStatus", "ctJuliaSetup",
      "outerLevelOnly"
    )
  )
  expect_identical(status$available, all(status$api) && status$JuliaConnectoR)
  if (!status$available) expect_type(status$reason, "character")
})

test_that("study effects do not require a subject-level random effect", {
  model <- list(
    subjectIDname = "id",
    groupIDnames = character(),
    TIpredNames = character(),
    pars = data.frame(
      param = c("cint_eta1", "diff_eta1"),
      indvarying = c(FALSE, FALSE),
      sdscale = c(1, 1),
      stringsAsFactors = FALSE
    )
  )
  setup <- structure(list(
    data = data.frame(id = rep(1:2, each = 2), y = 1:4),
    ctModel = model
  ), class = "CoTiMAFit")

  local_mocked_bindings(
    .ctmaRequireCtsemJuliaLaplace = function() list(available = TRUE),
    .ctma_is_unfitted_cotima = function(object) TRUE,
    .ctmaLaplaceStudyMembership = function(...) list(
      row_id = c(1L, 1L, 2L, 2L),
      labels = c("A", "B"),
      dummy_names = character()
    ),
    .ctma_generate_ctsem_setup = function(...) list(),
    .package = "CoTiMA2"
  )

  prepared <- .ctmaPrepareCtsemJuliaLaplace(
    setup,
    study_parameters = "cint_eta1",
    verbose = FALSE
  )

  expect_false(any(prepared$model$pars$indvarying))
  expect_true(prepared$model$pars$indvarying_study[1])
  expect_false(prepared$model$pars$indvarying_study[2])
  expect_identical(prepared$model$groupIDnames, "study")
})

test_that("matrix names select study and subject random effects", {
  model <- list(
    subjectIDname = "id",
    groupIDnames = character(),
    TIpredNames = character(),
    pars = data.frame(
      matrix = c("CINT", "DRIFT", "DIFFUSION"),
      param = c("cint_eta1", "drift_eta1", "diff_eta1"),
      indvarying = c(FALSE, FALSE, FALSE),
      sdscale = c(1, 1, 1),
      stringsAsFactors = FALSE
    )
  )
  setup <- structure(list(
    data = data.frame(id = rep(1:2, each = 2), y = 1:4),
    ctModel = model
  ), class = "CoTiMAFit")

  local_mocked_bindings(
    .ctmaRequireCtsemJuliaLaplace = function() list(available = TRUE),
    .ctma_is_unfitted_cotima = function(object) TRUE,
    .ctmaLaplaceStudyMembership = function(...) list(
      row_id = c(1L, 1L, 2L, 2L), labels = c("A", "B"),
      dummy_names = character()
    ),
    .ctma_generate_ctsem_setup = function(...) list(),
    .package = "CoTiMA2"
  )

  prepared <- .ctmaPrepareCtsemJuliaLaplace(
    setup,
    study_random_effects = c("CINT", "DRIFT"),
    subject_random_effects = "DIFFUSION",
    verbose = FALSE
  )

  expect_identical(
    prepared$selected_parameters,
    c("cint_eta1", "drift_eta1")
  )
  expect_identical(
    prepared$model$pars$indvarying,
    c(FALSE, FALSE, TRUE)
  )
  expect_identical(
    prepared$model$pars$indvarying_study,
    c(TRUE, TRUE, FALSE)
  )
  expect_error(
    .ctmaPrepareCtsemJuliaLaplace(
      setup, study_random_effects = "NOT_A_MATRIX", verbose = FALSE
    ),
    "Unknown study-level"
  )
})

test_that("the public ctsemjulia wrapper retains CoTiMA provenance", {
  setup <- structure(list(
    argumentList = list(fit = FALSE),
    marker = "source"
  ), class = "CoTiMAFit")
  prepared <- structure(list(marker = "prepared"),
    class = "ctmaCtsemJuliaLaplacePrep")

  local_mocked_bindings(
    .ctmaPrepareCtsemJuliaLaplace = function(...) prepared,
    .ctmaFitCtsemJuliaLaplace = function(...) list(
      estimate = list(loglik = -10, converged = TRUE)
    ),
    .package = "CoTiMA2"
  )

  result <- .ctmaFitCtsemJuliaBackend(
    setup = setup,
    backend_args = ctmaBackendDefaults("ctsemjulia"),
    cores = 2L,
    cotima_call = quote(ctmaFit(backend = "ctsemjulia"))
  )

  expect_identical(result$backend, "ctsemjulia")
  expect_identical(result$source_object, setup)
  expect_identical(result$ctsemjulia_prepared, prepared)
  expect_identical(result$argumentList$backend, "ctsemjulia")
  expect_true(result$argumentList$fit)
})

test_that("ctmaFit advertises and dispatches the ctsemjulia backend", {
  expect_true("ctsemjulia" %in% eval(formals(ctmaFit)$backend))
  expect_true(any(grepl(
    ".ctmaFitCtsemJuliaBackend",
    deparse(body(ctmaFit)),
    fixed = TRUE
  )))
})

test_that("explicit study membership must be strictly nested", {
  setup <- structure(list(
    data = data.frame(id = c(1, 1, 2, 2), y = 1:4),
    ctModel = list(subjectIDname = "id", TIpredNames = "TI1")
  ), class = "CoTiMAFit")

  expect_error(
    .ctmaLaplaceStudyMembership(
      setup,
      study_id = c("A", "B", "A", "A"),
      verbose = FALSE
    ),
    "exactly one study"
  )

  membership <- .ctmaLaplaceStudyMembership(
    setup,
    study_id = c("A", "A", "B", "B"),
    study_dummy_ti = "TI1",
    verbose = FALSE
  )
  expect_identical(membership$row_id, c(1L, 1L, 2L, 2L))
  expect_identical(membership$labels, c("A", "B"))
  expect_identical(membership$dummy_names, "TI1")
})
