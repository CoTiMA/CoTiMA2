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

test_that("missing native API reports isolated at-your-own-risk setup code", {
  local_mocked_bindings(
    .ctmaCtsemJuliaLaplaceStatus = function() structure(list(
      available = FALSE,
      reason = "The installed ctsem 3.11.1 is missing the native API."
    ), class = "ctmaCtsemJuliaLaplaceStatus"),
    .package = "CoTiMA2"
  )

  expect_error(
    .ctmaRequireCtsemJuliaLaplace(),
    "At your own risk"
  )
  install_code <- .ctmaCtsemJuliaInstallCode()
  expect_match(install_code, "cdriveraus/ctsem", fixed = TRUE)
  expect_match(install_code, "juliaFit", fixed = TRUE)
  expect_match(install_code, "ctsem-juliaFit-lib", fixed = TRUE)
  expect_match(install_code, "Restart R", fixed = TRUE)
})

test_that("ctsemjulia requests are checked centrally for every workflow", {
  calls <- 0L
  local_mocked_bindings(
    .ctmaRequireCtsemJuliaLaplace = function() {
      calls <<- calls + 1L
      list(available = TRUE)
    },
    .package = "CoTiMA2"
  )

  expect_silent(.ctmaCheckCtsemJuliaRequest("stan", "ctmaInit"))
  expect_identical(calls, 0L)
  expect_silent(.ctmaCheckCtsemJuliaRequest(
    "ctsemjulia", "ctmaFit", supported = TRUE
  ))
  expect_identical(calls, 1L)
  expect_silent(.ctmaCheckCtsemJuliaRequest(
    "ctsemjulia", "ctmaInit", supported = TRUE
  ))
  expect_identical(calls, 2L)
})

test_that("all backend entry points intercept unsupported ctsemjulia requests", {
  local_mocked_bindings(
    .ctmaRequireCtsemJuliaLaplace = function() list(available = TRUE),
    .package = "CoTiMA2"
  )
  entry_points <- list(
    ctmaAllInvFit = ctmaAllInvFit,
    ctmaEqual = ctmaEqual,
    ctmaPower = ctmaPower,
    ctmaOptimizeFit = ctmaOptimizeFit,
    ctmaReFit = ctmaReFit,
    ctmaFitConfig = ctmaFitConfig
  )
  for (name in names(entry_points)) {
    expect_error(
      do.call(entry_points[[name]], list(backend = "ctsemjulia")),
      paste0(name, "\\(\\) does not yet support")
    )
  }
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

test_that("native Laplace preparation accepts subject-only random effects", {
  model <- list(
    subjectIDname = "id",
    groupIDnames = character(),
    TIpredNames = character(),
    pars = data.frame(
      matrix = c("CINT", "DRIFT", "DIFFUSION"),
      param = c("cint_eta1", "drift_eta1", "diff_eta1"),
      indvarying = c(TRUE, FALSE, FALSE),
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
    .ctmaLaplaceStudyMembership = function(...) stop(
      "Study membership should not be constructed"
    ),
    .ctma_generate_ctsem_setup = function(...) list(),
    .package = "CoTiMA2"
  )

  prepared <- .ctmaPrepareCtsemJuliaLaplace(setup, verbose = FALSE)

  expect_length(prepared$selected_parameters, 0L)
  expect_identical(prepared$selected_subject_parameters, "cint_eta1")
  expect_false(any(prepared$model$pars$indvarying_study))
  expect_true(prepared$model$pars$indvarying[[1L]])
  expect_length(prepared$model$groupIDnames, 0L)
  expect_false("study" %in% names(prepared$data))
  expect_length(prepared$study_labels, 0L)
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

  backend_args <- ctmaBackendDefaults("ctsemjulia")
  backend_args$study_random_effects <- "CINT"
  result <- .ctmaFitCtsemJuliaBackend(
    setup = setup,
    backend_args = backend_args,
    cores = 2L,
    cotima_call = quote(ctmaFit(backend = "ctsemjulia"))
  )

  expect_identical(result$backend, "ctsemjulia")
  expect_identical(result$source_object, setup)
  expect_identical(result$ctsemjulia_prepared, prepared)
  expect_identical(result$ctsemjulia_mode, "laplace")
  expect_identical(result$argumentList$backend, "ctsemjulia")
  expect_true(result$argumentList$fit)
})

test_that("ctsem 3.12 optimizer diagnostics are normalized", {
  prepared <- structure(list(
    selected_parameters = "V1toV1",
    selected_subject_parameters = character()
  ), class = "ctmaCtsemJuliaLaplacePrep")
  result <- list(
    optim = list(
      converged = TRUE, gradient_norm = 2.3e-5, iterations = 24L,
      trace = data.frame(
        iteration = 0:2,
        inner_converged = c(25L, 25L, 25L)
      )
    ),
    estimate = list(subject_loglik = rep(-2.802063, 100L)),
    laplace = list(inner_converged = TRUE)
  )

  expect_message(
    annotated <- .ctmaAnnotateCtsemJuliaResult(
      result, prepared = prepared, use_laplace = TRUE
    ),
    "not a valid personwise decomposition"
  )
  diagnostic <- annotated$ctsemjulia_diagnostics
  expect_true(diagnostic$converged)
  expect_equal(diagnostic$gradient_norm, 2.3e-5)
  expect_identical(diagnostic$iterations, 24L)
  expect_identical(diagnostic$inner_modes$expected, 25L)
  expect_identical(diagnostic$inner_modes$converged, 25L)
  expect_true(diagnostic$inner_modes$all_converged)
  expect_false(diagnostic$casewise_loglik$usable_for_loo)
  expect_identical(
    diagnostic$casewise_loglik$scope, "study-marginal, not personwise"
  )
  expect_length(annotated$estimate$subject_loglik, 100L)
})

test_that("nonconstant personwise likelihood contributions remain usable", {
  result <- list(
    optim = list(converged = TRUE),
    estimate = list(subject_loglik = c(-1, -2, -3))
  )
  annotated <- .ctmaAnnotateCtsemJuliaResult(
    result, prepared = NULL, use_laplace = FALSE
  )
  expect_true(
    annotated$ctsemjulia_diagnostics$casewise_loglik$usable_for_loo
  )
})

test_that("an explicit failed native inner-mode flag is respected", {
  result <- list(
    optim = list(trace = data.frame(inner_converged = c(24L, 24L))),
    estimate = list(),
    laplace = list(inner_converged = FALSE)
  )
  annotated <- .ctmaAnnotateCtsemJuliaResult(
    result, prepared = NULL, use_laplace = TRUE
  )
  expect_false(annotated$ctsemjulia_diagnostics$inner_modes$all_converged)
  expect_false(
    annotated$ctsemjulia_diagnostics$casewise_loglik$usable_for_loo
  )
})

test_that("ctsemjulia uses direct native ML when no random effects exist", {
  setup <- structure(list(
    data = data.frame(id = rep(1:2, each = 2), y = 1:4),
    ctModel = list(
      pars = data.frame(
        param = "drift_eta1", indvarying = FALSE,
        stringsAsFactors = FALSE
      )
    ),
    argumentList = list(fit = FALSE)
  ), class = "CoTiMAFit")

  local_mocked_bindings(
    .ctmaFitCtsemJuliaDirect = function(...) list(
      estimate = list(loglik = -9, converged = TRUE)
    ),
    .ctmaPrepareCtsemJuliaLaplace = function(...) {
      stop("Laplace preparation must not run")
    },
    .package = "CoTiMA2"
  )

  result <- .ctmaFitCtsemJuliaBackend(
    setup = setup,
    backend_args = ctmaBackendDefaults("ctsemjulia"),
    cores = 2L,
    cotima_call = quote(ctmaFit(backend = "ctsemjulia"))
  )

  expect_identical(result$backend, "ctsemjulia")
  expect_identical(result$ctsemjulia_mode, "direct")
  expect_null(result$ctsemjulia_prepared)
  expect_identical(result$estimate$loglik, -9)
})

test_that("ctsemjulia detects explicit and existing random effects", {
  setup <- structure(list(ctModel = list(pars = data.frame(
    param = "drift_eta1", indvarying = FALSE,
    stringsAsFactors = FALSE
  ))), class = "CoTiMAFit")
  defaults <- ctmaBackendDefaults("ctsemjulia")

  expect_false(.ctmaCtsemJuliaRandomEffectsRequested(setup, defaults))
  defaults$study_parameters <- "drift_eta1"
  expect_true(.ctmaCtsemJuliaRandomEffectsRequested(setup, defaults))
  defaults$study_parameters <- NULL
  defaults$reDRIFT <- matrix(FALSE, 1L, 1L)
  expect_false(.ctmaCtsemJuliaRandomEffectsRequested(setup, defaults))
  defaults$reDRIFT[1L, 1L] <- TRUE
  expect_true(.ctmaCtsemJuliaRandomEffectsRequested(setup, defaults))

  setup$ctModel$pars$indvarying <- TRUE
  defaults$reDRIFT <- NULL
  expect_true(.ctmaCtsemJuliaRandomEffectsRequested(setup, defaults))
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
