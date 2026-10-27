test_that("Stan control overrides are merged by name and recursively", {
  defaults <- list(
    iter = 1000,
    verbose = 0,
    optimcontrol = list(is = FALSE, finishsamples = 1000)
  )
  overrides <- list(
    verbose = 2,
    iter = 500,
    optimcontrol = list(is = NULL)
  )

  result <- CoTiMA2:::ctmaMergeStanctArgs(defaults, overrides)

  expect_identical(result$iter, 500)
  expect_identical(result$verbose, 2)
  expect_true("is" %in% names(result$optimcontrol))
  expect_null(result$optimcontrol$is)
  expect_identical(result$optimcontrol$finishsamples, 1000)
})

test_that("default importance-sampling controls are NULL", {
  expect_null(CoTiMA2::CoTiMAStanctArgs$optimcontrol$is)
  expect_null(CoTiMA2::CoTiMAStanctArgs$optimcontrol$isitersize)
  expect_null(CoTiMA2::CoTiMAStanctArgs$optimcontrol$isESS)
})

test_that("Stan control overrides reject unnamed and unknown values", {
  defaults <- list(iter = 1000)

  expect_error(
    CoTiMA2:::ctmaMergeStanctArgs(defaults, list(500)),
    "fully named list"
  )
  expect_error(
    CoTiMA2:::ctmaMergeStanctArgs(defaults, list(unknown = 1)),
    "Unknown CoTiMAStanctArgs: unknown"
  )
})

test_that("removed ctsem variational-Bayes controls are handled explicitly", {
  sanitized <- CoTiMA2:::.ctmaSanitizeCtFitArgs(list(
    optimize = TRUE,
    vb = FALSE
  ))
  expect_identical(sanitized, list(optimize = TRUE))
  expect_error(
    CoTiMA2:::.ctmaSanitizeCtFitArgs(list(vb = TRUE)),
    "Variational Bayes is no longer supported by ctsem"
  )
})

test_that("ctmaInit resolves legacy Stan controls recursively by name", {
  result <- CoTiMA2:::.ctmaResolveInitStanArgs(
    overrides = list(
      verbose = 2,
      iter = 750,
      optimcontrol = list(finishsamples = 250)
    )
  )

  expect_identical(result$verbose, 2)
  expect_identical(result$iter, 750)
  expect_identical(result$optimcontrol$finishsamples, 250)
  expect_identical(
    result$optimcontrol$carefulfit,
    CoTiMAStanctArgs$optimcontrol$carefulfit
  )
})

test_that("explicit ctmaInit arguments override the legacy Stan list", {
  result <- CoTiMA2:::.ctmaResolveInitStanArgs(
    overrides = list(
      scaleTI = FALSE,
      iter = 400,
      verbose = 1,
      optimcontrol = list(finishsamples = 100)
    ),
    scaleTI = TRUE,
    iter = 900,
    verbose = 3,
    finishsamples = 500
  )

  expect_true(result$scaleTI)
  expect_identical(result$iter, 900)
  expect_identical(result$verbose, 3)
  expect_identical(result$optimcontrol$finishsamples, 500)
})

test_that("ctmaInit backend merge detects misordered and unknown controls", {
  result <- CoTiMA2:::.ctmaResolveInitStanArgs(list(
    verbose = 4,
    iter = 321
  ))
  expect_identical(result$verbose, 4)
  expect_identical(result$iter, 321)

  expect_error(
    CoTiMA2:::.ctmaResolveInitStanArgs(list(not_a_control = TRUE)),
    "Unknown CoTiMAStanctArgs: not_a_control"
  )
  expect_error(
    CoTiMA2:::.ctmaResolveInitStanArgs(list(
      optimcontrol = list(not_a_control = TRUE)
    )),
    "Unknown CoTiMAStanctArgs\\$optimcontrol: not_a_control"
  )
})

test_that("ctmaFit resolves partial legacy Stan controls recursively", {
  result <- CoTiMA2:::.ctmaResolveFitStanArgs(
    overrides = list(
      verbose = 2,
      control = list(adapt_delta = 0.95),
      optimcontrol = list(finishsamples = 300)
    )
  )

  expect_identical(result$verbose, 2)
  expect_identical(result$control$adapt_delta, 0.95)
  expect_identical(
    result$control$max_treedepth,
    CoTiMAStanctArgs$control$max_treedepth
  )
  expect_identical(result$optimcontrol$finishsamples, 300)
  expect_identical(
    result$optimcontrol$carefulfit,
    CoTiMAStanctArgs$optimcontrol$carefulfit
  )
})

test_that("ctmaFit direct controls have defined precedence", {
  inherited <- CoTiMA2:::.ctmaResolveFitStanArgs(
    overrides = list(verbose = 2, iter = 400, priors = TRUE),
    optimize = TRUE,
    priors = FALSE,
    iter = 900,
    verbose_was_supplied = FALSE
  )
  expect_identical(inherited$verbose, 2)
  expect_identical(inherited$iter, 900)
  expect_false(inherited$priors)

  explicit <- CoTiMA2:::.ctmaResolveFitStanArgs(
    overrides = list(verbose = 2),
    verbose = 0,
    verbose_was_supplied = TRUE
  )
  expect_identical(explicit$verbose, 0)
})

test_that("Stan warmup defaults to half of iter and can be overridden", {
  defaults <- CoTiMA2:::.ctmaResolveFitStanArgs(iter = 100L)
  expect_null(defaults$warmup)
  expect_identical(CoTiMA2:::.ctmaStanWarmup(defaults), 50L)

  explicit <- CoTiMA2:::.ctmaResolveFitStanArgs(
    overrides = list(warmup = 30L),
    iter = 100L
  )
  expect_identical(CoTiMA2:::.ctmaStanWarmup(explicit), 30L)
  expect_error(
    CoTiMA2:::.ctmaStanWarmup(list(iter = 100L, warmup = 100L)),
    "0 <= warmup < iter"
  )
})

test_that("ctmaFit rejects unknown legacy Stan controls", {
  expect_error(
    CoTiMA2:::.ctmaResolveFitStanArgs(list(not_a_control = TRUE)),
    "Unknown CoTiMAStanctArgs: not_a_control"
  )
  expect_error(
    CoTiMA2:::.ctmaResolveFitStanArgs(list(
      control = list(not_a_control = TRUE)
    )),
    "Unknown CoTiMAStanctArgs\\$control: not_a_control"
  )
})

test_that("ctmaFit applies backend validation before fitting", {
  expect_error(
    ctmaFit(
      ctmaInitFit = NULL,
      CoTiMAStanctArgs = list(control = list(adapt_delta = 0.9))
    ),
    "A fitted CoTiMA object has to be supplied"
  )
  expect_error(
    ctmaFit(
      ctmaInitFit = NULL,
      CoTiMAStanctArgs = list(not_a_control = TRUE)
    ),
    "Unknown CoTiMAStanctArgs: not_a_control"
  )
})

test_that("ctmaFit keeps Stan as its default backend", {
  expect_identical(
    eval(formals(ctmaFit)$backend),
    c("stan", "julia", "ctsemjulia")
  )

  resolved <- CoTiMA2:::.ctmaResolveFitBackend()
  expect_identical(resolved$backend, "stan")
  expect_null(resolved$stan_args)
})

test_that("ctmaFit accepts either new or legacy Stan controls", {
  modern <- CoTiMA2:::.ctmaResolveFitBackend(
    backendArgs = list(control = list(adapt_delta = 0.9))
  )
  legacy <- CoTiMA2:::.ctmaResolveFitBackend(
    legacy_stan_args = list(control = list(adapt_delta = 0.8))
  )
  expect_identical(modern$stan_args$control$adapt_delta, 0.9)
  expect_identical(legacy$stan_args$control$adapt_delta, 0.8)

  expect_error(
    CoTiMA2:::.ctmaResolveFitBackend(
      backendArgs = list(iter = 500),
      legacy_stan_args = list(iter = 1000)
    ),
    "either 'CoTiMAStanctArgs' or 'backendArgs'"
  )
  expect_error(
    CoTiMA2:::.ctmaResolveFitBackend(
      backendArgs = list(not_a_control = TRUE)
    ),
    "Unknown backendArgs: not_a_control"
  )
})

test_that("ctmaFit resolves Julia ML and Bayesian controls", {
  ml <- CoTiMA2:::.ctmaResolveFitBackend("julia", estimation = "ML")
  bayes <- CoTiMA2:::.ctmaResolveFitBackend(
    "julia",
    backendArgs = list(chains = 2L, metric = "dense"),
    estimation = "Bayes"
  )

  expect_identical(ml$backend, "julia")
  expect_identical(ml$estimation, "ml")
  expect_true(ml$backend_args$uncertainty)
  expect_true(ml$backend_args$transformed)
  expect_identical(bayes$estimation, "bayes")
  expect_identical(bayes$backend_args$chains, 2L)
  expect_identical(bayes$backend_args$metric, "dense")
  expect_error(
    CoTiMA2:::.ctmaResolveFitBackend(
      "julia", legacy_stan_args = list(iter = 1000)
    ),
    "cannot be used with backend = 'julia'"
  )
})

test_that("Julia backend requests point users to the installation guide", {
  expect_message(
    resolved <- CoTiMA2:::.ctmaResolveFitBackend(
      backend = "julia", estimation = "ML", workflow = "ctmaFit"
    ),
    'vignette\\("Julia_Install", package = "CoTiMA2"\\)'
  )
  expect_identical(resolved$backend, "julia")
  expect_silent(
    CoTiMA2:::.ctmaResolveFitBackend(
      backend = "stan", estimation = "ML", workflow = "ctmaFit"
    )
  )
})

test_that("ctmaInit exposes Julia backend defaults and rejects outer parallelism", {
  expect_identical(
    eval(formals(ctmaInit)$backend),
    c("stan", "julia", "ctsemjulia")
  )
  defaults <- ctmaBackendDefaults("julia", "ML", "ctmaInit")
  expect_true(defaults$uncertainty)
  expect_true(defaults$transformed)
  expect_error(
    ctmaInit(
      primaryStudies = NULL,
      backend = "julia",
      doPar = 2L
    ),
    "currently requires doPar = 1"
  )
})

test_that("ctmaInit native ctsem Julia dispatcher preserves ctFit compatibility", {
  captured <- NULL
  testthat::local_mocked_bindings(
    .ctmaRequireCtsemJuliaLaplace = function() list(available = TRUE),
    .ctmaNativeCtsemJuliaFit = function(args) {
      captured <<- args
      list(
        pars = data.frame(),
        standata = list(priors = 0L)
      )
    },
    .ctmaNativeCtsemJuliaSummary = function(fit, digits) {
      list(loglik = -4, npars = 2L, digits = digits)
    },
    .package = "CoTiMA2"
  )

  model <- list(pars = data.frame(
    param = "drift_eta1", indvarying = FALSE
  ))
  result <- CoTiMA2:::.ctmaInitCtsemJuliaStudyFit(
    datalong = data.frame(id = 1L, time = 0, V1 = 0),
    model = model,
    same_initial_times = FALSE,
    backend_args = list(
      inits = NULL, optimcontrol = NULL, priors = NULL, verbose = NULL
    ),
    start_values = 0.2,
    priors = FALSE,
    finishsamples = 25L,
    cores = 3L,
    verbose = 0,
    digits = 4L
  )

  expect_identical(captured$backend, "julia")
  expect_identical(captured$intoverpop, "auto")
  expect_identical(captured$cores, 3L)
  expect_identical(captured$inits, 0.2)
  expect_identical(captured$optimcontrol$finishsamples, 25L)
  expect_identical(result$backend, "ctsemjulia")
  expect_identical(result$ctstanmodelbase, model)
  expect_identical(result$ctsemjulia_mode, "direct")
  expect_identical(result$resultsSummary$npars, 2L)
})

test_that("ctmaInit converts native Julia random-effect summaries", {
  popsd <- data.frame(mean = c(2, 3, 4), sd = NA_real_, check.names = FALSE)
  rownames(popsd) <- c("a", "b", "c")
  popcor <- data.frame(mean = c(.1, .2, .3), sd = NA_real_)
  rownames(popcor) <- c("b__a", "c__a", "c__b")

  result <- CoTiMA2:::.ctmaInitCtsemJuliaRandomEffectSummary(
    list(popsd = popsd, rawpopcorr = popcor)
  )

  expect_equal(unname(diag(result$popcov_mean)), c(4, 9, 16))
  expect_equal(
    result$popcor_mean[lower.tri(result$popcor_mean)],
    c(.1, .2, .3)
  )
  expect_equal(result$popcov_mean["b", "a"], .6)
  expect_true(all(is.na(result$popcov_sd)))
})

test_that("ctModelText reconstructs the exact generated ctsem model", {
  model <- ctsem::ctModel(
    type = "ct",
    n.latent = 2L,
    n.manifest = 2L,
    manifestNames = c("binary", "continuous"),
    latentNames = c("eta1", "eta2"),
    manifesttype = c(1L, 0L),
    LAMBDA = diag(2),
    DRIFT = matrix(c("a11", "a21", "a12", "a22"), 2L, 2L),
    DIFFUSION = matrix(c("d1", "d21", 0, "d2"), 2L, 2L),
    T0VAR = matrix(c("t1", "t21", 0, "t2"), 2L, 2L),
    CINT = matrix(c("c1", "c2"), 2L, 1L),
    n.TIpred = 1L,
    TIpredNames = "study"
  )
  model$pars[model$pars$matrix == "CINT", "indvarying"] <- TRUE
  model$pars[model$pars$matrix == "DRIFT", "study_effect"] <- c(
    TRUE, FALSE, FALSE, TRUE
  )

  code <- CoTiMA2:::.ctmaModelText(model, object_name = "recreated")
  expect_type(code, "character")
  expect_length(code, 1L)
  expect_match(code, "ctsem::ctModel")

  target <- new.env(parent = globalenv())
  eval(parse(text = code), envir = target)
  recreated <- target$recreated

  expect_equal(recreated$matrices, model$matrices)
  expect_equal(recreated$pars, model$pars)
  expect_equal(recreated$manifesttype, model$manifesttype)
  expect_identical(recreated$TIpredNames, model$TIpredNames)
})

test_that("ctmaFit Julia dispatcher forwards the prepared setup", {
  testthat::local_mocked_bindings(
    ctmaJuliaFit = function(object, estimation, ...) {
      list(object = object, estimation = estimation, controls = list(...))
    },
    .package = "CoTiMA2"
  )
  setup <- structure(list(marker = "prepared"), class = "CoTiMAFit")
  result <- CoTiMA2:::.ctmaJuliaFitFromSetup(
    setup = setup,
    estimation = "ML",
    backend_args = list(maxiters = 321L, julia_threads = 2L)
  )

  expect_identical(result$object, setup)
  expect_identical(result$estimation, "ML")
  expect_identical(result$controls$maxiters, 321L)
  expect_identical(result$controls$julia_threads, 2L)
})

test_that("ctmaFit suppresses the no-fit banner for internal Julia setup", {
  expect_true(CoTiMA2:::.ctmaShouldShowNoFitBanner(FALSE, FALSE))
  expect_false(CoTiMA2:::.ctmaShouldShowNoFitBanner(FALSE, TRUE))
  expect_false(CoTiMA2:::.ctmaShouldShowNoFitBanner(TRUE, FALSE))
})

test_that("Julia resource note reports effective Julia and BLAS threads", {
  expect_message(
    CoTiMA2:::.ctmaJuliaThreadNote(8L, 1L),
    "8 Julia pool thread\\(s\\).*8 concurrent subject task.*1 BLAS"
  )
})

test_that("ctmaFit-style automatic Julia resources use all resolved cores", {
  configured <- CoTiMA2:::.ctmaInitJuliaResources(
    cores = 9L,
    estimation = "ML",
    backendArgs = list(),
    resolved_backend_args = ctmaBackendDefaults("julia", "ML", "ctmaFit"),
    announce = FALSE,
    workflow = "ctmaFit"
  )
  expect_true(configured$automatic)
  expected_pool <- CoTiMA2:::.ctmaResolveJuliaThreads("auto")
  expect_identical(configured$backend_args$julia_threads, expected_pool)
  expect_identical(configured$backend_args$subject_threads, expected_pool)
  expect_identical(configured$backend_args$blas_threads, 1L)

  explicit <- CoTiMA2:::.ctmaInitJuliaResources(
    cores = 9L,
    estimation = "ML",
    backendArgs = list(julia_threads = 3L),
    resolved_backend_args = utils::modifyList(
      ctmaBackendDefaults("julia", "ML", "ctmaFit"),
      list(julia_threads = 3L)
    ),
    announce = FALSE,
    workflow = "ctmaFit"
  )
  expect_false(explicit$automatic)
  expect_identical(explicit$backend_args$julia_threads, 3L)
})

test_that("automatic Julia resource notes identify their calling workflow", {
  expect_message(
    CoTiMA2:::.ctmaInitJuliaResources(
      cores = 6L,
      estimation = "ML",
      backendArgs = list(),
      resolved_backend_args = ctmaBackendDefaults("julia", "ML", "ctmaEqual"),
      workflow = "ctmaEqual"
    ),
    "ctmaEqual automatically selected.*Julia pool thread.*subject task.*1 BLAS"
  )
})

test_that("ctmaInit Julia study dispatcher restores ctsem compatibility", {
  transformed_summary <- list(loglik = -12, npars = 2L)
  testthat::local_mocked_bindings(
    ctmaJuliaFit = function(object, estimation, ...) {
      structure(
        list(source_object = object, controls = list(...)),
        class = "ctmaJuliaFit"
      )
    },
    ctAddJuliaSamples = function(fit, setup, nsamples, cores, seed, digits) {
      fit$ctsem_summary <- transformed_summary
      fit
    },
    .ctmaJuliaGenerateSetup = function(object) {
      list(standata = list(priors = 0L))
    },
    .package = "CoTiMA2"
  )

  result <- CoTiMA2:::.ctmaInitJuliaStudyFit(
    datalong = data.frame(id = 1, time = 0, V1 = 0),
    model = list(marker = "model"),
    same_initial_times = FALSE,
    stan_args = list(iter = 100L),
    backend_args = list(
      transformed_draws = 20L,
      transformed_cores = 1L,
      seed = 3L,
      inits = NULL
    ),
    estimation = "ML",
    digits = 4L
  )

  expect_s3_class(result, "ctmaJuliaFit")
  expect_identical(result$resultsSummary, transformed_summary)
  expect_identical(result$standata$priors, 0L)
  expect_identical(result$ctstanmodelbase$marker, "model")
  expect_true(result$controls$uncertainty)
  expect_false(result$controls$transformed)
  expect_identical(result$backend, "julia")
})

test_that("ctmaFit backendArgs work through the public entry point", {
  expect_error(
    ctmaFit(
      ctmaInitFit = NULL,
      backend = "stan",
      backendArgs = list(control = list(adapt_delta = 0.9))
    ),
    "A fitted CoTiMA object has to be supplied"
  )
  expect_error(
    ctmaFit(
      ctmaInitFit = NULL,
      CoTiMAStanctArgs = list(iter = 1000),
      backendArgs = list(iter = 500)
    ),
    "either 'CoTiMAStanctArgs' or 'backendArgs'"
  )
  expect_error(
    ctmaFit(ctmaInitFit = NULL, backend = "julia"),
    "A fitted CoTiMA object has to be supplied"
  )
})

test_that("ctmaAllInvFit resolves partial nested Stan controls", {
  result <- CoTiMA2:::.ctmaResolveAllInvStanArgs(
    overrides = list(
      verbose = 2,
      iter = 400,
      control = list(adapt_delta = 0.95),
      optimcontrol = list(finishsamples = 250)
    ),
    scaleTime = 0.5,
    iter = 900,
    finishsamples = 500
  )

  expect_identical(result$verbose, 2)
  expect_identical(result$scaleTime, 0.5)
  expect_identical(result$iter, 900)
  expect_identical(result$control$adapt_delta, 0.95)
  expect_identical(
    result$control$max_treedepth,
    CoTiMAStanctArgs$control$max_treedepth
  )
  expect_identical(result$optimcontrol$finishsamples, 500)
  expect_identical(
    result$optimcontrol$carefulfit,
    CoTiMAStanctArgs$optimcontrol$carefulfit
  )
})

test_that("ctmaAllInvFit explicit verbose overrides the legacy list", {
  inherited <- CoTiMA2:::.ctmaResolveAllInvStanArgs(
    overrides = list(verbose = 2),
    verbose_was_supplied = FALSE
  )
  explicit <- CoTiMA2:::.ctmaResolveAllInvStanArgs(
    overrides = list(verbose = 2),
    verbose = 0,
    verbose_was_supplied = TRUE
  )
  expect_identical(inherited$verbose, 2)
  expect_identical(explicit$verbose, 0)
})

test_that("ctmaAllInvFit validates controls before requiring a fit", {
  expect_error(
    ctmaAllInvFit(
      ctmaInitFit = NULL,
      CoTiMAStanctArgs = list(control = list(adapt_delta = 0.9))
    ),
    "A fitted CoTiMA object has to be supplied"
  )
  expect_error(
    ctmaAllInvFit(
      ctmaInitFit = NULL,
      CoTiMAStanctArgs = list(not_a_control = TRUE)
    ),
    "Unknown CoTiMAStanctArgs: not_a_control"
  )
})

test_that("ctmaPower resolves partial nested Stan controls", {
  result <- CoTiMA2:::.ctmaResolvePowerStanArgs(
    overrides = list(
      scaleTime = 0.5,
      verbose = 2,
      control = list(adapt_delta = 0.95),
      optimcontrol = list(finishsamples = 250)
    ),
    iter = 900,
    finishsamples = 500
  )

  expect_identical(result$scaleTime, 0.5)
  expect_identical(result$verbose, 2)
  expect_identical(result$iter, 900)
  expect_identical(result$control$adapt_delta, 0.95)
  expect_identical(
    result$control$max_treedepth,
    CoTiMAStanctArgs$control$max_treedepth
  )
  expect_identical(result$optimcontrol$finishsamples, 500)
  expect_identical(
    result$optimcontrol$carefulfit,
    CoTiMAStanctArgs$optimcontrol$carefulfit
  )
})

test_that("ctmaPower direct scale and verbose controls take precedence", {
  result <- CoTiMA2:::.ctmaResolvePowerStanArgs(
    overrides = list(scaleTime = 0.5, verbose = 2),
    scaleTime = 2,
    verbose = 0,
    verbose_was_supplied = TRUE
  )
  expect_identical(result$scaleTime, 2)
  expect_identical(result$verbose, 0)
})

test_that("ctmaPower validates controls before requiring a fit", {
  expect_error(
    ctmaPower(
      ctmaInitFit = NULL,
      statisticalPower = 0.8,
      CoTiMAStanctArgs = list(control = list(adapt_delta = 0.9))
    ),
    "A fitted CoTiMA object.*has to be supplied"
  )
  expect_error(
    ctmaPower(
      ctmaInitFit = NULL,
      statisticalPower = 0.8,
      CoTiMAStanctArgs = list(not_a_control = TRUE)
    ),
    "Unknown CoTiMAStanctArgs: not_a_control"
  )
})

test_that("ctmaPower exposes and validates backend controls", {
  expect_identical(eval(formals(ctmaPower)$backend), c("stan", "julia"))
  expect_true("backendArgs" %in% names(formals(ctmaPower)))

  expect_error(
    ctmaPower(
      ctmaInitFit = NULL,
      statisticalPower = 0.8,
      backend = "julia",
      backendArgs = list(not_a_control = TRUE)
    ),
    "Unknown backendArgs: not_a_control"
  )
  expect_error(
    ctmaPower(
      ctmaInitFit = NULL,
      statisticalPower = 0.8,
      backend = "julia",
      CoTiMAStanctArgs = list()
    ),
    "cannot be used with backend = 'julia'"
  )
})

test_that("ctmaEqual exposes and validates backend controls", {
  expect_identical(eval(formals(ctmaEqual)$backend), c("stan", "julia"))
  expect_true("backendArgs" %in% names(formals(ctmaEqual)))
  expect_true("ctmaEqual" %in% eval(formals(ctmaBackendDefaults)$workflow))

  fixture <- list(
    model.type = "ct",
    argumentList = list(CoTiMAStanctArgs = NULL, estimation = "ML"),
    modelResults = list(DRIFT = c(`a (invariant)` = 0, `b (invariant)` = 0))
  )
  expect_error(
    ctmaEqual(
      ctmaInvariantFit = fixture,
      backend = "julia",
      backendArgs = list(not_a_control = TRUE)
    ),
    "Unknown backendArgs: not_a_control"
  )
})

test_that("ctmaOptimizeFit merges defaults, inherited controls, and overrides", {
  result <- CoTiMA2:::.ctmaResolveOptimizeStanArgs(
    inherited = list(
      verbose = 2,
      iter = 700,
      control = list(adapt_delta = 0.9),
      optimcontrol = list(finishsamples = 300)
    ),
    overrides = list(
      iter = 900,
      control = list(max_treedepth = 12L)
    ),
    finishsamples = 500
  )

  expect_identical(result$verbose, 2)
  expect_identical(result$iter, 900)
  expect_identical(result$control$adapt_delta, 0.9)
  expect_identical(result$control$max_treedepth, 12L)
  expect_identical(result$optimcontrol$finishsamples, 500)
  expect_identical(
    result$optimcontrol$carefulfit,
    CoTiMAStanctArgs$optimcontrol$carefulfit
  )
})

test_that("ctmaOptimizeFit rejects unknown stored and supplied controls", {
  expect_error(
    CoTiMA2:::.ctmaResolveOptimizeStanArgs(
      inherited = list(not_a_control = TRUE)
    ),
    "Unknown stored CoTiMAStanctArgs: not_a_control"
  )
  expect_error(
    CoTiMA2:::.ctmaResolveOptimizeStanArgs(
      overrides = list(not_a_control = TRUE)
    ),
    "Unknown CoTiMAStanctArgs: not_a_control"
  )
})

test_that("ctmaOptimizeFit recognizes only finite fit criteria", {
  valid <- list(summary = list(minus2ll = 123.4))
  nonfinite <- list(summary = list(minus2ll = Inf))
  missing <- list(summary = list())
  failed <- simpleError("failed")

  expect_identical(CoTiMA2:::.ctmaOptimizeFitValue(valid), 123.4)
  expect_true(is.na(CoTiMA2:::.ctmaOptimizeFitValue(nonfinite)))
  expect_true(is.na(CoTiMA2:::.ctmaOptimizeFitValue(missing)))
  expect_true(is.na(CoTiMA2:::.ctmaOptimizeFitValue(failed)))
})

test_that("stored Julia routing fixture has a CoTiMA-compatible criterion", {
  fixture_path <- testthat::test_path("fixtures", "julia-routing-fixture.R")
  fixture_environment <- new.env(parent = baseenv())
  sys.source(fixture_path, envir = fixture_environment)
  fixture <- fixture_environment$ctma_julia_routing_fixture()

  result <- CoTiMA2:::.ctmaFinalizeJuliaCoTiMAFit(
    fit = fixture$fit,
    setup = fixture$setup,
    cotima_call = quote(ctmaFit(backend = "julia")),
    backend_args = list(maxiters = 25L)
  )

  expect_s3_class(result, "ctmaJuliaFit")
  expect_identical(result$summary$minus2ll, 247)
  expect_identical(result$summary$loglik, -123.5)
  expect_identical(result$summary$n.parameters, 3L)
  expect_identical(result$argumentList$backend, "julia")
  expect_identical(result$argumentList$backendArgs$maxiters, 25L)
  expect_true(result$argumentList$fit)
  expect_identical(CoTiMA2:::.ctmaOptimizeFitValue(result), 247)
  expect_identical(result$source_object, fixture$setup)
})

test_that("ctmaOptimizeFit validates backend controls before input objects", {
  expect_identical(
    eval(formals(ctmaOptimizeFit)$backend),
    c("stan", "julia")
  )
  expect_error(
    ctmaOptimizeFit(
      CoTiMAStanctArgs = list(control = list(adapt_delta = 0.9))
    ),
    "ctmaFitFit or ctmaIniFit"
  )
  expect_error(
    ctmaOptimizeFit(CoTiMAStanctArgs = list(not_a_control = TRUE)),
    "Unknown CoTiMAStanctArgs: not_a_control"
  )
  expect_error(
    ctmaOptimizeFit(backend = "julia", reFits = 2L),
    "requires the complete-model 'ctmaFitFit' branch"
  )
  expect_error(
    ctmaOptimizeFit(
      backend = "stan",
      backendArgs = list(not_a_control = TRUE)
    ),
    "Unknown backendArgs: not_a_control"
  )
})

test_that("parallel fit selection uses the greatest finite log-posterior", {
  attempt <- function(value) {
    list(
      ok = TRUE,
      fit = list(stanfit = list(optimfit = list(value = value))),
      warn_hessinv = FALSE,
      warnings = character(),
      error = NA_character_
    )
  }

  result <- CoTiMA2:::ctmaSelectBestParallelFit(list(
    attempt(-100),
    attempt(-10),
    attempt(NA_real_)
  ))

  expect_identical(result$logposterior, -10)
  expect_identical(result$fit$stanfit$optimfit$value, -10)
})

test_that("parallel fit selection ignores errors and reports total failure", {
  failed <- list(
    ok = FALSE,
    fit = NULL,
    warn_hessinv = FALSE,
    warnings = character(),
    error = "fit failed"
  )
  successful <- list(
    ok = TRUE,
    fit = list(stanfit = list(optimfit = list(value = -5))),
    warn_hessinv = FALSE,
    warnings = character(),
    error = NA_character_
  )

  expect_identical(
    CoTiMA2:::ctmaSelectBestParallelFit(list(failed, successful))$logposterior,
    -5
  )
  expect_error(
    CoTiMA2:::ctmaSelectBestParallelFit(list(failed)),
    "All parallel fitting attempts failed"
  )
})

test_that("ctmaInit maps Julia individual variation to pseudo processes only", {
  defaults <- list(gradient_tolerance = 0.005)
  manifest <- CoTiMA2:::.ctmaInitJuliaIndividualVariation(
    "MANIFEST", FALSE, "julia", list(), defaults, announce = FALSE
  )
  expect_true(manifest$translated)
  expect_identical(manifest$indVarying, FALSE)
  expect_identical(manifest$randomIntercepts, "MANIFEST")
  expect_identical(manifest$backend_args$gradient_tolerance, 0.05)

  cint <- CoTiMA2:::.ctmaInitJuliaIndividualVariation(
    "CINT", FALSE, "julia", list(gradient_tolerance = 0.01),
    list(gradient_tolerance = 0.01), announce = FALSE
  )
  expect_identical(cint$randomIntercepts, "CINT")
  expect_identical(cint$backend_args$gradient_tolerance, 0.01)

  stan <- CoTiMA2:::.ctmaInitJuliaIndividualVariation(
    "MANIFEST", FALSE, "stan", list(), defaults, announce = FALSE
  )
  expect_false(stan$translated)
  expect_identical(stan$indVarying, "MANIFEST")
  expect_identical(stan$randomIntercepts, FALSE)
})

test_that("ctmaInit automatically allocates Julia resources from coresToUse", {
  ml <- CoTiMA2:::.ctmaInitJuliaResources(
    63L, "ML", list(),
    list(multistart = 1L, julia_threads = 1L), announce = FALSE
  )
  expect_true(ml$automatic)
  expect_identical(ml$backend_args$julia_threads, 63L)

  multistart <- CoTiMA2:::.ctmaInitJuliaResources(
    63L, "ML", list(),
    list(multistart = 4L, julia_threads = 1L), announce = FALSE
  )
  expect_identical(multistart$backend_args$parallel_starts, 4L)
  expect_identical(multistart$backend_args$threads_per_start, 15L)
  expect_identical(multistart$backend_args$julia_threads, 63L)

  bayes <- CoTiMA2:::.ctmaInitJuliaResources(
    63L, "Bayes", list(),
    list(chains = 4L, julia_threads = 1L), announce = FALSE
  )
  expect_identical(bayes$backend_args$parallel_chains, 4L)
  expect_identical(bayes$backend_args$threads_per_chain, 15L)
  expect_identical(bayes$backend_args$julia_threads, 1L)

  explicit <- CoTiMA2:::.ctmaInitJuliaResources(
    63L, "ML", list(threads_per_start = 3L),
    list(multistart = 4L, julia_threads = 1L,
         parallel_starts = 1L, threads_per_start = 3L),
    announce = FALSE
  )
  expect_false(explicit$automatic)
  expect_identical(explicit$backend_args$julia_threads, 1L)
  expect_identical(explicit$backend_args$parallel_starts, 1L)
})

test_that("ctmaInit maps finishsamples to Julia transformed draws", {
  defaults <- list(transformed_draws = 1000L, seed = NULL)

  mapped <- CoTiMA2:::.ctmaInitJuliaFinishsamples(
    finishsamples = 10000,
    backendArgs = list(),
    resolved_backend_args = defaults
  )
  expect_identical(mapped$transformed_draws, 10000L)

  explicit <- CoTiMA2:::.ctmaInitJuliaFinishsamples(
    finishsamples = 10000,
    backendArgs = list(transformed_draws = 25000L),
    resolved_backend_args = list(transformed_draws = 25000L, seed = NULL)
  )
  expect_identical(explicit$transformed_draws, 25000L)

  unchanged <- CoTiMA2:::.ctmaInitJuliaFinishsamples(
    finishsamples = NULL,
    backendArgs = list(),
    resolved_backend_args = defaults
  )
  expect_identical(unchanged$transformed_draws, 1000L)

  expect_error(
    CoTiMA2:::.ctmaInitJuliaFinishsamples(
      finishsamples = 1,
      backendArgs = list(),
      resolved_backend_args = defaults
    ),
    "integer of at least 2"
  )
})

test_that("ctmaInit warnings identify their primary study", {
  messages <- character()
  value <- withCallingHandlers(
    CoTiMA2:::.ctmaWithInitStudyWarnings(
      {
        warning("Hessian repaired", call. = FALSE)
        warning("wide uncertainty intervals", call. = FALSE)
        42L
      },
      index = 2L,
      total = 3L,
      study = list(originalStudyNo = "B7", source = "Smith et al. (2024)")
    ),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_identical(value, 42L)
  expect_length(messages, 2L)
  expect_match(
    messages[[1L]],
    "^\\[ctmaInit primary study 2 of 3; Study: B7; Source: Smith et al\\. \\(2024\\)\\] Hessian repaired$"
  )
  expect_match(messages[[2L]], "wide uncertainty intervals$")

  no_reference <- CoTiMA2:::.ctmaInitStudyWarningLabel(
    1L, 3L,
    list(originalStudyNo = 1L, source = "Reference not provided")
  )
  expect_identical(no_reference, "[ctmaInit primary study 1 of 3; Study: 1] ")
})
