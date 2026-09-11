# Internal helpers used by ctmaInit().

.ctmaInitJuliaIndividualVariation <- function(
    indVarying, randomIntercepts, backend, backendArgs,
    resolved_backend_args, announce = TRUE) {
  translated <- identical(backend, "julia") &&
    identical(randomIntercepts, FALSE) &&
    indVarying %in% c("MANIFEST", "CINT")
  if (!translated) {
    return(list(
      indVarying = indVarying,
      randomIntercepts = randomIntercepts,
      backend_args = resolved_backend_args,
      translated = FALSE
    ))
  }
  randomIntercepts <- indVarying
  indVarying <- FALSE
  if (!"gradient_tolerance" %in% names(backendArgs)) {
    resolved_backend_args$gradient_tolerance <- 0.05
  }
  if (isTRUE(announce)) {
    message(
      "Note: For ctmaInit(backend='julia'), indVarying='",
      randomIntercepts,
      "' is fitted through CoTiMA's equivalent latent pseudo-process ",
      "random-intercept formulation. This automatic translation applies ",
      "only to ctmaInit()."
    )
  }
  list(
    indVarying = indVarying,
    randomIntercepts = randomIntercepts,
    backend_args = resolved_backend_args,
    translated = TRUE
  )
}

.ctmaInitJuliaResources <- function(
    cores, estimation, backendArgs, resolved_backend_args,
    announce = TRUE, workflow = "ctmaInit") {
  cores <- suppressWarnings(as.integer(cores)[1L])
  if (is.na(cores) || cores < 1L) cores <- 1L
  explicit_controls <- c(
    "julia_threads", "subject_threads", "parallel_starts", "threads_per_start",
    "parallel_chains", "threads_per_chain"
  )
  if (is.character(resolved_backend_args$julia_threads) &&
      identical(tolower(resolved_backend_args$julia_threads), "auto")) {
    pool <- .ctmaResolveJuliaThreads("auto")
    resolved_backend_args$julia_threads <- pool
    resolved_backend_args$subject_threads <- .ctmaResolveSubjectThreads(
      resolved_backend_args$subject_threads, pool
    )
    description <- sprintf(
      "%d Julia pool thread(s) with %d concurrent subject task(s)",
      pool, resolved_backend_args$subject_threads
    )
    if (isTRUE(announce)) message(
      "Note: ", workflow, " automatically selected ", description,
      " and ", resolved_backend_args$blas_threads, " BLAS thread(s)."
    )
    return(list(
      backend_args = resolved_backend_args, automatic = TRUE, cores = pool
    ))
  }
  if (any(explicit_controls %in% names(backendArgs))) {
    return(list(
      backend_args = resolved_backend_args,
      automatic = FALSE,
      cores = cores
    ))
  }

  estimation <- tolower(estimation)
  if (identical(estimation, "bayes")) {
    chains <- max(1L, as.integer(resolved_backend_args$chains)[1L])
    parallel_chains <- min(chains, cores)
    threads_per_chain <- max(1L, cores %/% parallel_chains)
    resolved_backend_args$julia_threads <- 1L
    resolved_backend_args$parallel_chains <- parallel_chains
    resolved_backend_args$threads_per_chain <- threads_per_chain
    description <- sprintf(
      "%d concurrent chain(s) with %d subject thread(s) per chain",
      parallel_chains, threads_per_chain
    )
  } else {
    starts <- max(1L, as.integer(resolved_backend_args$multistart)[1L])
    resolved_backend_args$julia_threads <- cores
    if (starts > 1L) {
      parallel_starts <- min(starts, cores)
      threads_per_start <- max(1L, cores %/% parallel_starts)
      resolved_backend_args$parallel_starts <- parallel_starts
      resolved_backend_args$threads_per_start <- threads_per_start
      description <- sprintf(
        "%d concurrent start(s) with %d subject thread(s) per start; %d main-process threads",
        parallel_starts, threads_per_start, cores
      )
    } else {
      description <- sprintf("%d Julia subject thread(s)", cores)
    }
  }
  if (isTRUE(announce)) {
    message(
      "Note: ", workflow, " automatically translated coresToUse into ",
      description,
      " and ", resolved_backend_args$blas_threads, " BLAS thread(s)",
      ". Supply an explicit Julia thread control in backendArgs to override ",
      "this allocation."
    )
  }
  list(
    backend_args = resolved_backend_args,
    automatic = TRUE,
    cores = cores
  )
}

.ctmaJuliaThreadNote <- function(
    julia_threads, blas_threads, subject_threads = julia_threads) {
  message(
    "Note: Effective Julia resource configuration: ", julia_threads,
    " Julia pool thread(s), ", subject_threads,
    " concurrent subject task(s), and ", blas_threads,
    " BLAS thread(s). Parallel likelihoods and AD gradients are enabled ",
    "independently only after exact serial-consistency checks; Hessians use ",
    "the serial-safe path."
  )
}

.ctmaInitJuliaFinishsamples <- function(
    finishsamples, backendArgs, resolved_backend_args) {
  if (is.null(finishsamples) || "transformed_draws" %in% names(backendArgs)) {
    return(resolved_backend_args)
  }
  finishsamples <- suppressWarnings(as.integer(finishsamples)[1L])
  if (is.na(finishsamples) || finishsamples < 2L) {
    stop("'finishsamples' must be an integer of at least 2.", call. = FALSE)
  }
  resolved_backend_args$transformed_draws <- finishsamples
  resolved_backend_args
}

ctmaMergeStanctArgs <- function(defaults, overrides = NULL) {
  if (is.null(overrides)) {
    return(defaults)
  }

  if (!is.list(overrides) || is.null(names(overrides)) || any(names(overrides) == "")) {
    stop("CoTiMAStanctArgs must be a fully named list.", call. = FALSE)
  }

  unknown <- setdiff(names(overrides), names(defaults))
  if (length(unknown) > 0L) {
    stop(
      "Unknown CoTiMAStanctArgs: ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  utils::modifyList(defaults, overrides, keep.null = TRUE)
}

# Resolve ctmaInit()'s legacy Stan list through the backend infrastructure.
# Explicit top-level ctmaInit() arguments retain precedence over values in the
# legacy list, matching the historical interface.
.ctmaResolveInitStanArgs <- function(overrides = NULL, scaleTI = NULL,
    scaleTime = NULL, optimize = NULL, priors = NULL, finishsamples = NULL,
    chains = NULL, iter = NULL, verbose = NULL) {
  resolved <- .ctmaMergeBackendArgs(
    ctmaBackendDefaults(backend = "stan"),
    overrides,
    path = "CoTiMAStanctArgs"
  )

  direct <- list(
    scaleTI = scaleTI,
    scaleTime = scaleTime,
    optimize = optimize,
    priors = priors,
    chains = chains,
    iter = iter,
    verbose = verbose
  )
  supplied <- !vapply(direct, is.null, logical(1L))
  resolved[names(direct)[supplied]] <- direct[supplied]
  if (!is.null(finishsamples)) {
    resolved$optimcontrol$finishsamples <- finishsamples
  }
  resolved
}

# Resolve ctmaFit() controls while retaining its historical defaults. Unlike
# optimize and priors, verbose can be inherited from CoTiMAStanctArgs when the
# direct ctmaFit() argument was omitted; `verbose_was_supplied` preserves that
# distinction without changing the public function signature.
.ctmaResolveFitStanArgs <- function(overrides = NULL, optimize = TRUE,
    priors = FALSE, finishsamples = NULL, chains = NULL, iter = NULL,
    verbose = 0, verbose_was_supplied = FALSE) {
  resolved <- .ctmaMergeBackendArgs(
    ctmaBackendDefaults(backend = "stan"),
    overrides,
    path = "CoTiMAStanctArgs"
  )

  direct <- list(
    optimize = optimize,
    priors = priors,
    chains = chains,
    iter = iter
  )
  supplied <- !vapply(direct, is.null, logical(1L))
  resolved[names(direct)[supplied]] <- direct[supplied]
  if (isTRUE(verbose_was_supplied)) resolved$verbose <- verbose
  if (!is.null(finishsamples)) {
    resolved$optimcontrol$finishsamples <- finishsamples
  }
  resolved
}

.ctmaResolveAllInvStanArgs <- function(overrides = NULL, scaleTime = NULL,
    optimize = TRUE, priors = FALSE, finishsamples = NULL, chains = NULL,
    iter = NULL, verbose = NULL, verbose_was_supplied = FALSE) {
  resolved <- .ctmaMergeBackendArgs(
    ctmaBackendDefaults(backend = "stan"),
    overrides,
    path = "CoTiMAStanctArgs"
  )

  direct <- list(
    scaleTime = scaleTime,
    optimize = optimize,
    priors = priors,
    chains = chains,
    iter = iter
  )
  supplied <- !vapply(direct, is.null, logical(1L))
  resolved[names(direct)[supplied]] <- direct[supplied]
  if (isTRUE(verbose_was_supplied) && !is.null(verbose)) {
    resolved$verbose <- verbose
  }
  if (!is.null(finishsamples)) {
    resolved$optimcontrol$finishsamples <- finishsamples
  }
  resolved
}

.ctmaResolvePowerStanArgs <- function(overrides = NULL, scaleTime = NULL,
    optimize = TRUE, priors = FALSE, finishsamples = NULL, chains = NULL,
    iter = NULL, verbose = NULL, verbose_was_supplied = FALSE) {
  resolved <- .ctmaMergeBackendArgs(
    ctmaBackendDefaults(backend = "stan"),
    overrides,
    path = "CoTiMAStanctArgs"
  )

  direct <- list(
    scaleTime = scaleTime,
    optimize = optimize,
    priors = priors,
    chains = chains,
    iter = iter
  )
  supplied <- !vapply(direct, is.null, logical(1L))
  resolved[names(direct)[supplied]] <- direct[supplied]
  if (isTRUE(verbose_was_supplied) && !is.null(verbose)) {
    resolved$verbose <- verbose
  }
  if (!is.null(finishsamples)) {
    resolved$optimcontrol$finishsamples <- finishsamples
  }
  resolved
}

.ctmaResolveOptimizeStanArgs <- function(overrides = NULL,
    inherited = NULL, finishsamples = NULL) {
  defaults <- ctmaBackendDefaults(backend = "stan")
  base <- .ctmaMergeBackendArgs(
    defaults,
    inherited,
    path = "stored CoTiMAStanctArgs"
  )
  resolved <- .ctmaMergeBackendArgs(
    base,
    overrides,
    path = "CoTiMAStanctArgs"
  )
  if (!is.null(finishsamples)) {
    resolved$optimcontrol$finishsamples <- finishsamples
  }
  resolved
}

.ctmaOptimizeFitValue <- function(fit) {
  if (!is.list(fit) || inherits(fit, "condition")) return(NA_real_)
  value <- fit$summary$minus2ll
  if (is.null(value) && inherits(fit, "ctmaJuliaFit")) {
    loglik <- if (!is.null(fit$optimization)) {
      fit$optimization$loglik
    } else {
      fit$likelihood$julia
    }
    value <- -2 * loglik
  }
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    return(NA_real_)
  }
  as.numeric(value)
}

.ctmaFinalizeJuliaCoTiMAFit <- function(fit, setup, cotima_call,
    backend_args) {
  loglik <- if (!is.null(fit$optimization)) {
    fit$optimization$loglik
  } else {
    fit$likelihood$julia
  }
  n_parameters <- length(fit$coefficients$raw)
  fit$summary <- list(
    minus2ll = -2 * as.numeric(loglik),
    loglik = as.numeric(loglik),
    n.parameters = n_parameters
  )
  fit$argumentList <- setup$argumentList
  fit$argumentList$backend <- "julia"
  fit$argumentList$backendArgs <- backend_args
  fit$argumentList$fit <- TRUE
  fit$cotima_call <- cotima_call
  fit$backend <- "julia"
  fit
}

.ctmaJuliaInstallationGuideNote <- function() {
  message(
    "Note: If Julia is not installed or you encounter problems setting up or ",
    "using Julia, open the step-by-step installation guide by entering this ",
    "command in the R console:\n  ",
    "vignette(\"Julia_Install\", package = \"CoTiMA2\")"
  )
  invisible(NULL)
}

.ctmaJuliaCheckProjectDependencies <- function() {
  required <- c("AdvancedHMC", "ForwardDiff", "Optim")
  missing <- as.character(JuliaCall::julia_eval(paste0(
    "filter(name -> Base.find_package(name) === nothing, ",
    "[\"AdvancedHMC\", \"ForwardDiff\", \"Optim\"])"
  )))
  if (length(missing)) {
    stop(
      "The bundled CoTiMA2 Julia project is missing required package(s): ",
      paste(missing, collapse = ", "), ". Restart R, run ",
      "CoTiMA2::ctmaJuliaSetup(), and then retry the fit. Do not install ",
      "these packages in Julia's global environment.",
      call. = FALSE
    )
  }
  invisible(required)
}

.ctmaResolveFitBackend <- function(
    backend = c("stan", "julia", "ctsemjulia"),
    backendArgs = list(), legacy_stan_args = NULL,
    estimation = c("ML", "Bayes"), workflow = "ctmaFit") {
  backend <- match.arg(
    tolower(backend)[1L], c("stan", "julia", "ctsemjulia")
  )
  estimation <- match.arg(tolower(estimation), c("ml", "bayes"))
  if (!is.null(legacy_stan_args)) {
    .Deprecated("backendArgs", package = "CoTiMA2",
      old = "CoTiMAStanctArgs")
  }
  if (is.null(backendArgs)) backendArgs <- list()
  if (!is.list(backendArgs) ||
      (length(backendArgs) && (is.null(names(backendArgs)) ||
        anyNA(names(backendArgs)) || any(names(backendArgs) == "")))) {
    stop("backendArgs must be a fully named list.", call. = FALSE)
  }

  if (!is.null(legacy_stan_args) && length(backendArgs)) {
    stop(
      "Supply either 'CoTiMAStanctArgs' or 'backendArgs', not both. ",
      "Use backendArgs for new code.",
      call. = FALSE
    )
  }

  if (identical(backend, "julia")) {
    .ctmaJuliaInstallationGuideNote()
    if (!is.null(legacy_stan_args)) {
      stop(
        "'CoTiMAStanctArgs' cannot be used with backend = 'julia'.",
        call. = FALSE
      )
    }
    julia_args <- .ctmaMergeBackendArgs(
      ctmaBackendDefaults(
        backend = "julia",
        estimation = estimation,
        workflow = workflow
      ),
      backendArgs,
      path = "backendArgs"
    )
    return(list(
      backend = backend,
      estimation = estimation,
      stan_args = NULL,
      backend_args = julia_args
    ))
  }

  if (identical(backend, "ctsemjulia")) {
    .ctmaCheckCtsemJuliaRequest(
      backend, workflow = workflow,
      supported = identical(workflow, "ctmaFit")
    )
    .ctmaJuliaInstallationGuideNote()
    if (!is.null(legacy_stan_args)) stop(
      "'CoTiMAStanctArgs' cannot be used with backend = 'ctsemjulia'.",
      call. = FALSE
    )
    if (!identical(estimation, "ml")) stop(
      "The experimental ctsemjulia backend currently supports ML/Laplace ",
      "estimation only (use optimize = TRUE).",
      call. = FALSE
    )
    native_args <- .ctmaMergeBackendArgs(
      ctmaBackendDefaults(
        backend = "ctsemjulia", estimation = estimation,
        workflow = workflow
      ),
      backendArgs,
      path = "backendArgs"
    )
    return(list(
      backend = backend,
      estimation = estimation,
      stan_args = NULL,
      backend_args = native_args
    ))
  }

  selected <- if (length(backendArgs)) {
    .ctmaMergeBackendArgs(
      ctmaBackendDefaults(backend = "stan"),
      backendArgs,
      path = "backendArgs"
    )
  } else legacy_stan_args
  list(
    backend = backend,
    estimation = estimation,
    stan_args = selected,
    backend_args = selected
  )
}

.ctmaJuliaFitFromSetup <- function(setup, estimation, backend_args) {
  do.call(
    ctmaJuliaFit,
    c(list(object = setup, estimation = estimation), backend_args)
  )
}

.ctmaInitJuliaStudyFit <- function(
    datalong, model, same_initial_times, stan_args, backend_args,
    estimation, start_values = NULL, digits = 4L) {
  wrapper <- structure(
    list(
      data = datalong,
      ctModel = model,
      argumentList = list(
        fit = FALSE,
        sameInitialTimes = same_initial_times,
        CoTiMAStanctArgs = stan_args,
        backend = "julia",
        backendArgs = backend_args
      )
    ),
    class = "CoTiMAFit"
  )

  args <- backend_args
  if (!is.null(start_values) && length(start_values) &&
      all(is.finite(start_values)) && is.null(args$inits)) {
    args$inits <- as.numeric(start_values)
  }
  # ctmaInit's ordinary result requires transformed uncertainty summaries.
  # ML therefore always computes the Hessian; Bayesian fits use retained NUTS
  # draws. Transformation is performed once below through ctAddJuliaSamples().
  args$uncertainty <- identical(tolower(estimation), "ml")
  args$transformed <- FALSE

  julia_fit <- .ctmaJuliaFitFromSetup(
    setup = wrapper,
    estimation = estimation,
    backend_args = args
  )
  nsamples <- as.integer(args$transformed_draws)[1L]
  if (is.na(nsamples) || nsamples < 2L) nsamples <- 1000L
  transformed_cores <- as.integer(args$transformed_cores)[1L]
  if (is.na(transformed_cores) || transformed_cores < 1L) {
    transformed_cores <- 1L
  }
  julia_fit <- ctAddJuliaSamples(
    fit = julia_fit,
    setup = wrapper,
    nsamples = nsamples,
    cores = transformed_cores,
    seed = args$seed,
    digits = digits
  )
  setup <- .ctmaJuliaGenerateSetup(wrapper)
  julia_fit$standata <- setup$standata
  julia_fit$ctstanmodelbase <- model
  julia_fit$ctsem_setup <- setup
  julia_fit$resultsSummary <- julia_fit$ctsem_summary
  julia_fit$backend <- "julia"
  julia_fit
}

ctmaSelectBestParallelFit <- function(attempts) {
  successful <- vapply(attempts, function(x) isTRUE(x$ok), logical(1))

  if (!any(successful)) {
    errors <- vapply(
      attempts,
      function(x) {
        if (is.null(x$error) || is.na(x$error)) "Unknown fitting error" else x$error
      },
      character(1)
    )
    stop(
      "All parallel fitting attempts failed:\n",
      paste(errors, collapse = "\n"),
      call. = FALSE
    )
  }

  valid <- attempts[successful]
  objectives <- vapply(
    valid,
    function(x) x$fit$stanfit$optimfit$value,
    numeric(1)
  )
  finite <- is.finite(objectives)

  if (!any(finite)) {
    stop("All successful fitting attempts returned non-finite log-posterior values.", call. = FALSE)
  }

  valid <- valid[finite]
  objectives <- objectives[finite]
  selected <- which.max(objectives)

  list(
    fit = valid[[selected]]$fit,
    diagnostics = valid[[selected]][c("warn_hessinv", "warnings", "error")],
    logposterior = objectives[[selected]]
  )
}
