# Keep historical defaults usable when newer ctsem removes inert controls.
# Explicit requests for retired behaviour are rejected rather than ignored.
.ctmaRetiredOptimControls <- function(args, checker = NULL) {
  controls <- args$optimcontrol
  if (!is.list(controls)) return(args)
  retired <- list(parstepsAutoModel = FALSE, groupFreeThreshold = 0.5)
  present <- intersect(names(retired), names(controls))
  if (!length(present)) return(args)
  if (is.null(checker)) {
    ns <- asNamespace("ctsem")
    if (!exists(".ctFitCheckControls", ns, inherits = FALSE)) return(args)
    checker <- get(".ctFitCheckControls", ns, inherits = FALSE)
  }
  backend <- if (is.null(args$backend)) "stan" else args$backend
  for (name in present) {
    probe <- tryCatch(checker(stats::setNames(list(controls[[name]]), name),
      backend = backend), error = function(e) e)
    removed <- inherits(probe, "error") &&
      grepl("Unrecognised optimcontrol name", conditionMessage(probe), fixed = TRUE)
    if (!removed) next
    if (!isTRUE(all.equal(controls[[name]], retired[[name]]))) stop(
      "This ctsem version no longer supports optimcontrol$", name,
      ". Remove the retired control; use supported explicit parsteps if needed. ",
      "A non-default request was not silently discarded.", call. = FALSE)
    controls[[name]] <- NULL
  }
  args$optimcontrol <- controls
  args
}

# Internal helpers used by ctmaInit().

.ctmaInitRecordIntervalProblem <- function(problems, study, family, indices,
    parameter_names) {
  if (!length(indices)) return(problems)
  labels <- vapply(indices, function(index) {
    if (index <= length(parameter_names) &&
        !is.na(parameter_names[index]) && nzchar(parameter_names[index])) {
      parameter_names[index]
    } else paste0(family, "[", index, "]")
  }, character(1))
  problems[[length(problems) + 1L]] <- paste0(
    "Possible problems for Study ", study, " in estimating ", family,
    ": ", paste(unique(labels), collapse = " ")
  )
  problems
}

.ctmaInitStudyWarningLabel <- function(index, total, study = NULL) {
  clean_text <- function(x) {
    if (is.null(x) || length(x) == 0L || all(is.na(x))) return(NULL)
    x <- trimws(paste(x, collapse = ", "))
    if (!nzchar(x) || identical(x, "Reference not provided")) return(NULL)
    x
  }

  details <- character()
  study_number <- clean_text(study$originalStudyNo)
  study_source <- clean_text(study$source)
  if (!is.null(study_number)) details <- c(details, paste0("Study: ", study_number))
  if (!is.null(study_source)) details <- c(details, paste0("Source: ", study_source))

  label <- paste0("ctmaInit primary study ", index, " of ", total)
  if (length(details)) label <- paste(label, paste(details, collapse = "; "), sep = "; ")
  paste0("[", label, "] ")
}

.ctmaWithInitStudyWarnings <- function(expr, index, total, study = NULL) {
  captured <- list()
  value <- withCallingHandlers(
    force(expr),
    warning = function(w) {
      captured[[length(captured) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )

  prefix <- .ctmaInitStudyWarningLabel(index, total, study)
  for (w in captured) {
    warning(paste0(prefix, conditionMessage(w)), call. = FALSE)
  }
  value
}

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
  if (any(explicit_controls %in% names(backendArgs))) {
    a <- resolved_backend_args
    a$julia_threads <- if (is.null(a$julia_threads) ||
      (is.character(a$julia_threads) && length(a$julia_threads) == 1L &&
       identical(tolower(a$julia_threads), "auto"))) cores else
        .ctmaResourceInteger(a$julia_threads, "julia_threads")
    bayes <- identical(tolower(estimation), "bayes")
    workers <- if (bayes) a$parallel_chains else a$parallel_starts
    per <- if (bayes) a$threads_per_chain else a$threads_per_start
    if (is.null(workers)) workers <- 1L
    if (is.null(per)) per <- a$julia_threads
    .ctmaCheckJuliaBudget(a$julia_threads, workers, per,
      if (is.null(a$blas_threads)) 1L else a$blas_threads, cores)
    a$subject_threads <- .ctmaResolveSubjectThreads(a$subject_threads, a$julia_threads)
    if (!is.null(a$transformed_cores) && a$transformed_cores > cores)
      stop("transformed_cores exceeds the resolved CPU budget.", call. = FALSE)
    return(list(backend_args = a, automatic = FALSE, cores = cores))
  }

  blas <- if (is.null(resolved_backend_args$blas_threads)) 1L else
    .ctmaResourceInteger(resolved_backend_args$blas_threads, "blas_threads")
  if (blas > cores) stop("blas_threads exceeds the resolved CPU budget.", call. = FALSE)
  thread_budget <- max(1L, cores %/% blas)
  estimation <- tolower(estimation)
  if (identical(estimation, "bayes")) {
    chains <- max(1L, as.integer(resolved_backend_args$chains)[1L])
    parallel_chains <- min(chains, thread_budget)
    threads_per_chain <- max(1L, thread_budget %/% parallel_chains)
    resolved_backend_args$julia_threads <- if (parallel_chains == 1L) threads_per_chain else 1L
    resolved_backend_args$parallel_chains <- parallel_chains
    resolved_backend_args$threads_per_chain <- threads_per_chain
    description <- sprintf(
      "%d concurrent chain(s) with %d subject thread(s) per chain",
      parallel_chains, threads_per_chain
    )
  } else {
    starts <- max(1L, as.integer(resolved_backend_args$multistart)[1L])
    resolved_backend_args$julia_threads <- thread_budget
    if (starts > 1L) {
      parallel_starts <- min(starts, thread_budget)
      threads_per_start <- max(1L, thread_budget %/% parallel_starts)
      resolved_backend_args$parallel_starts <- parallel_starts
      resolved_backend_args$threads_per_start <- threads_per_start
      description <- sprintf(
        "%d concurrent start(s) with %d subject thread(s) per start; %d main-process threads",
        parallel_starts, threads_per_start, thread_budget
      )
    } else {
      description <- sprintf("%d Julia subject thread(s)", thread_budget)
    }
  }
  resolved_backend_args$subject_threads <- .ctmaResolveSubjectThreads(
    resolved_backend_args$subject_threads, resolved_backend_args$julia_threads)
  if (!is.null(resolved_backend_args$transformed_cores) &&
      resolved_backend_args$transformed_cores > cores)
    stop("transformed_cores exceeds the resolved CPU budget.", call. = FALSE)
  if (isTRUE(announce)) {
    message(
      "Note: ", workflow, " automatically translated coresToUse into ",
      description,
      " and ", resolved_backend_args$blas_threads, " BLAS thread(s)",
      ". Supply an explicit Julia thread control in backendArgs to override ",
      "this allocation within the resolved CPU budget."
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

.ctmaSanitizeCtFitArgs <- function(args) {
  if (!is.list(args)) {
    stop("Internal ctFit arguments must be supplied as a list.", call. = FALSE)
  }
  if ("vb" %in% names(args)) {
    if (isTRUE(args$vb)) {
      stop(
        "Variational Bayes is no longer supported by ctsem. Use ",
        "optimize = TRUE for maximum likelihood/MAP estimation or ",
        "optimize = FALSE for Bayesian sampling.",
        call. = FALSE
      )
    }
    args$vb <- NULL
  }
  .ctmaRetiredOptimControls(args)
}

.ctmaCtFitCompat <- function(...) {
  fit <- ctsem::ctFit
  args <- .ctmaBoundCtFitResources(
    .ctmaSanitizeCtFitArgs(.ctmaStationaryFitArgs(list(...)))
  )
  # Older ctsem exposes only Stan and does not accept a backend argument.
  if ("backend" %in% names(formals(fit))) args$backend <- "stan"
  do.call(fit, args)
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
  .ctmaRetiredOptimControls(resolved)
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
  .ctmaRetiredOptimControls(resolved)
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
  .ctmaWrapPooledJuliaFit(fit, setup)
}

.ctmaJuliaInstallationGuideNote <- function() {
  .ctma_print_message(paste0(
    "For help installing or configuring Julia, open the installation guide in R.\n",
    'vignette("Julia_Install", package = "CoTiMA2")'
  ))
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
      "CoTiMA2::ctmaJuliaSetup(agree = TRUE), and then retry the fit. Do not install ",
      "these packages in Julia's global environment.",
      call. = FALSE
    )
  }
  invisible(required)
}

.ctmaCtsemJuliaSampleControl <- function(defaults, overrides = NULL) {
  if (is.null(overrides)) overrides <- list()
  if (!is.list(overrides) || (length(overrides) &&
      (is.null(names(overrides)) || anyNA(names(overrides)) ||
        any(names(overrides) == "")))) {
    stop("backendArgs$sampleControl must be a fully named list.", call. = FALSE)
  }
  accepted <- c(
    "sampler", "placement", "iter", "chains", "warmup", "draws", "maxDraws",
    "processes", "seed", "saveEffects", "minESS", "meanESS",
    "rhatTarget", "maxdepth", "max_treedepth", "target_accept",
    "adapt_delta", "maxdelta", "init_scale", "adapt_metric",
    "adapt_effects", "settleTol", "stepsize", "callback", "particles", "coreset"
  )
  unknown <- setdiff(names(overrides), accepted)
  if (length(unknown)) stop(
    "Unknown backendArgs$sampleControl: ", paste(unknown, collapse = ", "),
    ". Accepted names: ", paste(accepted, collapse = ", "), ".",
    call. = FALSE
  )
  result <- utils::modifyList(defaults, overrides, keep.null = TRUE)
  if (!is.null(result$sampler)) {
    result$sampler <- match.arg(
      tolower(as.character(result$sampler)[1L]), c("nuts", "saem", "map", "smc")
    )
  }
  if (!is.null(result$placement)) {
    result$placement <- match.arg(
      tolower(as.character(result$placement)[1L]), c("fit", "saem")
    )
  }
  if (identical(result$sampler, "smc")) {
    conflicts <- intersect(c("chains", "warmup", "stepsize"),
      names(result)[!vapply(result, is.null, logical(1))])
    if (length(conflicts)) stop(
      "SMC does not use ", paste(conflicts, collapse = ", "),
      "; remove these controls and use particles instead.", call. = FALSE)
  }
  result[!vapply(result, is.null, logical(1))]
}

# Keep the sampling target independent from the kernel that draws from it.
.ctmaCtsemJuliaRoute <- function(backend_args, estimation, has_random_effects) {
  integration <- backend_args$intoverpop
  if (is.null(integration)) integration <- "auto"
  if (is.logical(integration) && length(integration) == 1L && !is.na(integration)) {
    integration <- if (integration) "laplace" else "none"
  }
  integration <- match.arg(integration, c("auto", "none", "laplace", "augmented", "quadrature"))
  sampling <- identical(tolower(estimation), "bayes")
  if (identical(integration, "quadrature") &&
      (sampling || !has_random_effects)) {
    stop("intoverpop = 'quadrature' refines an optimized Laplace fit with ",
      "random effects. Use optimize = TRUE and select random effects; ",
      "for posterior sampling choose another integration route.", call. = FALSE)
  }
  if (sampling) .ctmaCtsemJuliaSampleControl(list(), backend_args$sampleControl)
  sampler <- backend_args$sampleControl$sampler
  if (is.null(sampler)) sampler <- ""
  placement <- backend_args$sampleControl$placement
  if (sampling && sampler %in% c("saem", "map") && !has_random_effects) {
    stop("sampler = '", sampler, "' requires explicit random effects. ",
      "Use automatic sampling or sampler = 'nuts' for this model.", call. = FALSE)
  }
  if (sampling && sampler %in% c("saem", "map") &&
      integration %in% c("laplace", "augmented")) {
    stop("sampler = '", sampler, "' requires intoverpop = 'none' ",
      "(joint sampling). Choose the target explicitly.", call. = FALSE)
  }
  if (sampling && identical(sampler, "smc") && integration == "none") {
    stop("sampler = 'smc' requires a marginal target: intoverpop = ",
      "'laplace' or 'augmented'.", call. = FALSE)
  }
  if (sampling && identical(placement, "saem") &&
      (!has_random_effects || integration %in% c("laplace", "augmented"))) {
    stop("placement = 'saem' requires joint random-effect sampling; ",
      "use placement = 'fit' or omit placement.", call. = FALSE)
  }
  if (sampling && identical(sampler, "map") && identical(placement, "saem")) {
    stop("sampler = 'map' requires placement = 'fit' or automatic placement.",
      call. = FALSE)
  }
  # Explicit kernels must not turn 'auto' into a different posterior target.
  if (sampling && identical(integration, "auto") && identical(sampler, "smc") &&
      has_random_effects) {
    stop("Set intoverpop = 'laplace' or 'augmented' explicitly for SMC ",
      "with random effects.", call. = FALSE)
  }
  integration
}

.ctmaCtsemJuliaMode <- function(result, integration, sampling, has_random_effects) {
  resolved <- result$args$resolved$intoverpop
  if (!is.null(resolved)) integration <- resolved
  if (!has_random_effects) return("direct")
  if (identical(integration, "none")) return("joint")
  if (identical(integration, "auto")) return(if (sampling) "auto" else "laplace")
  integration
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
      supported = workflow %in% c("ctmaFit", "ctmaInit")
    )
    .ctmaJuliaInstallationGuideNote()
    if (!is.null(legacy_stan_args)) stop(
      "'CoTiMAStanctArgs' cannot be used with backend = 'ctsemjulia'.",
      call. = FALSE
    )
    native_defaults <- ctmaBackendDefaults(
      backend = "ctsemjulia", estimation = estimation,
      workflow = workflow
    )
    # ctsem's optimcontrol is an extensible native interface. Validate CoTiMA2's
    # outer argument names strictly, but do not reject a newly added ctsem
    # optimizer control merely because it is absent from this package's visible
    # defaults. Merge the explicit Hessian default recursively so partial user
    # controls retain it.
    optim_overrides <- backendArgs$optimcontrol
    sample_overrides <- backendArgs$sampleControl
    outer_overrides <- backendArgs
    outer_overrides$optimcontrol <- NULL
    outer_overrides$sampleControl <- NULL
    outer_defaults <- native_defaults
    outer_defaults$optimcontrol <- NULL
    outer_defaults$sampleControl <- NULL
    native_args <- .ctmaMergeBackendArgs(
      outer_defaults, outer_overrides, path = "backendArgs"
    )
    if (!is.null(optim_overrides) && (!is.list(optim_overrides) ||
        (length(optim_overrides) && (is.null(names(optim_overrides)) ||
          anyNA(names(optim_overrides)) || any(names(optim_overrides) == ""))))) {
      stop("backendArgs$optimcontrol must be a fully named list.", call. = FALSE)
    }
    native_args$optimcontrol <- if (is.null(optim_overrides)) {
      native_defaults$optimcontrol
    } else {
      utils::modifyList(
        native_defaults$optimcontrol, optim_overrides, keep.null = TRUE
      )
    }
    native_args$intoverpop <- match.arg(
      native_args$intoverpop, c("auto", "none", "laplace", "augmented", "quadrature")
    )
    if (identical(estimation, "bayes")) {
      native_args$sampleControl <- .ctmaCtsemJuliaSampleControl(
        native_defaults$sampleControl, sample_overrides
      )
    } else if (!is.null(sample_overrides)) {
      stop(
        "backendArgs$sampleControl is available only when optimize = FALSE.",
        call. = FALSE
      )
    }
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

# Frame informational engine notices; leave live progress and warnings alone.
.ctmaWithNativeNotices <- function(expr) {
  withCallingHandlers(expr, message = function(condition) {
    text <- trimws(conditionMessage(condition))
    notices <- c("^Bayesian estimation via ",
      "^Maximum (likelihood|a posteriori) estimation requested",
      "^T0VAR rows/columns for latents with individually varying T0MEANS disabled:",
      "^Sampling: placing the sampler first,", "^Starting Julia [0-9]")
    if (any(vapply(notices, grepl, logical(1), x = text))) {
      .ctma_print_message(text)
      invokeRestart("muffleMessage")
    }
  })
}

.ctmaNativeCtsemJuliaFit <- function(args) {
  .ctmaWithNativeNotices(
    do.call(ctsem::ctFit, .ctmaBoundCtFitResources(args)))
}

.ctmaNativeCtsemJuliaSummary <- function(fit, digits) {
  result <- summary(
    fit, digits = 2 * digits, parmatrices = TRUE,
    residualcov = FALSE
  )
  .ctmaInitMaskNativeSummary(fit, result)
}

.ctmaInitCtsemJuliaStudyFit <- function(
    datalong, model, same_initial_times, backend_args,
    start_values = NULL, priors = FALSE, finishsamples = NULL,
    cores = 1L, verbose = 0, digits = 4L, estimation = "ML") {
  .ctmaRequireCtsemJuliaLaplace()
  cores <- as.integer(cores)[1L]
  if (is.na(cores) || cores < 1L) stop(
    "cores must be one positive integer.", call. = FALSE
  )

  fit_args <- list(
    datalong = datalong,
    model = model,
    fit = TRUE,
    sameInitialTimes = same_initial_times,
    stationary = if (is.null(backend_args$stationary)) FALSE else backend_args$stationary,
    backend = "julia",
    intoverpop = .ctmaCtsemJuliaRoute(
      backend_args, estimation,
      is.data.frame(model$pars) && any(model$pars$indvarying %in% TRUE)
    ),
    optimize = identical(tolower(estimation), "ml"),
    cores = cores,
    priors = if (is.null(backend_args$priors)) priors else {
      backend_args$priors
    },
    verbose = if (is.null(backend_args$verbose)) verbose else {
      backend_args$verbose
    }
  )
  init <- if (is.null(backend_args$inits)) start_values else {
    backend_args$inits
  }
  if (!is.null(init) && length(init) && all(is.finite(init))) {
    fit_args$inits <- as.numeric(init)
  }
  optimcontrol <- backend_args$optimcontrol
  if (is.null(optimcontrol)) optimcontrol <- list()
  if (!is.list(optimcontrol)) stop(
    "backendArgs$optimcontrol must be NULL or a named list.",
    call. = FALSE
  )
  if (identical(tolower(estimation), "ml") && !is.null(finishsamples) &&
      is.null(optimcontrol$finishsamples)) {
    optimcontrol$finishsamples <- as.integer(finishsamples)[1L]
  }
  if (length(optimcontrol)) fit_args$optimcontrol <- optimcontrol
  if (!is.null(backend_args$sampleControl)) {
    fit_args$sampleControl <- backend_args$sampleControl
  }

  result <- .ctmaNativeCtsemJuliaFit(fit_args)
  result$resultsSummary <- .ctmaNativeCtsemJuliaSummary(result, digits)
  # ctmaInit's legacy extraction restores `currentModel` from this field.
  # Native ctsem Julia fits retain the model under `$model`, so expose the
  # established compatibility alias used by Stan and CoTiMA2's Julia adapter.
  result$ctstanmodelbase <- model
  result$backend <- "ctsemjulia"
  result$ctsemjulia_mode <- .ctmaCtsemJuliaMode(
    result, fit_args$intoverpop, !fit_args$optimize,
    is.data.frame(model$pars) && any(model$pars$indvarying %in% TRUE)
  )
  .ctmaAnnotateCtsemJuliaResult(result,
    use_laplace = result$ctsemjulia_mode %in% c("laplace", "quadrature"))
}

.ctmaInitCtsemJuliaRandomEffectSummary <- function(results_summary, digits = 4L) {
  popsd <- results_summary$popsd
  popcor <- results_summary$rawpopcorr
  if (!is.data.frame(popsd) || !"mean" %in% names(popsd) ||
      !is.data.frame(popcor) || !"mean" %in% names(popcor)) {
    stop(
      "The native ctsem Julia summary does not contain the expected popsd and rawpopcorr tables.",
      call. = FALSE
    )
  }

  parameter_names <- rownames(popsd)
  dimension <- length(parameter_names)
  correlation_mean <- diag(1, dimension)
  dimnames(correlation_mean) <- list(parameter_names, parameter_names)
  expected_correlations <- dimension * (dimension - 1L) / 2L
  if (NROW(popcor) != expected_correlations) {
    stop(
      "The native ctsem Julia random-effect correlation table has an unexpected size.",
      call. = FALSE
    )
  }
  correlation_mean[lower.tri(correlation_mean)] <- popcor$mean
  correlation_mean[upper.tri(correlation_mean)] <-
    t(correlation_mean)[upper.tri(correlation_mean)]

  standard_deviations <- popsd$mean
  covariance_mean <-
    diag(standard_deviations, dimension) %*%
    correlation_mean %*%
    diag(standard_deviations, dimension)
  dimnames(covariance_mean) <- list(parameter_names, parameter_names)

  unavailable <- matrix(
    NA_real_, dimension, dimension,
    dimnames = list(parameter_names, parameter_names)
  )
  list(
    popsd = popsd,
    popcov_mean = round(covariance_mean, digits),
    popcov_sd = unavailable,
    popcov_T = unavailable,
    popcov_025 = unavailable,
    popcov_50 = unavailable,
    popcov_975 = unavailable,
    popcor_mean = round(correlation_mean, digits),
    popcor_sd = unavailable,
    popcor_T = unavailable,
    popcor_025 = unavailable,
    popcor_50 = unavailable,
    popcor_975 = unavailable
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

# The bundled optimizer has a likelihood target; its NUTS sampler has priors.
# Never silently present an unpenalized optimization as MAP.
.ctmaCheckBundledPriorRequest <- function(backend, optimize, priors) {
  if (identical(tolower(as.character(backend)[1L]), "julia") &&
      isTRUE(as.logical(optimize)[1L]) && isTRUE(as.logical(priors)[1L])) {
    stop("backend = 'julia' does not support MAP optimization (optimize = TRUE, ",
      "priors = TRUE). Use backend = 'ctsemjulia' or 'stan' for MAP; ",
      "use priors = FALSE only if you intend maximum likelihood. ",
      "Bundled Julia NUTS sampling is a separate supported Bayesian method.",
      call. = FALSE)
  }
  invisible(NULL)
}
