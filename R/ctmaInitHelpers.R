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

.ctmaResolveFitBackend <- function(backend = c("stan", "julia"),
    backendArgs = list(), legacy_stan_args = NULL,
    estimation = c("ML", "Bayes"), workflow = "ctmaFit") {
  backend <- match.arg(tolower(backend), c("stan", "julia"))
  estimation <- match.arg(tolower(estimation), c("ml", "bayes"))
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
