# Internal helpers used by ctmaInit().

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
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    return(NA_real_)
  }
  as.numeric(value)
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
