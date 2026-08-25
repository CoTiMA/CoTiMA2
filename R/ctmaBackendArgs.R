#' Backend-specific default arguments
#'
#' Return a fresh named list of backend controls. Stan remains CoTiMA's
#' default backend. During the compatibility phase, the Stan result is
#' identical in value to the deprecated \code{CoTiMAStanctArgs} snapshot;
#' Julia defaults are grouped by the public
#' workflow and estimation method that uses them.
#'
#' Model-construction controls such as time scaling and TI-predictor scaling
#' remain top-level CoTiMA arguments. \code{backendArgs} is intended only for
#' controls of the selected estimation engine.
#'
#' @param backend Estimation backend, either \code{"stan"} (the default) or
#'   \code{"julia"}.
#' @param estimation Estimation method. For
#'   \code{workflow = "ctmaJuliaFit"}, use \code{"ML"} or \code{"Bayes"};
#'   \code{"none"} returns the common Julia evaluation controls.
#'   \code{workflow = "ctmaReFit"} currently supports Julia defaults only for
#'   \code{"Bayes"}. This argument is ignored for the compatibility Stan list.
#' @param workflow Public workflow for which defaults are requested. The
#'   current Julia workflows are \code{"ctmaJuliaFit"}, \code{"ctmaInit"},
#'   \code{"ctmaFit"}, \code{"ctmaAllInvFit"}, \code{"ctmaEqual"}, and
#'   \code{"ctmaReFit"}.
#'
#' @return A fresh named list of backend-specific defaults.
#' @export
#'
#' @examples
#' ctmaBackendDefaults()
#' ctmaBackendDefaults("julia", "ML", "ctmaJuliaFit")
#' ctmaBackendDefaults("julia", "Bayes", "ctmaReFit")
ctmaBackendDefaults <- function(
    backend = c("stan", "julia"),
    estimation = NULL,
    workflow = c("ctmaJuliaFit", "ctmaFit", "ctmaInit", "ctmaAllInvFit",
                 "ctmaEqual", "ctmaReFit")) {
  backend <- match.arg(tolower(backend), c("stan", "julia"))

  if (identical(backend, "stan")) {
    return(.ctmaStanBackendDefaults())
  }

  workflow <- match.arg(workflow)
  if (is.null(estimation)) {
    estimation <- if (identical(workflow, "ctmaReFit")) "bayes" else "ml"
  }
  estimation <- match.arg(tolower(estimation), c("none", "ml", "bayes"))

  if (identical(workflow, "ctmaReFit")) {
    if (!identical(estimation, "bayes")) {
      stop(
        "Julia backend defaults for ctmaReFit currently support only ",
        "estimation = 'Bayes'.",
        call. = FALSE
      )
    }
    return(list(
      julia_file = system.file("julia", "CtsemJulia.jl", package = "CoTiMA2"),
      initialize = TRUE,
      julia_home = NULL,
      julia_threads = 1L,
      blas_threads = 1L,
      warmup = 1000L,
      samples = 1000L,
      target_accept = 0.9,
      max_depth = 8L,
      metric = "dense",
      chunksize = 7L,
      parallel_chains = 1L,
      threads_per_chain = NULL,
      nuts_progress = TRUE,
      progress_every = 25L
    ))
  }

  common <- list(
    julia_file = NULL,
    initialize = TRUE,
    julia_home = NULL,
    validate = TRUE,
    tolerance = 1e-6,
    chunksize = 7L,
    inits = NULL,
    seed = NULL,
    performance_warning = TRUE,
    show_trace = TRUE,
    julia_threads = "auto",
    subject_threads = "auto",
    blas_threads = 1L
  )
  if (identical(estimation, "none")) return(common)

  if (identical(estimation, "ml")) {
    result <- c(common, list(
      uncertainty = TRUE,
      maxiters = 5000L,
      gradient_tolerance = 1e-3,
      function_relative_tolerance = 0,
      parameter_relative_tolerance = 0,
      optimizer_initial_step = 0.01,
      multistart = 1L,
      parallel_starts = 1L,
      threads_per_start = NULL,
      multistart_jitter = 0.1,
      carefulFit = TRUE,
      confidence = 0.95,
      transformed = TRUE,
      transformed_draws = 1000L,
      transformed_cores = 1L
    ))
    if (workflow %in% c("ctmaFit", "ctmaInit", "ctmaAllInvFit", "ctmaEqual")) {
      result$uncertainty <- TRUE
      result$transformed <- TRUE
    }
    if (identical(workflow, "ctmaInit")) {
      result$maxiters <- 5000L
      result$gradient_tolerance <- 0.005
    }
    return(result)
  }

  c(common, list(
    chains = 4L,
    warmup = 1000L,
    samples = 1000L,
    target_accept = 0.8,
    max_depth = 10L,
    metric = "diagonal",
    init_jitter = 0.01,
    nuts_progress = TRUE,
    parallel_chains = 1L,
    threads_per_chain = NULL,
    transformed = TRUE,
    transformed_draws = 1000L,
    transformed_cores = 1L
  ))
}

.ctmaAvailableCores <- function(omit = 0L) {
  detected <- if (requireNamespace("parallelly", quietly = TRUE)) {
    suppressWarnings(as.integer(parallelly::availableCores(omit = omit)[1L]))
  } else {
    suppressWarnings(as.integer(parallel::detectCores(logical = TRUE)[1L])) -
      as.integer(omit)
  }
  if (is.na(detected) || detected < 1L) detected <- 1L
  detected
}

.ctmaResolveJuliaThreads <- function(value, available = .ctmaAvailableCores()) {
  if (is.character(value) && length(value) == 1L &&
      identical(tolower(value), "auto")) {
    return(max(1L, as.integer(available) - 1L))
  }
  resolved <- suppressWarnings(as.integer(value)[1L])
  if (is.na(resolved) || resolved < 1L) {
    stop("'julia_threads' must be 'auto' or a positive integer.", call. = FALSE)
  }
  resolved
}

.ctmaResolveSubjectThreads <- function(value, julia_threads) {
  if (is.null(value) || (is.character(value) && length(value) == 1L &&
      identical(tolower(value), "auto"))) return(as.integer(julia_threads))
  resolved <- suppressWarnings(as.integer(value)[1L])
  if (is.na(resolved) || resolved < 1L) {
    stop("'subject_threads' must be 'auto' or a positive integer.", call. = FALSE)
  }
  min(resolved, as.integer(julia_threads))
}

.ctmaResolveCores <- function(
    coresToUse, backend = c("stan", "julia"), task_count = NULL,
    max_stan_workers = 32L, announce = TRUE, workflow = "CoTiMA") {
  backend <- match.arg(tolower(backend), c("stan", "julia"))
  detected <- if (requireNamespace("parallelly", quietly = TRUE)) {
    suppressWarnings(as.integer(parallelly::availableCores()[1L]))
  } else {
    suppressWarnings(as.integer(parallel::detectCores(logical = TRUE)[1L]))
  }
  if (is.na(detected) || detected < 1L) detected <- 1L
  requested <- suppressWarnings(as.integer(coresToUse)[1L])
  if (is.na(requested)) stop("'coresToUse' must be one integer.", call. = FALSE)
  if (requested < 0L) requested <- detected + requested
  if (requested == 0L) requested <- 1L
  resolved <- min(max(1L, requested), detected)
  reason <- character()
  if (!is.null(task_count)) {
    task_count <- suppressWarnings(as.integer(task_count)[1L])
    if (!is.na(task_count) && task_count > 0L && identical(backend, "stan") &&
        resolved > task_count) {
      resolved <- task_count
      reason <- c(reason, "number of parallel tasks")
    }
  }
  max_stan_workers <- suppressWarnings(as.integer(max_stan_workers)[1L])
  if (identical(backend, "stan") && !is.na(max_stan_workers) &&
      max_stan_workers > 0L && resolved > max_stan_workers) {
    resolved <- max_stan_workers
    reason <- c(reason, "safe R-worker limit")
  }
  if (isTRUE(announce) && (resolved != requested || length(reason))) {
    message(
      "Note: ", workflow, " resolved coresToUse=", coresToUse, " to ",
      resolved, " ", if (identical(backend, "stan")) "Stan/R worker" else "Julia thread",
      "(s)", if (length(reason)) paste0(" (limited by ", paste(reason, collapse = " and "), ")") else "", "."
    )
  }
  list(
    requested = as.integer(coresToUse)[1L], resolved = resolved,
    available = detected, backend = backend, reasons = reason
  )
}

.ctmaValidateBackendArgs <- function(defaults, overrides = NULL,
    path = "backendArgs") {
  if (is.null(overrides)) return(invisible(TRUE))
  if (!is.list(overrides) || (length(overrides) &&
      (is.null(names(overrides)) || anyNA(names(overrides)) ||
        any(names(overrides) == "")))) {
    stop(path, " must be a fully named list.", call. = FALSE)
  }

  unknown <- setdiff(names(overrides), names(defaults))
  if (length(unknown)) {
    stop(
      "Unknown ", path, ": ", paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  nested <- intersect(names(overrides), names(defaults))
  for (name in nested) {
    if (is.list(defaults[[name]]) && !is.null(overrides[[name]])) {
      .ctmaValidateBackendArgs(
        defaults[[name]], overrides[[name]],
        path = paste0(path, "$", name)
      )
    }
  }
  invisible(TRUE)
}

.ctmaMergeBackendArgs <- function(defaults, overrides = NULL,
    path = "backendArgs") {
  if (!is.list(defaults) || is.null(names(defaults)) ||
      anyNA(names(defaults)) || any(names(defaults) == "")) {
    stop("Backend defaults must be a fully named list.", call. = FALSE)
  }
  .ctmaValidateBackendArgs(defaults, overrides, path = path)
  if (is.null(overrides)) return(defaults)
  utils::modifyList(defaults, overrides, keep.null = TRUE)
}
