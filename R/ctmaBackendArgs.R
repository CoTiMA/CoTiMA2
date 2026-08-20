#' Backend-specific default arguments
#'
#' Return a fresh named list of backend controls. Stan remains CoTiMA's
#' default backend. During the compatibility phase, the Stan result is
#' identical to \code{CoTiMAStanctArgs}; Julia defaults are grouped by the public
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
#'   current Julia workflows are \code{"ctmaJuliaFit"} and \code{"ctmaReFit"}.
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
    workflow = c("ctmaJuliaFit", "ctmaReFit")) {
  backend <- match.arg(tolower(backend), c("stan", "julia"))

  # Preserve the established Stan defaults exactly during migration. Assigning
  # the object to a local variable also makes the intended copy-on-modify
  # behavior explicit to callers.
  if (identical(backend, "stan")) {
    defaults <- CoTiMAStanctArgs
    return(defaults)
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
      julia_file = system.file("julia", "CtsemJulia.jl", package = "CoTiMA"),
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
    julia_threads = 1L,
    blas_threads = 1L
  )
  if (identical(estimation, "none")) return(common)

  if (identical(estimation, "ml")) {
    return(c(common, list(
      uncertainty = FALSE,
      maxiters = 100L,
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
      transformed = FALSE,
      transformed_draws = 1000L,
      transformed_cores = 1L
    )))
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

.ctmaValidateBackendArgs <- function(defaults, overrides = NULL,
    path = "backendArgs") {
  if (is.null(overrides)) return(invisible(TRUE))
  if (!is.list(overrides) || is.null(names(overrides)) ||
      anyNA(names(overrides)) || any(names(overrides) == "")) {
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
