#' Benchmark a CoTiMA fitting configuration
#'
#' `ctmaFitConfig()` accepts the model-building arguments of [ctmaFit()] via
#' `...`, creates the identical unfitted setup, and performs a bounded backend
#' pilot. For Julia ML it benchmarks subject concurrency and ForwardDiff chunk
#' sizes at one common raw starting vector. For Stan and ctsem's native Julia
#' backend (`"ctsemjulia"`) it validates and returns the resolved controls and
#' reports setup cost; neither path exposes CoTiMA2's subject/chunk tuning grid.
#'
#' @param ctmaInitFit Object accepted by `ctmaFit(ctmaInitFit=...)`.
#' @param ... Model-building arguments accepted by [ctmaFit()]. Do not supply
#'   `fit`; configuration always uses `fit=FALSE`.
#' @param backend Estimation backend: `"julia"`, `"stan"`, or experimental
#'   `"ctsemjulia"`. The latter requires Charles Driver's compatible ctsem
#'   `juliaFit` development version and currently supports ML only.
#' @param backendArgs Named backend controls. Explicit values are retained as
#'   constraints and the recommended values are merged into this list. Calling
#'   [ctmaBackendDefaults()] without arguments returns Stan defaults only.
#' @param time_limit Approximate pilot budget in seconds. An individual Julia
#'   likelihood or gradient evaluation cannot be interrupted midway.
#' @param candidate_subject_threads Optional positive integer candidates.
#' @param candidate_chunksizes Optional positive ForwardDiff chunk candidates.
#' @param candidate_carefulFit Logical values selecting joint (`FALSE`) and
#'   staged (`TRUE`) short optimizer pilots.
#' @param pilot_iterations Maximum optimizer iterations per pilot stage.
#' @param benchmark_repetitions Number of repeated likelihood and gradient
#'   evaluations per candidate. Medians reduce selection by timing noise.
#' @param warmup Whether to run one untimed Julia gradient to exclude most
#'   compilation cost from the comparison.
#' @param verbose Print progress and the recommendation.
#'
#' @return An object of class `ctmaFitConfig` containing the matching setup,
#'   benchmark table, resolved controls, and recommended `backendArgs`.
#' @export
ctmaFitConfig <- function(
    ctmaInitFit = NULL,
    ...,
    backend = c("julia", "stan", "ctsemjulia"),
    backendArgs = list(),
    time_limit = 600,
    candidate_subject_threads = NULL,
    candidate_chunksizes = NULL,
    candidate_carefulFit = c(FALSE, TRUE),
    pilot_iterations = 5L,
    benchmark_repetitions = 2L,
    warmup = TRUE,
    verbose = TRUE) {
  .ctmaVersionNote()
  .ctma_timer <- .ctmaRegisterExecutionTimer("ctmaFitConfig")
  on.exit(.ctma_timer(), add = TRUE)

  .ctmaCheckCtsemJuliaRequest(
    backend, "ctmaFitConfig", supported = TRUE
  )

  backend <- match.arg(
    tolower(backend), c("julia", "stan", "ctsemjulia")
  )
  time_limit <- as.numeric(time_limit)[1L]
  if (!is.finite(time_limit) || time_limit <= 0) {
    stop("'time_limit' must be a positive number of seconds.", call. = FALSE)
  }
  dots <- list(...)
  forbidden <- intersect(names(dots), c("fit", "backend", "backendArgs"))
  if (length(forbidden)) stop(
    "Supply ", paste(forbidden, collapse = ", "),
    " through ctmaFitConfig's explicit arguments; 'fit' is always FALSE.",
    call. = FALSE
  )
  unknown <- setdiff(names(dots), names(formals(ctmaFit)))
  if (length(unknown)) stop(
    "Unknown ctmaFit argument(s): ", paste(unknown, collapse = ", "),
    call. = FALSE
  )
  if (!is.logical(warmup) || length(warmup) != 1L || is.na(warmup) ||
      !is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("'warmup' and 'verbose' must be TRUE or FALSE.", call. = FALSE)
  }
  candidate_carefulFit <- unique(as.logical(candidate_carefulFit))
  candidate_carefulFit <- candidate_carefulFit[!is.na(candidate_carefulFit)]
  if (!length(candidate_carefulFit)) stop(
    "'candidate_carefulFit' must contain TRUE and/or FALSE.", call. = FALSE
  )
  pilot_iterations <- suppressWarnings(as.integer(pilot_iterations)[1L])
  if (is.na(pilot_iterations) || pilot_iterations < 1L) stop(
    "'pilot_iterations' must be a positive integer.", call. = FALSE
  )
  benchmark_repetitions <- suppressWarnings(
    as.integer(benchmark_repetitions)[1L]
  )
  if (is.na(benchmark_repetitions) || benchmark_repetitions < 1L) stop(
    "'benchmark_repetitions' must be a positive integer.", call. = FALSE
  )

  started <- proc.time()[["elapsed"]]
  setup_timing <- system.time({
    setup <- do.call(
      ctmaFit,
      c(list(
        ctmaInitFit = ctmaInitFit,
        fit = FALSE,
        backend = backend,
        backendArgs = backendArgs
      ), dots)
    )
  })

  estimation <- if (!is.null(dots$optimize) && !isTRUE(dots$optimize)) {
    "Bayes"
  } else "ML"
  resolved <- .ctmaResolveFitBackend(
    backend = backend, backendArgs = backendArgs,
    estimation = estimation, workflow = "ctmaFit"
  )$backend_args

  if (backend %in% c("stan", "ctsemjulia")) {
    note <- if (identical(backend, "stan")) {
      paste0(
        "The pooled Stan/ctsem backend does not expose Julia-style subject ",
        "concurrency or ForwardDiff chunk controls. The supplied Stan controls ",
        "were validated and retained; setup time is reported separately."
      )
    } else {
      paste0(
        "The native ctsem Julia/Laplace backend does not expose CoTiMA2's ",
        "subject-thread or ForwardDiff-chunk benchmarking interface. Its ",
        "backendArgs were validated and retained; setup time is reported ",
        "separately. Control worker capacity with ctmaFit's top-level ",
        "coresToUse argument."
      )
    }
    out <- list(
      call = match.call(), backend = backend, setup = setup,
      setup_timing = setup_timing, time_limit = time_limit,
      time_used = proc.time()[["elapsed"]] - started,
      benchmarks = data.frame(), resolved_backend_args = resolved,
      optimizer_pilots = data.frame(),
      recommended = resolved,
      note = note
    )
    class(out) <- "ctmaFitConfig"
    if (isTRUE(verbose)) print(out)
    return(out)
  }

  if (!identical(tolower(estimation), "ml")) stop(
    "Julia configuration benchmarking currently supports optimize=TRUE/ML. ",
    "Bayesian chain allocation should be configured with parallel_chains and ",
    "threads_per_chain.", call. = FALSE
  )

  pool <- .ctmaResolveJuliaThreads(resolved$julia_threads)
  initial_subjects <- .ctmaResolveSubjectThreads(
    resolved$subject_threads, pool
  )
  initialize_args <- resolved
  initialize_args$estimation <- "none"
  initialize_args$optimize <- NULL
  initialize_args$uncertainty <- FALSE
  initialize_args$transformed <- FALSE
  initialize_args$validate <- FALSE
  initialize_args$performance_warning <- FALSE
  initialize_args$show_trace <- FALSE
  initialize_args$julia_threads <- pool
  initialize_args$subject_threads <- initial_subjects
  initialize_args$inits <- if (is.null(resolved$inits)) NULL else resolved$inits
  initialized <- do.call(
    ctmaJuliaFit,
    c(list(object = setup), initialize_args)
  )
  dimensions <- initialized$dimensions
  nbase <- as.integer(dimensions$nparams)[1L]
  parallel_backend_check <- NULL
  if (pool > 1L && initial_subjects > 1L) {
    JuliaCall::julia_eval(sprintf(
      "CtsemJulia.set_subject_threads!(%d)", initial_subjects
    ))
    parallel_backend_check <- JuliaCall::julia_eval(sprintf(
      paste0(
        "CtsemJulia.validate_parallel_backend!(ctma_rawest,ctma_sd;",
        "chunksize=%d,likelihood_tolerance=1e-12,",
        "gradient_tolerance=1e-10)"
      ),
      as.integer(resolved$chunksize)
    ))
    if (!isTRUE(parallel_backend_check$fully_parallel) && isTRUE(verbose)) {
      message(
        "One or more concurrent backend checks failed; unsafe components use their serial-safe paths."
      )
    }
  }
  parallel_likelihood_safe <- !is.null(parallel_backend_check) &&
    isTRUE(parallel_backend_check$likelihood$safe)
  parallel_gradient_safe <- !is.null(parallel_backend_check) &&
    isTRUE(parallel_backend_check$gradient$safe)

  if (is.null(candidate_subject_threads)) {
    candidate_subject_threads <- unique(c(
      1L, 2L, 4L, 8L, 16L, 32L, initial_subjects, pool
    ))
  }
  candidate_subject_threads <- sort(unique(as.integer(
    candidate_subject_threads
  )))
  candidate_subject_threads <- candidate_subject_threads[
    !is.na(candidate_subject_threads) & candidate_subject_threads >= 1L &
      candidate_subject_threads <= pool
  ]
  if (!parallel_likelihood_safe && !parallel_gradient_safe) {
    candidate_subject_threads <- 1L
    if (isTRUE(verbose) && pool > 1L) message(
      "Concurrent likelihood and gradient validation failed; benchmarking ",
      "the deterministic serial path only."
    )
  }
  if (!length(candidate_subject_threads)) stop(
    "No candidate_subject_threads value lies between 1 and julia_threads.",
    call. = FALSE
  )

  if (is.null(candidate_chunksizes)) {
    divisors <- which(nbase %% seq_len(nbase) == 0L)
    candidate_chunksizes <- unique(c(
      divisors, resolved$chunksize, min(7L, nbase), nbase
    ))
  }
  candidate_chunksizes <- sort(unique(as.integer(candidate_chunksizes)))
  candidate_chunksizes <- candidate_chunksizes[
    !is.na(candidate_chunksizes) & candidate_chunksizes >= 1L &
      candidate_chunksizes <= nbase
  ]
  if (!length(candidate_chunksizes)) stop(
    "No candidate_chunksizes value lies between 1 and the number of base parameters.",
    call. = FALSE
  )

  if (isTRUE(warmup)) {
    if (isTRUE(verbose)) message("Warming the Julia subject gradient...")
    JuliaCall::julia_eval(sprintf(
      "CtsemJulia.set_subject_threads!(%d)", initial_subjects
    ))
    invisible(JuliaCall::julia_eval(sprintf(
      paste0("CtsemJulia.ctsem_gradient_by_subject(ctma_rawest,ctma_sd;",
             "chunksize=%d,progress_every=0,threaded=%s)"),
      as.integer(resolved$chunksize),
      if (initial_subjects > 1L) "true" else "false"
    )))
  }

  grid <- expand.grid(
    subject_threads = candidate_subject_threads,
    chunksize = candidate_chunksizes,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  # Start near the user's/default configuration, then broaden the search.
  grid$priority <- abs(log2(grid$subject_threads / initial_subjects)) +
    abs(grid$chunksize - as.integer(resolved$chunksize)) / max(1L, nbase)
  grid <- grid[order(grid$priority), , drop = FALSE]
  rows <- list()
  for (index in seq_len(nrow(grid))) {
    elapsed <- proc.time()[["elapsed"]] - started
    if (elapsed >= time_limit && length(rows)) break
    subjects <- grid$subject_threads[index]
    chunk <- grid$chunksize[index]
    if (isTRUE(verbose)) message(
      "Benchmarking subject_threads=", subjects,
      ", chunksize=", chunk, "..."
    )
    JuliaCall::julia_eval(sprintf(
      "CtsemJulia.set_subject_threads!(%d)", subjects
    ))
    likelihood_times <- gradient_times <- numeric(benchmark_repetitions)
    for (repetition in seq_len(benchmark_repetitions)) {
      likelihood_times[repetition] <- system.time(
        JuliaCall::julia_eval(sprintf(
          paste0("CtsemJulia.ctsem_loglik_by_subject(ctma_rawest,ctma_sd;",
                 "threaded=%s)"),
          if (subjects > 1L) "true" else "false"
        ))
      )[["elapsed"]]
      gradient_times[repetition] <- system.time(
        JuliaCall::julia_eval(sprintf(
          paste0("CtsemJulia.ctsem_gradient_by_subject(ctma_rawest,ctma_sd;",
                 "chunksize=%d,progress_every=0,threaded=%s)"),
          chunk, if (subjects > 1L) "true" else "false"
        ))
      )[["elapsed"]]
    }
    likelihood <- stats::median(likelihood_times)
    gradient <- stats::median(gradient_times)
    rows[[length(rows) + 1L]] <- data.frame(
      subject_threads = subjects, chunksize = chunk,
      likelihood_seconds = unname(likelihood),
      gradient_seconds = unname(gradient),
      weighted_seconds = unname(gradient + 2 * likelihood),
      repetitions = benchmark_repetitions,
      likelihood_parallel = subjects > 1L && parallel_likelihood_safe,
      gradient_parallel = subjects > 1L && parallel_gradient_safe
    )
  }
  benchmarks <- do.call(rbind, rows)
  best <- benchmarks[which.min(benchmarks$weighted_seconds), , drop = FALSE]
  best_subjects <- as.integer(best$subject_threads)
  best_chunk <- as.integer(best$chunksize)
  JuliaCall::julia_eval(sprintf(
    "CtsemJulia.set_subject_threads!(%d)", best_subjects
  ))
  initial_objective <- -as.numeric(JuliaCall::julia_eval(sprintf(
    paste0("CtsemJulia.ctsem_loglik_by_subject(ctma_rawest,ctma_sd;",
           "threaded=%s)"),
    if (best_subjects > 1L) "true" else "false"
  )))[1L]
  initial_maximum_absolute_gradient <- as.numeric(JuliaCall::julia_eval(
    sprintf(
      paste0("maximum(abs,CtsemJulia.ctsem_gradient_by_subject(",
             "ctma_rawest,ctma_sd;chunksize=%d,progress_every=0,threaded=%s))"),
      best_chunk, if (best_subjects > 1L) "true" else "false"
    )
  ))[1L]
  nti <- as.integer(dimensions$ntipredeffects)[1L]
  if (nti < 1L) candidate_carefulFit <- FALSE
  pilot_rows <- list()
  for (careful in candidate_carefulFit) {
    if (isTRUE(verbose)) message(
      "Running ", pilot_iterations, "-iteration ",
      if (careful) "staged" else "joint", " optimizer pilot..."
    )
    optimizer <- if (careful) {
      "ctsem_optimize_staged"
    } else "ctsem_optimize_auto"
    pilot <- tryCatch({
      timing <- system.time(JuliaCall::julia_eval(sprintf(
        paste0(
          "ctma_config_pilot=CtsemJulia.%s(ctma_rawest,ctma_sd;",
          "maxiters=%d,chunksize=%d,show_trace=false,threaded=%s,",
          "initial_raw_step=%.17g,gradient_tolerance=%.17g,",
          "function_relative_tolerance=%.17g,",
          "parameter_relative_tolerance=%.17g)"
        ),
        optimizer, pilot_iterations, best_chunk,
        if (best_subjects > 1L) "true" else "false",
        resolved$optimizer_initial_step, resolved$gradient_tolerance,
        resolved$function_relative_tolerance,
        resolved$parameter_relative_tolerance
      )))
      diagnostics <- JuliaCall::julia_eval(paste0(
        "(objective=ctma_config_pilot.best_objective,",
        "iterations=ctma_config_pilot.total_iterations,",
        "maximum_absolute_gradient=last(",
        "ctma_config_pilot.stage_maximum_absolute_gradients))"
      ))
      data.frame(
        carefulFit = careful, success = TRUE,
        elapsed = unname(timing[["elapsed"]]),
        iterations = as.integer(diagnostics$iterations),
        objective = as.numeric(diagnostics$objective),
        objective_improvement = initial_objective -
          as.numeric(diagnostics$objective),
        maximum_absolute_gradient =
          as.numeric(diagnostics$maximum_absolute_gradient),
        error = ""
      )
    }, error = function(error) data.frame(
      carefulFit = careful, success = FALSE, elapsed = NA_real_,
      iterations = NA_integer_, objective = NA_real_,
      objective_improvement = NA_real_,
      maximum_absolute_gradient = NA_real_,
      error = conditionMessage(error)
    ))
    pilot_rows[[length(pilot_rows) + 1L]] <- pilot
  }
  optimizer_pilots <- do.call(rbind, pilot_rows)
  optimizer_pilots$improvement_per_second <- with(
    optimizer_pilots,
    ifelse(success & elapsed > 0, objective_improvement / elapsed, NA_real_)
  )
  optimizer_pilots$gradient_reduction_factor <- with(
    optimizer_pilots,
    ifelse(
      success & is.finite(maximum_absolute_gradient) &
        maximum_absolute_gradient > 0,
      initial_maximum_absolute_gradient / maximum_absolute_gradient,
      NA_real_
    )
  )
  optimizer_pilots$convergence_adjusted_score <- with(
    optimizer_pilots,
    improvement_per_second * sqrt(pmax(gradient_reduction_factor, 0))
  )
  successful_pilots <- which(
    optimizer_pilots$success &
      is.finite(optimizer_pilots$convergence_adjusted_score)
  )
  recommended_careful <- if (length(successful_pilots)) {
    optimizer_pilots$carefulFit[successful_pilots[
      which.max(optimizer_pilots$convergence_adjusted_score[successful_pilots])
    ]]
  } else resolved$carefulFit
  recommended <- resolved
  recommended$julia_threads <- pool
  recommended$subject_threads <- best_subjects
  recommended$chunksize <- best_chunk
  recommended$carefulFit <- isTRUE(recommended_careful)

  out <- list(
    call = match.call(), backend = backend, setup = setup,
    setup_timing = setup_timing, time_limit = time_limit,
    time_used = proc.time()[["elapsed"]] - started,
    dimensions = dimensions, benchmarks = benchmarks,
    parallel_backend_check = parallel_backend_check,
    optimizer_pilots = optimizer_pilots,
    resolved_backend_args = resolved, recommended = recommended,
    ranking_rule = "gradient_seconds + 2 * likelihood_seconds",
    optimizer_ranking_rule = paste0(
      "objective improvement per elapsed second, multiplied by the square ",
      "root of the raw-gradient reduction factor"
    ),
    note = paste0(
      "The recommendation ranks evaluation throughput at one common starting ",
      "vector. Complete-fit runtime also depends on optimizer conditioning and ",
      "the number of iterations."
    )
  )
  class(out) <- "ctmaFitConfig"
  if (isTRUE(verbose)) print(out)
  out
}

#' @export
print.ctmaFitConfig <- function(x, ...) {
  cat("ctmaFitConfig\n")
  cat("  backend:                 ", x$backend, "\n", sep = "")
  cat(sprintf("  setup elapsed:           %.3f seconds\n",
              unname(x$setup_timing[["elapsed"]])))
  cat(sprintf("  pilot elapsed:           %.3f seconds\n", x$time_used))
  if (nrow(x$benchmarks)) {
    cat(
      "  recommended backendArgs = list(\n",
      "    julia_threads = ", x$recommended$julia_threads, "L,\n",
      "    subject_threads = ", x$recommended$subject_threads, "L,\n",
      "    chunksize = ", x$recommended$chunksize, "L,\n",
      "    carefulFit = ", if (isTRUE(x$recommended$carefulFit)) "TRUE" else "FALSE", "\n",
      "  )\n",
      sep = ""
    )
    print(x$benchmarks[order(x$benchmarks$weighted_seconds), ], row.names = FALSE)
    if (nrow(x$optimizer_pilots)) {
      print(x$optimizer_pilots, row.names = FALSE)
    }
  }
  cat("  note: ", x$note, "\n", sep = "")
  invisible(x)
}
