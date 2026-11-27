#' Validate or fit a ctsem model with the experimental Julia backend
#'
#' `ctmaJuliaFit()` accepts either a ctsem fit, a CoTiMA object containing a
#' ctsem fit in `studyFitList`, an unfitted CoTiMA object returned by
#' `ctmaFit(fit=FALSE)`, or the path to an RDS file containing one of these
#' objects. It transfers the generated ctsem `standata` and raw estimates to
#' Julia, checks whether the model uses implemented branches, optionally
#' compares the Julia and Stan log-likelihoods, and can run L-BFGS using the
#' subject-decomposed gradient.
#'
#' @param object A ctsem fit, a fitted or unfitted CoTiMA object, or a path to
#'   an RDS file containing one of these objects.
#' @param study Which element of `studyFitList` to use when it is a list of
#'   ctsem fits. Ignored when `studyFitList` is itself a ctsem fit.
#' @param julia_file Path to `CtsemJulia.jl`. In an installed CoTiMA2 package,
#'   the default looks for `inst/julia/CtsemJulia.jl` via `system.file()`.
#' @param initialize Whether to initialise JuliaCall when it is not already
#'   active. An embedded Julia runtime cannot be restarted inside the same R
#'   process; if its thread count differs from `julia_threads`, restart R.
#' @param julia_home Optional Julia executable directory passed to
#'   `julia_setup()`. When omitted, CoTiMA uses the same platform-aware Julia
#'   discovery as `ctmaJuliaSetup()`.
#' @param validate Whether to compare the Julia likelihood at the retained raw
#'   estimates with `stanfit$optimfit$value`, when available.
#' @param tolerance Maximum absolute Julia-minus-Stan difference accepted by
#'   the validation check.
#' @param optimize Whether to run fast subject-gradient L-BFGS.
#' @param estimation Estimation method: `"none"`, `"ML"`, or `"Bayes"`.
#'   When omitted, `optimize=TRUE` selects ML and otherwise no estimation.
#' @param uncertainty Whether to calculate the observed Hessian, its inverse
#'   covariance matrix, and a normal-theory raw-parameter summary.
#' @param inits Optional raw-scale starting vector. A fitted object's retained
#'   `rawest` is used by default. For `ctmaFit(fit=FALSE)` input, the default is
#'   a correctly sized zero vector.
#' @param maxiters Maximum L-BFGS iterations.
#' @param chunksize ForwardDiff chunk size for subject gradients.
#' @param gradient_tolerance Maximum absolute raw-scale gradient required for
#'   convergence.
#' @param function_relative_tolerance Relative objective-change tolerance.
#'   Zero (the default) disables this secondary rule so it cannot stop before
#'   the raw-gradient criterion.
#' @param parameter_relative_tolerance Relative scaled-parameter tolerance.
#'   Zero (the default) disables this secondary rule.
#' @param optimizer_initial_step Target maximum raw-coordinate displacement in
#'   the automatically preconditioned first L-BFGS direction.
#' @param multistart Number of independent ML optimization starts. The first
#'   uses `inits` exactly; later starts use reproducible relative jitter.
#' @param parallel_starts Number of ML starts run concurrently in separate
#'   local Julia worker processes.
#' @param threads_per_start Subject-likelihood threads assigned to each ML
#'   worker. By default this equals `julia_threads`.
#' @param multistart_jitter Relative raw-scale jitter for starts after the
#'   first. A value of 0.1 means a 10 percent perturbation, with a small floor
#'   for coordinates initialized at zero.
#' @param warm_start \code{"auto"} (the default), \code{TRUE}, or
#'   \code{FALSE}. Before full-data ML, an automatic warm start may optimize
#'   only population/base parameters on a stratified subset of complete
#'   subjects. Explicit \code{inits} disable this pass.
#' @param warm_start_fraction Initial fraction of subjects requested for the
#'   warm-start pilot. Structural minimum rules may enlarge it.
#' @param warm_start_subjects_per_parameter Minimum informative subjects per
#'   population/base parameter.
#' @param warm_start_minimum_per_study Minimum subjects retained from each
#'   study, or all subjects when a study is smaller.
#' @param warm_start_maximum_subjects Soft maximum pilot size. Identification,
#'   study, and binary-category minima take precedence.
#' @param warm_start_binary_minimum Minimum subjects representing each
#'   available zero/one category of every binary manifest.
#' @param warm_start_maxiters Maximum L-BFGS iterations for the subset pilot.
#' @param carefulFit Whether to use staged optimization for models with TI/WEC
#'   effects: optimize base parameters first, then all effects, then perform a
#'   jointly rescaled refinement if needed. The default is `TRUE`, and
#'   `maxiters` applies to each stage.
#' @param staged_optimization Deprecated backward-compatible alias for
#'   `carefulFit`. Leave as `NULL` and use `carefulFit` in new code.
#' @param confidence Confidence level for the raw-parameter intervals.
#' @param transformed Whether to transform approximate raw-parameter draws
#'   through ctsem and construct its population, matrix, and TI-predictor
#'   summary. By default this follows `uncertainty`.
#' @param transformed_draws Number of multivariate-normal raw draws used for
#'   nonlinear ctsem transformations.
#' @param transformed_cores Number of cores passed to ctsem's transformation.
#' @param seed Optional random seed for reproducible transformed summaries.
#' @param chains Number of independent NUTS chains for Bayesian estimation.
#' @param warmup Number of adaptation iterations per NUTS chain.
#' @param samples Number of retained posterior draws per NUTS chain.
#' @param target_accept Target NUTS acceptance probability.
#' @param max_depth Maximum NUTS tree depth.
#' @param metric NUTS mass-matrix structure: `"diagonal"` (the default) or
#'   experimental `"dense"`. A dense metric adapts parameter correlations and
#'   may shorten trajectories, but requires more warmup draws.
#' @param init_jitter Standard deviation of chain-specific perturbations around
#'   `inits` on the raw parameter scale.
#' @param nuts_progress Whether to print AdvancedHMC sampling progress.
#' @param parallel_chains Number of Bayesian chains to run concurrently in
#'   separate local Julia worker processes. One preserves sequential behavior.
#' @param threads_per_chain Julia subject-likelihood threads assigned to each
#'   concurrent chain worker. By default this equals `julia_threads`.
#' @param performance_warning Whether to issue a heuristic advisory when
#'   ctsem's compiled Stan backend may be faster than the experimental Julia
#'   backend for maximum-likelihood estimation.
#' @param show_trace Whether Julia should print the Optim iteration trace.
#' @param julia_threads Number of Julia threads used to evaluate independent
#'   subject likelihoods and gradients. This must be set when Julia starts;
#'   consequently, restart R when changing an active thread pool.
#'   The requested Julia/worker/BLAS combination must fit the available CPU budget.
#' @param subject_threads Maximum number of concurrent subject tasks used by
#'   this fit. Unlike `julia_threads`, this can be changed between fits in the
#'   same R session. `"auto"` uses the complete Julia thread pool. Numeric
#'   likelihoods use these tasks directly; gradients use them only after an
#'   automatic serial-versus-parallel AD consistency check passes. Hessians use
#'   the serial-safe path.
#' @param compute Optional execution plan from \code{ctmaCompute()}.
#' @param blas_threads Number of BLAS threads. One avoids competing with other
#'   jobs and is generally appropriate for the small ctsem matrices here.
#'
#' @return An object of class `ctmaJuliaFit` containing compatibility results,
#'   likelihood validation, timing, model dimensions, and (if requested) the
#'   optimized raw estimates and optimization diagnostics.
#' @importFrom stats coef vcov
#' @importFrom utils getFromNamespace tail
#' @export
ctmaJuliaFit <- function(
    object,
    study = 1L,
    julia_file = NULL,
    initialize = TRUE,
    julia_home = NULL,
    validate = TRUE,
    tolerance = 1e-6,
    optimize = FALSE,
    estimation = NULL,
    uncertainty = TRUE,
    inits = NULL,
    maxiters = 5000L,
    chunksize = 7L,
    gradient_tolerance = 1e-3,
    function_relative_tolerance = 0,
    parameter_relative_tolerance = 0,
    optimizer_initial_step = 0.01,
    multistart = 1L,
    parallel_starts = 1L,
    threads_per_start = NULL,
    multistart_jitter = 0.1,
    warm_start = "auto",
    warm_start_fraction = 0.10,
    warm_start_subjects_per_parameter = 20L,
    warm_start_minimum_per_study = 25L,
    warm_start_maximum_subjects = 5000L,
    warm_start_binary_minimum = 20L,
    warm_start_maxiters = 250L,
    carefulFit = TRUE,
    staged_optimization = NULL,
    confidence = 0.95,
    transformed = uncertainty,
    transformed_draws = 1000L,
    transformed_cores = 1L,
    seed = NULL,
    chains = 4L,
    warmup = 1000L,
    samples = 1000L,
    target_accept = 0.8,
    max_depth = 10L,
    metric = c("diagonal", "dense"),
    init_jitter = 0.01,
    nuts_progress = TRUE,
    parallel_chains = 1L,
    threads_per_chain = NULL,
    performance_warning = TRUE,
    show_trace = TRUE,
    julia_threads = "auto",
    subject_threads = "auto",
    blas_threads = 1L,
    compute = NULL) {
  if (!is.null(compute)) {
    if (!inherits(compute, "ctmaCompute")) stop("Invalid compute plan")
    julia_threads <- subject_threads <- blas_threads <- 1L
    parallel_starts <- threads_per_start <- 1L
  }
  .ctmaVersionNote()
  .ctma_timer <- .ctmaRegisterExecutionTimer("ctmaJuliaFit")
  on.exit(.ctma_timer(), add = TRUE)

  julia_threads_was_supplied <- !missing(julia_threads)
  uncertainty_was_supplied <- !missing(uncertainty)
  transformed_was_supplied <- !missing(transformed)
  user_supplied_inits <- !is.null(inits)

  if (!is.null(staged_optimization)) {
    warning(
      "'staged_optimization' is deprecated; use 'carefulFit' instead.",
      call. = FALSE
    )
    carefulFit <- staged_optimization
  }
  if (!is.logical(carefulFit) || length(carefulFit) != 1L ||
      is.na(carefulFit)) {
    stop("'carefulFit' must be TRUE or FALSE.", call. = FALSE)
  }
  carefulFit <- isTRUE(carefulFit)
  if (!is.logical(performance_warning) || length(performance_warning) != 1L ||
      is.na(performance_warning)) {
    stop("'performance_warning' must be TRUE or FALSE.", call. = FALSE)
  }

  if (is.null(estimation)) {
    estimation <- if (isTRUE(optimize)) "ML" else "none"
  }
  estimation <- match.arg(tolower(estimation), c("none", "ml", "bayes"))
  if (identical(estimation, "none") && !uncertainty_was_supplied) {
    uncertainty <- FALSE
  }
  if (isTRUE(optimize) && estimation != "ml") {
    stop("'optimize=TRUE' is only compatible with estimation='ML'.", call. = FALSE)
  }
  optimize <- identical(estimation, "ml")
  if (!is.null(compute) && !identical(estimation, "ml")) stop("compute currently supports maximum likelihood only")

  if (identical(estimation, "ml")) {
    message("Maximum likelihood estimation requested.")
  } else if (identical(estimation, "bayes")) {
    message("Bayesian NUTS estimation requested.")
  }

  if (!requireNamespace("JuliaCall", quietly = TRUE)) {
    stop("Package 'JuliaCall' is required.", call. = FALSE)
  }

  source_object <- .ctmaJuliaReadObject(object)
  prepared_input <- .ctmaJuliaPrepareInput(
    source_object,
    study = study,
    inits = inits
  )
  fit <- prepared_input$fit
  .ctmaJuliaCheckCovarianceCode(fit$standata)
  setup_only <- prepared_input$setup_only
  if (setup_only && !transformed_was_supplied) {
    transformed <- FALSE
  } else if (identical(estimation, "bayes") && !transformed_was_supplied) {
    transformed <- TRUE
  }

  if (is.null(julia_file)) {
    julia_file <- system.file("julia", "CtsemJulia.jl", package = "CoTiMA2")
  }
  if (length(julia_file) != 1L || !nzchar(julia_file) ||
      !file.exists(julia_file)) {
    stop(
      "Cannot find CtsemJulia.jl. Supply its path through 'julia_file'.",
      call. = FALSE
    )
  }

  julia_file <- .ctmaJuliaRuntimeFile(julia_file)

  requested_julia_threads <- julia_threads
  julia_threads <- .ctmaResolveJuliaThreads(julia_threads)
  cpu_budget <- .ctmaAvailableCores()
  blas_threads <- .ctmaResourceInteger(blas_threads, "blas_threads")
  concurrent <- if (isTRUE(optimize)) parallel_starts else parallel_chains
  per_worker <- if (isTRUE(optimize)) threads_per_start else threads_per_chain
  if (is.null(per_worker)) per_worker <- julia_threads
  .ctmaCheckJuliaBudget(julia_threads, concurrent, per_worker,
    blas_threads, cpu_budget)
  transformed_cores <- .ctmaResourceInteger(transformed_cores, "transformed_cores")
  if (transformed_cores > cpu_budget) stop(
    "transformed_cores exceeds the available CPU budget.", call. = FALSE)


  # This assessment deliberately precedes Julia initialization and the first
  # likelihood evaluation, so the user can still stop an unattractive run.
  dimensions <- .ctmaJuliaDimensionsFromFit(fit)
  warm_start_plan <- list(
    use = FALSE,
    reason = if (user_supplied_inits) "explicit inits were supplied" else
      "not requested for this estimation method",
    subjects = integer(),
    nsubjects = dimensions$nsubjects
  )
  if (isTRUE(optimize) && !user_supplied_inits) {
    warm_start_plan <- .ctmaWarmStartPlan(
      fit = fit, source_object = source_object, mode = warm_start,
      fraction = warm_start_fraction,
      subjects_per_parameter = warm_start_subjects_per_parameter,
      minimum_per_study = warm_start_minimum_per_study,
      maximum_subjects = warm_start_maximum_subjects,
      binary_minimum = warm_start_binary_minimum,
      seed = seed
    )
  } else {
    .ctmaWarmStartMode(warm_start)
  }
  if (isTRUE(optimize) && user_supplied_inits &&
      .ctmaWarmStartMode(warm_start) != "never") {
    message(
      "Note: The subset warm start is ignored because explicit 'inits' ",
      "were supplied."
    )
  } else if (isTRUE(optimize) && isTRUE(warm_start_plan$use)) {
    message(sprintf(
      paste0(
        "Note: Subset warm start selected %d of %d complete subjects ",
        "(%.1f%%; %d study stratum/strata) after structural checks."
      ),
      warm_start_plan$selected, warm_start_plan$nsubjects,
      100 * warm_start_plan$fraction, warm_start_plan$studies
    ))
  } else if (isTRUE(optimize) &&
      .ctmaWarmStartMode(warm_start) != "never") {
    message(
      "Note: Subset warm start was skipped because ",
      warm_start_plan$reason, "."
    )
  }
  performance_assessment <- .ctmaJuliaPerformanceAssessment(
    dimensions = dimensions,
    julia_threads = julia_threads,
    chunksize = chunksize,
    carefulFit = carefulFit
  )
  if (isTRUE(optimize) && isTRUE(performance_warning) &&
      identical(performance_assessment$classification, "stan_may_be_faster")) {
    warning(
      performance_assessment$message,
      call. = FALSE,
      immediate. = TRUE
    )
  }

  julia_is_active <- .ctmaJuliaSessionInitialized()
  if (isTRUE(initialize) && !julia_is_active) {
    julia_home <- dirname(.ctmaJuliaFindExecutable(julia_home = julia_home))
    old_julia_threads <- Sys.getenv("JULIA_NUM_THREADS", unset = NA_character_)
    Sys.setenv(JULIA_NUM_THREADS = as.character(julia_threads))
    on.exit({
      if (is.na(old_julia_threads)) {
        Sys.unsetenv("JULIA_NUM_THREADS")
      } else {
        Sys.setenv(JULIA_NUM_THREADS = old_julia_threads)
      }
    }, add = TRUE)
    JuliaCall::julia_setup(
      JULIA_HOME = julia_home,
      installJulia = FALSE,
      install = FALSE,
      useRCall = FALSE,
      force = FALSE
    )
    julia_is_active <- TRUE
  } else if (!julia_is_active) {
    stop(
      paste0(
        "JuliaCall is not initialized. Use 'initialize=TRUE' for the first ",
        "ctmaJuliaFit() call in this R session."
      ),
      call. = FALSE
    )
  }

  active_julia_threads <- as.integer(
    JuliaCall::julia_eval("Threads.nthreads()")
  )[1L]
  if (julia_threads_was_supplied && julia_threads != active_julia_threads) {
    stop(
      sprintf(
        paste0(
          "The active Julia session has %d thread(s), but ",
          "'julia_threads=%d' was requested. Julia's thread count cannot be ",
          "changed inside an active R process. Restart R and make this the ",
          "first ctmaJuliaFit() call."
        ),
        active_julia_threads,
        julia_threads
      ),
      call. = FALSE
    )
  }
  .ctmaCheckJuliaBudget(active_julia_threads, concurrent, per_worker,
    blas_threads, cpu_budget)
  subject_threads <- .ctmaResolveSubjectThreads(
    subject_threads, active_julia_threads
  )
  use_subject_threads <- subject_threads > 1L

  JuliaCall::julia_assign("ctma_julia_project", dirname(julia_file))
  JuliaCall::julia_eval(
    "import Pkg; Pkg.activate(ctma_julia_project; io=devnull)"
  )
  .ctmaJuliaCheckProjectDependencies()
  JuliaCall::julia_source(julia_file)
  JuliaCall::julia_eval(sprintf(
    "CtsemJulia.set_subject_threads!(%d)", subject_threads
  ))
  JuliaCall::julia_assign("ctma_sd", fit$standata)
  JuliaCall::julia_assign("ctma_rawest", as.numeric(fit$stanfit$rawest))
  JuliaCall::julia_eval(sprintf(
    "LinearAlgebra.BLAS.set_num_threads(%d)",
    as.integer(blas_threads)[1L]
  ))
  .ctmaJuliaThreadNote(
    julia_threads = active_julia_threads,
    subject_threads = subject_threads,
    blas_threads = as.integer(blas_threads)[1L]
  )

  unsupported <- JuliaCall::julia_eval(
    "CtsemJulia.check_supported(ctma_sd; throw_error=false)"
  )
  unsupported <- as.character(unsupported)
  if (length(unsupported)) {
    stop(
      paste(
        "The model uses ctsem branches not yet supported by CtsemJulia:",
        paste0("  - ", unsupported, collapse = "\n"),
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  # Retain the actual active thread count in the stored assessment. The early
  # advisory used the requested count because Julia had not yet been started.
  performance_assessment <- .ctmaJuliaPerformanceAssessment(
    dimensions = dimensions,
    julia_threads = active_julia_threads,
    chunksize = chunksize,
    carefulFit = carefulFit
  )

  parallel_backend_check <- NULL
  parallel_backend_check_timing <- NULL
  if (use_subject_threads && estimation %in% c("ml", "bayes")) {
    parallel_backend_check_timing <- system.time({
      parallel_backend_check <- JuliaCall::julia_eval(sprintf(
        paste0(
          "CtsemJulia.validate_parallel_backend!(ctma_rawest,ctma_sd;",
          "chunksize=%d,likelihood_tolerance=1e-12,",
          "gradient_tolerance=1e-10)"
        ),
        as.integer(chunksize)[1L]
      ))
    })
    if (isTRUE(parallel_backend_check$likelihood$safe)) {
      message(sprintf(
        paste0(
          "Note: Concurrent Julia likelihood validation passed ",
          "(maximum relative difference %.3g)."
        ),
        parallel_backend_check$likelihood$maximum_relative_error
      ))
    } else {
      warning(sprintf(
        paste0(
          "Concurrent Julia likelihood validation failed ",
          "(maximum relative difference %.3g); likelihood evaluations will ",
          "use the deterministic serial path."
        ),
        parallel_backend_check$likelihood$maximum_relative_error
      ), call. = FALSE, immediate. = TRUE)
    }
    if (isTRUE(parallel_backend_check$gradient$safe)) {
      message(sprintf(
        paste0(
          "Note: Concurrent Julia AD gradient validation passed ",
          "(relative difference %.3g); subject-gradient parallelism is enabled."
        ),
        parallel_backend_check$gradient$relative_error
      ))
    } else {
      warning(sprintf(
        paste0(
          "Concurrent Julia AD gradients differed from the exact serial path ",
          "(relative difference %.3g). Gradients and Hessians will use the ",
          "serial-safe path."
        ),
        parallel_backend_check$gradient$relative_error
      ), call. = FALSE, immediate. = TRUE)
    }
  }

  likelihood_timing <- system.time({
    julia_loglik <- JuliaCall::julia_eval(
      paste0(
        "CtsemJulia.ctsem_loglik_by_subject(ctma_rawest,ctma_sd;threaded=",
        if (use_subject_threads) "true" else "false", ")"
      )
    )
  })
  julia_loglik <- as.numeric(julia_loglik)[1L]

  stan_loglik <- NA_real_
  if (!is.null(fit$stanfit$optimfit$value)) {
    stan_loglik <- as.numeric(fit$stanfit$optimfit$value)[1L]
  }
  difference <- julia_loglik - stan_loglik
  validation_passed <- is.na(difference) || abs(difference) <= tolerance

  if (isTRUE(validate) && !is.na(difference) && !validation_passed) {
    warning(
      sprintf(
        paste0(
          "Julia and retained Stan log-likelihoods differ by %.8g ",
          "(tolerance %.8g)."
        ),
        difference,
        tolerance
      ),
      call. = FALSE
    )
  }

  optimization <- NULL
  optimization_timing <- NULL
  warm_start_diagnostics <- warm_start_plan
  warm_start_timing <- NULL
  optimizer_controls <- NULL
  bayes <- NULL
  bayes_timing <- NULL

  if (isTRUE(optimize)) {
    optimizer_controls <- c(
      gradient_tolerance = gradient_tolerance,
      function_relative_tolerance = function_relative_tolerance,
      parameter_relative_tolerance = parameter_relative_tolerance,
      optimizer_initial_step = optimizer_initial_step
    )
    if (any(!is.finite(optimizer_controls)) ||
        gradient_tolerance <= 0 || function_relative_tolerance < 0 ||
        parameter_relative_tolerance < 0 || optimizer_initial_step <= 0) {
      stop(
        paste0(
          "Optimizer tolerances must be finite and non-negative; ",
          "'gradient_tolerance' and 'optimizer_initial_step' must be positive."
        ),
        call. = FALSE
      )
    }
    if (is.null(inits)) {
      inits <- fit$stanfit$rawest
    }
    inits <- as.numeric(inits)
    if (length(inits) != length(fit$stanfit$rawest)) {
      stop(
        sprintf(
          "'inits' has length %d; expected %d.",
          length(inits),
          length(fit$stanfit$rawest)
        ),
        call. = FALSE
      )
    }
    if (any(!is.finite(inits))) {
      stop("'inits' contains non-finite values.", call. = FALSE)
    }

    warm_start_maxiters <- suppressWarnings(
      as.integer(warm_start_maxiters)[1L]
    )
    if (is.na(warm_start_maxiters) || warm_start_maxiters < 1L) stop(
      "'warm_start_maxiters' must be a positive integer.", call. = FALSE
    )
    if (isTRUE(warm_start_plan$use)) {
      warm_start_threaded <- use_subject_threads &&
        !is.null(parallel_backend_check) &&
        isTRUE(parallel_backend_check$fully_parallel)
      JuliaCall::julia_assign(
        "ctma_warm_subjects", as.integer(warm_start_plan$subjects)
      )
      JuliaCall::julia_assign("ctma_warm_original", inits)
      warm_error <- NULL
      warm_result <- NULL
      warm_start_timing <- system.time({
        warm_result <- tryCatch(
          JuliaCall::julia_eval(sprintf(
            paste0(
              "CtsemJulia.ctsem_optimize_warm_start(",
              "ctma_warm_original,ctma_sd,ctma_warm_subjects;",
              "maxiters=%d,show_trace=%s,threaded=%s,",
              "initial_raw_step=%.17g,gradient_tolerance=%.17g,",
              "function_relative_tolerance=%.17g,",
              "parameter_relative_tolerance=%.17g)"
            ),
            warm_start_maxiters,
            if (isTRUE(show_trace)) "true" else "false",
            if (warm_start_threaded) "true" else "false",
            optimizer_initial_step, gradient_tolerance,
            function_relative_tolerance, parameter_relative_tolerance
          )),
          error = function(e) {
            warm_error <<- conditionMessage(e)
            NULL
          }
        )
      })
      if (!is.null(warm_result)) {
        JuliaCall::julia_assign(
          "ctma_warm_candidate", as.numeric(warm_result$estimates)
        )
        full_objectives <- JuliaCall::julia_eval(paste0(
          "begin\n",
          "ordinary=-CtsemJulia.ctsem_loglik_by_subject(",
          "ctma_warm_original,ctma_sd;threaded=",
          if (use_subject_threads) "true" else "false", ")\n",
          "candidate=-CtsemJulia.ctsem_loglik_by_subject(",
          "ctma_warm_candidate,ctma_sd;threaded=",
          if (use_subject_threads) "true" else "false", ")\n",
          "(ordinary=ordinary,candidate=candidate)\nend"
        ))
        accepted <- is.finite(full_objectives$candidate) &&
          (!is.finite(full_objectives$ordinary) ||
             full_objectives$candidate <= full_objectives$ordinary)
        if (accepted) {
          inits <- as.numeric(warm_result$estimates)
          message(sprintf(
            paste0(
              "Note: Subset warm start accepted after full-data validation ",
              "(objective %.6g versus %.6g for the ordinary start)."
            ),
            full_objectives$candidate, full_objectives$ordinary
          ))
        } else {
          message(sprintf(
            paste0(
              "Note: Subset warm start rejected after full-data validation ",
              "(objective %.6g versus %.6g for the ordinary start)."
            ),
            full_objectives$candidate, full_objectives$ordinary
          ))
        }
        warm_start_diagnostics <- c(warm_start_plan, list(
          attempted = TRUE, accepted = accepted,
          subset_starting_objective = warm_result$starting_objective,
          subset_final_objective = warm_result$best_objective,
          subset_iterations = warm_result$iterations,
          subset_initial_maximum_absolute_gradient =
            warm_result$initial_maximum_absolute_gradient,
          subset_maximum_absolute_gradient =
            warm_result$maximum_absolute_gradient,
          full_objective_ordinary = full_objectives$ordinary,
          full_objective_candidate = full_objectives$candidate,
          error = NULL
        ))
      } else {
        message(
          "Note: Subset warm start failed and ordinary starting values will ",
          "be used."
        )
        warm_start_diagnostics <- c(warm_start_plan, list(
          attempted = TRUE, accepted = FALSE, error = warm_error
        ))
      }
    } else {
      warm_start_diagnostics$attempted <- FALSE
      warm_start_diagnostics$accepted <- FALSE
    }

    multistart <- as.integer(multistart)[1L]
    parallel_starts <- as.integer(parallel_starts)[1L]
    if (is.null(threads_per_start)) threads_per_start <- julia_threads
    threads_per_start <- as.integer(threads_per_start)[1L]
    multistart_jitter <- as.numeric(multistart_jitter)[1L]
    if (is.na(multistart) || multistart < 1L ||
        is.na(parallel_starts) || parallel_starts < 1L ||
        is.na(threads_per_start) || threads_per_start < 1L ||
        !is.finite(multistart_jitter) || multistart_jitter < 0) {
      stop(
        paste0(
          "'multistart', 'parallel_starts', and 'threads_per_start' must be ",
          "positive; 'multistart_jitter' must be finite and non-negative."
        ),
        call. = FALSE
      )
    }
    if (parallel_starts > multistart) {
      warning(
        "'parallel_starts' exceeds 'multistart'; using ", multistart,
        " concurrent workers.", call. = FALSE
      )
      parallel_starts <- multistart
    }
    requested_cpus <- parallel_starts * threads_per_start
    available_cpus <- .ctmaAvailableCores()
    if (is.finite(available_cpus) && requested_cpus > available_cpus) warning(
      "Julia ML multistart requests ", requested_cpus, " CPU threads (",
      parallel_starts, " parallel starts x ", threads_per_start,
      " threads per start), but only ", available_cpus,
      " logical cores were detected.", call. = FALSE
    )
    if (multistart == 1L && threads_per_start != active_julia_threads) stop(
      "Single-start ML uses the active Julia thread pool (",
      active_julia_threads, " thread(s)); set 'threads_per_start' to that ",
      "value or restart R.", call. = FALSE
    )

    JuliaCall::julia_assign("ctma_start", inits)
    use_staged_optimization <- carefulFit &&
      as.integer(dimensions$ntipredeffects)[1L] > 0L
    optimizer_function <- if (use_staged_optimization) {
      "ctsem_optimize_staged"
    } else {
      "ctsem_optimize_auto"
    }
    if (!is.null(compute)) {
      optimization_timing <- system.time({
        optimization <- .ctmaComputeOptimize(compute, julia_file, inits,
          maxiters, chunksize, gradient_tolerance, optimizer_initial_step,
          multistart, multistart_jitter, seed, function_relative_tolerance,
          parameter_relative_tolerance, show_trace)
      })
    } else if (multistart > 1L) {
      JuliaCall::julia_assign(
        "ctma_ml_julia_file", normalizePath(julia_file)
      )
      multistart_seed <- if (is.null(seed)) 1L else as.integer(seed)[1L]
      if (is.na(multistart_seed)) stop(
        "'seed' must be coercible to an integer for ML multistart.",
        call. = FALSE
      )
      optimization_timing <- system.time(JuliaCall::julia_eval(sprintf(
        paste0(
          "ctma_ms=CtsemJulia.ctsem_optimize_multistart(ctma_start,ctma_sd;",
          "starts=%d,parallel_starts=%d,threads_per_start=%d,",
          "blas_threads=%d,start_jitter=%.17g,seed=%d,staged=%s,",
          "maxiters=%d,chunksize=%d,initial_raw_step=%.17g,",
          "gradient_tolerance=%.17g,function_relative_tolerance=%.17g,",
          "parameter_relative_tolerance=%.17g,progress=%s,",
          "source_file=ctma_ml_julia_file,project_dir=ctma_julia_project)"
        ),
        multistart, parallel_starts, threads_per_start,
        as.integer(blas_threads)[1L], multistart_jitter, multistart_seed,
        if (use_staged_optimization) "true" else "false",
        as.integer(maxiters)[1L], as.integer(chunksize)[1L],
        optimizer_initial_step, gradient_tolerance,
        function_relative_tolerance, parameter_relative_tolerance,
        if (isTRUE(show_trace)) "true" else "false"
      )))
      optimization <- JuliaCall::julia_eval("ctma_ms.best")
      optimization$staged <- use_staged_optimization
      optimization$selected_start <- as.integer(
        JuliaCall::julia_eval("ctma_ms.best_start")
      )[1L]
      optimization$multistart <- JuliaCall::julia_eval("ctma_ms.diagnostics")
      optimization$multistart_controls <- list(
        starts = multistart,
        parallel_starts = parallel_starts,
        threads_per_start = threads_per_start,
        jitter = multistart_jitter,
        seed = multistart_seed
      )
    } else {
      optimization_timing <- system.time({
        JuliaCall::julia_eval(sprintf(
          paste0(
            "begin\n",
            "ctma_opt=CtsemJulia.%s(",
            paste0(
              "ctma_start,ctma_sd;maxiters=%d,chunksize=%d,",
              "show_trace=%s,threaded=%s,initial_raw_step=%.17g,",
              "gradient_tolerance=%.17g,function_relative_tolerance=%.17g,",
              "parameter_relative_tolerance=%.17g)\n"
            ),
            "nothing\nend"
          ),
          optimizer_function,
          as.integer(maxiters)[1L],
          as.integer(chunksize)[1L],
          if (isTRUE(show_trace)) "true" else "false",
          if (use_subject_threads) "true" else "false",
          optimizer_initial_step,
          gradient_tolerance,
          function_relative_tolerance,
          parameter_relative_tolerance
        ))
      })

      optimization <- JuliaCall::julia_eval(
        paste0(
          "begin\n",
          "x=ctma_opt.estimates\n",
          "g=CtsemJulia.ctsem_gradient_by_subject(",
          "x,ctma_sd;chunksize=", as.integer(chunksize)[1L],
          ",threaded=", if (use_subject_threads) "true" else "false", ")\n",
          "maxg=maximum(abs.(g))\n",
          "(estimates=x,loglik=-ctma_opt.best_objective,",
          "iterations=ctma_opt.total_iterations,",
          "stage_iterations=ctma_opt.stage_iterations,",
          "stage_objectives=ctma_opt.stage_objectives,",
          "stage_maximum_absolute_gradients=",
          "ctma_opt.stage_maximum_absolute_gradients,",
          "staged=", if (use_staged_optimization) "true" else "false", ",",
          "optimizer_converged=CtsemJulia.Optim.converged(ctma_opt.result),",
          "termination_code=string(ctma_opt.result.termination_code),",
          "stopped_by=string(ctma_opt.result.stopped_by),",
          "converged=maxg<=", format(gradient_tolerance, scientific = TRUE), ",",
          "initial_maximum_absolute_gradient=",
          "ctma_opt.initial_maximum_absolute_gradient,",
          "maximum_absolute_gradient=maximum(abs.(g)),",
          "gradient_L2_norm=sqrt(sum(g.^2)),",
          "maximum_parameter_change=maximum(abs.(x .- ctma_start)),",
          "parameter_change_L2=sqrt(sum((x .- ctma_start).^2)),",
          "scale_minimum=minimum(ctma_opt.scale),",
          "scale_maximum=maximum(ctma_opt.scale))\n",
          "end"
        )
      )
    }
    if (!isTRUE(optimization$converged)) {
      reason <- if (tail(optimization$stage_iterations, 1L) >=
                    as.integer(maxiters)[1L]) {
        sprintf(
          "the final stage reached its maximum of %d iterations",
          maxiters
        )
      } else {
        sprintf(
          "Optim stopped with termination code '%s'",
          optimization$termination_code
        )
      }
      warning(
        sprintf(
          paste0(
            "Julia optimization has not met the raw-gradient tolerance: %s; ",
            "maximum absolute gradient is %.6g (target %.6g)."
          ),
          reason,
          optimization$maximum_absolute_gradient,
          gradient_tolerance
        ),
        call. = FALSE
      )
    }
  }

  if (identical(estimation, "bayes")) {
    if (!requireNamespace("posterior", quietly = TRUE)) {
      stop("Package 'posterior' is required for Bayesian diagnostics.", call. = FALSE)
    }
    metric <- match.arg(tolower(metric), c("diagonal", "dense"))
    controls <- c(
      chains = chains, warmup = warmup, samples = samples,
      target_accept = target_accept, max_depth = max_depth,
      init_jitter = init_jitter
    )
    if (any(!is.finite(controls)) || chains < 1 || warmup < 1 || samples < 1 ||
        target_accept <= 0 || target_accept >= 1 || max_depth < 1 ||
        init_jitter < 0) {
      stop("Invalid Bayesian sampling controls.", call. = FALSE)
    }
    if (warmup < 150L) {
      warning(
        paste0(
          "Fewer than 150 warmup iterations do not permit the full Stan-style ",
          "mass-matrix adaptation schedule. This is suitable only for an ",
          "execution/timing test."
        ),
        call. = FALSE
      )
    }
    chains <- as.integer(chains)[1L]
    warmup <- as.integer(warmup)[1L]
    samples <- as.integer(samples)[1L]
    max_depth <- as.integer(max_depth)[1L]
    parallel_chains <- as.integer(parallel_chains)[1L]
    if (is.null(threads_per_chain)) threads_per_chain <- julia_threads
    threads_per_chain <- as.integer(threads_per_chain)[1L]
    if (is.na(parallel_chains) || parallel_chains < 1L ||
        is.na(threads_per_chain) || threads_per_chain < 1L) {
      stop("'parallel_chains' and 'threads_per_chain' must be positive integers.",
           call. = FALSE)
    }
    if (parallel_chains > chains) {
      warning(
        "'parallel_chains' exceeds 'chains'; using ", chains,
        " concurrent chain workers.", call. = FALSE
      )
      parallel_chains <- chains
    }
    requested_cpus <- parallel_chains * threads_per_chain
    available_cpus <- .ctmaAvailableCores()
    if (is.finite(available_cpus) && requested_cpus > available_cpus) warning(
      "Julia chain parallelism requests ", requested_cpus, " CPU threads (",
      parallel_chains, " parallel chains x ", threads_per_chain,
      " threads per chain), but only ", available_cpus,
      " logical cores were detected.", call. = FALSE
    )
    if (parallel_chains == 1L && threads_per_chain != active_julia_threads) {
      stop(
        "Sequential chains use the active Julia thread pool (",
        active_julia_threads, " thread(s)); set 'threads_per_chain' to that ",
        "value or restart R.", call. = FALSE
      )
    }
    bayes_seed <- if (is.null(seed)) 1L else as.integer(seed)[1L]
    if (is.na(bayes_seed)) {
      stop("'seed' must be coercible to an integer.", call. = FALSE)
    }
    bayes_inits <- if (is.null(inits)) fit$stanfit$rawest else inits
    bayes_inits <- as.numeric(bayes_inits)
    if (length(bayes_inits) != length(fit$stanfit$rawest) ||
        any(!is.finite(bayes_inits))) {
      stop("'inits' has the wrong length or contains non-finite values.", call. = FALSE)
    }
    # A dense D-dimensional metric estimates D * (D + 1) / 2 covariance
    # elements. The former 2 * D advisory was much too optimistic (for
    # example, it recommended only 160 warmup iterations for 80 parameters).
    # This remains a pragmatic diagnostic threshold rather than a guarantee of
    # adequate adaptation; inferential adequacy must be judged from multiple
    # chains and sampler diagnostics.
    dense_warmup_advisory <- max(500L, 5L * length(bayes_inits))
    if (identical(metric, "dense") && warmup < dense_warmup_advisory) {
      warning(
        sprintf(
          paste0(
            "Dense-metric adaptation has %d parameters but only %d warmup ",
            "iterations. This is suitable for a timing/geometry test; consider ",
            "at least %d warmup iterations for a subsequent inferential run."
          ),
          length(bayes_inits), warmup, dense_warmup_advisory
        ),
        call. = FALSE
      )
    }
    JuliaCall::julia_assign("ctma_bayes_inits", bayes_inits)
    JuliaCall::julia_assign("ctma_bayes_julia_file", normalizePath(julia_file))
    bayes_timing <- system.time({
      JuliaCall::julia_eval(sprintf(
        paste0(
          "begin\n",
          "ctma_bayes=CtsemJulia.ctsem_nuts_parallel(",
          "ctma_bayes_inits,ctma_sd;",
          "samples=%d,warmup=%d,chains=%d,target_accept=%.17g,",
          "max_depth=%d,chunksize=%d,seed=%d,",
          "init_jitter=%.17g,progress=%s,metric=:%s,",
          "parallel_chains=%d,threads_per_chain=%d,blas_threads=%d,",
          "source_file=ctma_bayes_julia_file,",
          "project_dir=ctma_julia_project)\n",
          "nothing\nend"
        ),
        samples, warmup, chains, target_accept, max_depth,
        as.integer(chunksize)[1L],
        bayes_seed, init_jitter,
        if (isTRUE(nuts_progress)) "true" else "false",
        metric, parallel_chains, threads_per_chain,
        as.integer(blas_threads)[1L]
      ))
    })
    raw_draws_array <- JuliaCall::julia_eval("ctma_bayes.draws")
    dimnames(raw_draws_array) <- list(
      iteration = seq_len(samples),
      parameter = paste0("raw_", seq_len(dim(raw_draws_array)[2L])),
      chain = paste0("chain_", seq_len(chains))
    )
    posterior_array <- aperm(raw_draws_array, c(1L, 3L, 2L))
    posterior_draws <- posterior::as_draws_array(posterior_array)
    diagnostics <- posterior::summarise_draws(posterior_draws)
    sampler_diagnostics <- list(
      acceptance_rate = JuliaCall::julia_eval("ctma_bayes.acceptance_rate"),
      tree_depth = JuliaCall::julia_eval("ctma_bayes.tree_depth"),
      n_steps = JuliaCall::julia_eval("ctma_bayes.n_steps"),
      divergent = JuliaCall::julia_eval("ctma_bayes.divergent"),
      warmup_divergent = JuliaCall::julia_eval(
        "ctma_bayes.warmup_divergent"
      ),
      step_size = JuliaCall::julia_eval("ctma_bayes.step_size")
    )
    bayes <- list(
      draws = raw_draws_array,
      posterior_draws = posterior_draws,
      summary = diagnostics,
      sampler_diagnostics = sampler_diagnostics,
      controls = list(
        chains = chains, warmup = warmup, samples = samples,
        target_accept = target_accept, max_depth = max_depth,
        metric = metric, init_jitter = init_jitter, seed = bayes_seed,
        parallel_chains = parallel_chains,
        threads_per_chain = threads_per_chain
      ),
      divergences = sum(sampler_diagnostics$divergent),
      warmup_divergences = sum(sampler_diagnostics$warmup_divergent),
      maximum_tree_depth_hits = sum(
        sampler_diagnostics$tree_depth >= max_depth
      )
    )
  }

  information <- NULL
  information_timing <- NULL
  raw_summary <- NULL
  transformed_result <- NULL
  cotima_summary <- NULL

  if (isTRUE(uncertainty) && !is.null(optimization) &&
      !isTRUE(optimization$converged)) {
    warning(
      paste0(
        "Raw-gradient convergence was not achieved, so Hessian uncertainty ",
        "was skipped. Increase 'maxiters', improve 'inits', and inspect ",
        "result$optimization before requesting uncertainty again."
      ),
      call. = FALSE
    )
    uncertainty <- FALSE
    transformed <- FALSE
  }

  # Posterior uncertainty comes from retained NUTS draws. A likelihood Hessian
  # at the starting point is neither needed nor a posterior covariance estimate.
  if (isTRUE(uncertainty) && is.null(bayes)) {
    confidence <- as.numeric(confidence)[1L]
    if (!is.finite(confidence) || confidence <= 0 || confidence >= 1) {
      stop("'confidence' must lie strictly between zero and one.", call. = FALSE)
    }
    uncertainty_point <- if (is.null(optimization)) {
      as.numeric(fit$stanfit$rawest)
    } else {
      as.numeric(optimization$estimates)
    }
    JuliaCall::julia_assign("ctma_uncertainty_point", uncertainty_point)
    information_timing <- system.time({
      JuliaCall::julia_eval(sprintf(
        paste0(
          "begin\n",
          "ctma_information=CtsemJulia.ctsem_vcov_by_subject(",
          "ctma_uncertainty_point,ctma_sd;chunksize=%d,",
          "progress_every=1000,threaded=%s)\n",
          "nothing\nend"
        ),
        as.integer(chunksize)[1L],
        if (use_subject_threads) "true" else "false"
      ))
    })
    covariance <- JuliaCall::julia_eval("ctma_information.covariance")
    hessian <- JuliaCall::julia_eval("ctma_information.hessian")
    curvature <- JuliaCall::julia_eval(
      paste0(
        "(minimum_eigenvalue=ctma_information.minimum_eigenvalue,",
        "maximum_eigenvalue=ctma_information.maximum_eigenvalue,",
        "condition_number=ctma_information.condition_number,",
        "number_nonpositive=ctma_information.number_nonpositive,",
        "number_regularized=ctma_information.number_regularized,",
        "eigenvalue_tolerance=ctma_information.eigenvalue_tolerance,",
        "regularized_minimum_eigenvalue=",
        "ctma_information.regularized_minimum_eigenvalue)"
      )
    )
    if (curvature$number_regularized > 0L) {
      warning(
        sprintf(
          paste0(
            "%d near-zero Hessian eigenvalue(s) were regularized at %.6g ",
            "before inversion; the unregularized minimum was %.6g."
          ),
          curvature$number_regularized,
          curvature$eigenvalue_tolerance,
          curvature$minimum_eigenvalue
        ),
        call. = FALSE
      )
    }
    covariance <- as.matrix(covariance)
    hessian <- as.matrix(hessian)
    standard_error <- sqrt(diag(covariance))
    critical <- stats::qnorm(1 - (1 - confidence) / 2)
    raw_summary <- data.frame(
      parameter_index = seq_along(uncertainty_point),
      estimate = uncertainty_point,
      standard_error = standard_error,
      lower = uncertainty_point - critical * standard_error,
      median = uncertainty_point,
      upper = uncertainty_point + critical * standard_error,
      z = uncertainty_point / standard_error,
      row.names = NULL
    )
    information <- list(
      hessian = hessian,
      covariance = covariance,
      curvature = curvature,
      confidence = confidence,
      raw_summary = raw_summary
    )
  }

  if (isTRUE(transformed)) {
    if (setup_only) {
      stop(
        paste0(
          "ctsem transformed summaries require a compiled ctsem fit and are ",
          "not available directly from ctmaFit(fit=FALSE) input. Use ",
          "'transformed=FALSE'; Julia optimization and raw-scale uncertainty ",
          "remain available."
        ),
        call. = FALSE
      )
    }
    if (is.null(information) && is.null(bayes)) {
      stop(
        "Transformed summaries require Hessian uncertainty or Bayesian draws.",
        call. = FALSE
      )
    }
    transformed_draws <- as.integer(transformed_draws)[1L]
    transformed_cores <- as.integer(transformed_cores)[1L]
    if (is.na(transformed_draws) || transformed_draws < 2L) {
      stop("'transformed_draws' must be an integer of at least 2.", call. = FALSE)
    }
    if (is.na(transformed_cores) || transformed_cores < 1L) {
      stop("'transformed_cores' must be a positive integer.", call. = FALSE)
    }
    if (!is.null(bayes)) {
      bayes_raw_matrix <- do.call(
        rbind,
        lapply(seq_len(chains), function(chain) bayes$draws[, , chain])
      )
      transformed_result <- .ctmaJuliaTransformSummary(
        fit = fit,
        estimates = colMeans(bayes_raw_matrix),
        covariance = stats::cov(bayes_raw_matrix),
        loglik = julia_loglik,
        draws = nrow(bayes_raw_matrix),
        cores = transformed_cores,
        seed = seed,
        raw_draws = bayes_raw_matrix
      )
    } else {
      transformed_result <- .ctmaJuliaTransformSummary(
        fit = fit,
        estimates = uncertainty_point,
        covariance = covariance,
        loglik = if (is.null(optimization)) julia_loglik else optimization$loglik,
        draws = transformed_draws,
        cores = transformed_cores,
        seed = seed
      )
    }
    cotima_summary <- .ctmaJuliaCoTiMAWEC(
      object = source_object,
      fit = fit,
      raw_draws = transformed_result$draws
    )
  }

  out <- list(
    call = match.call(),
    # Retain the originating CoTiMA/ctsem object so downstream CoTiMA
    # workflows can recover the model and original data while using Julia's
    # fitted raw parameters as starting values. This is especially important
    # for ctmaRePrep(), because a ctmaJuliaFit object is an estimation wrapper
    # rather than a complete ctsem model setup by itself.
    source_object = source_object,
    estimation = estimation,
    compatible = TRUE,
    setup_only = setup_only,
    starting_values = prepared_input$starting_values,
    unsupported = unsupported,
    dimensions = dimensions,
    performance_assessment = performance_assessment,
    parallel_backend_check = parallel_backend_check,
    parallel_backend_check_timing = parallel_backend_check_timing,
    likelihood = list(
      julia = julia_loglik,
      stan = stan_loglik,
      difference = difference,
      tolerance = tolerance,
      passed = validation_passed
    ),
    likelihood_timing = likelihood_timing,
    optimization = optimization,
    optimization_timing = optimization_timing,
    warm_start = warm_start_diagnostics,
    warm_start_timing = warm_start_timing,
    optimizer_controls = optimizer_controls,
    bayes = bayes,
    bayes_timing = bayes_timing,
    information = information,
    information_timing = information_timing,
    coefficients = list(
      raw = if (is.null(optimization)) {
        if (is.null(bayes)) {
          as.numeric(fit$stanfit$rawest)
        } else {
          apply(bayes$draws, 2L, mean)
        }
      } else {
        as.numeric(optimization$estimates)
      },
      raw_summary = raw_summary,
      transformed = if (is.null(transformed_result)) {
        NULL
      } else {
        transformed_result$parameters
      },
      population = if (is.null(transformed_result)) {
        NULL
      } else {
        transformed_result$summary$popmeans
      },
      tipred = if (is.null(transformed_result)) {
        NULL
      } else {
        transformed_result$summary$tipreds
      },
      matrices = if (is.null(transformed_result)) {
        NULL
      } else {
        transformed_result$summary$parmatrices
      },
      study = if (is.null(cotima_summary)) {
        NULL
      } else {
        cotima_summary$study_specific_parameters
      }
    ),
    transformed_draws = if (is.null(transformed_result)) {
      NULL
    } else {
      transformed_result$draws
    },
    ctsem_summary = if (is.null(transformed_result)) {
      NULL
    } else {
      transformed_result$summary
    },
    cotima_summary = cotima_summary,
    threads = list(
      julia = active_julia_threads,
      subject = subject_threads,
      blas = as.integer(blas_threads)[1L],
      subject_parallel = use_subject_threads
    ),
    julia_file = normalizePath(julia_file, mustWork = TRUE)
  )
  class(out) <- "ctmaJuliaFit"
  out
}

#' @export
print.ctmaJuliaFit <- function(x, ...) {
  cat("ctmaJuliaFit\n")
  cat("  compatible: yes\n")
  cat(sprintf("  estimation:             %s\n", x$estimation))
  if (isTRUE(x$setup_only)) {
    cat(sprintf(
      "  input:                  ctmaFit(fit=FALSE), %s initialization\n",
      x$starting_values
    ))
  }
  cat(sprintf("  Julia log-likelihood: %.10f\n", x$likelihood$julia))
  if (!is.na(x$likelihood$stan)) {
    cat(sprintf("  Stan log-likelihood:  %.10f\n", x$likelihood$stan))
    cat(sprintf("  difference:            %.3g\n", x$likelihood$difference))
    cat(sprintf(
      "  validation:            %s\n",
      if (isTRUE(x$likelihood$passed)) "passed" else "FAILED"
    ))
  }
  cat(sprintf(
    "  likelihood elapsed:    %.3f seconds\n",
    unname(x$likelihood_timing[["elapsed"]])
  ))
  cat(sprintf(
    "  Julia pool threads:     %d\n",
    x$threads$julia
  ))
  cat(sprintf(
    "  Subject tasks:          %d%s\n",
    if (is.null(x$threads$subject)) x$threads$julia else x$threads$subject,
    if (isTRUE(x$threads$subject_parallel)) " (parallel)" else ""
  ))
  if (!is.null(x$optimization)) {
    if (!is.null(x$optimization$multistart_controls)) cat(sprintf(
      "  ML multistart:            %d starts, %d concurrent (%d thread(s) each)\n",
      x$optimization$multistart_controls$starts,
      x$optimization$multistart_controls$parallel_starts,
      x$optimization$multistart_controls$threads_per_start
    ))
    if (!is.null(x$optimization$selected_start)) cat(sprintf(
      "  selected ML start:        %d\n",
      x$optimization$selected_start
    ))
    cat(sprintf(
      "  optimized log-likelihood: %.10f\n",
      x$optimization$loglik
    ))
    cat(sprintf("  iterations:               %d\n", x$optimization$iterations))
    if (isTRUE(x$optimization$staged)) {
      cat(sprintf(
        "  stage iterations:          %s\n",
        paste(x$optimization$stage_iterations, collapse = " + ")
      ))
      cat(sprintf(
        "  stage maximum gradients:   %s\n",
        paste(
          format(x$optimization$stage_maximum_absolute_gradients, digits = 5),
          collapse = ", "
        )
      ))
    }
    cat(sprintf(
      "  converged:                %s\n",
      if (isTRUE(x$optimization$converged)) "yes" else "no"
    ))
    cat(sprintf(
      "  maximum absolute gradient: %.6g\n",
      x$optimization$maximum_absolute_gradient
    ))
    cat(sprintf(
      "  gradient tolerance:         %.6g\n",
      x$optimizer_controls[["gradient_tolerance"]]
    ))
    cat(sprintf(
      "  automatic scale range:     %.3g to %.3g\n",
      x$optimization$scale_minimum,
      x$optimization$scale_maximum
    ))
    cat(sprintf(
      "  optimizer termination:     %s\n",
      x$optimization$termination_code
    ))
  }
  if (!is.null(x$information)) {
    cat(sprintf(
      "  Hessian condition number: %.6g\n",
      x$information$curvature$condition_number
    ))
    cat(sprintf(
      "  uncertainty elapsed:      %.3f seconds\n",
      unname(x$information_timing[["elapsed"]])
    ))
    if (x$information$curvature$number_regularized > 0L) {
      cat(sprintf(
        "  Hessian eigenvalues regularized: %d (floor %.6g)\n",
        x$information$curvature$number_regularized,
        x$information$curvature$eigenvalue_tolerance
      ))
    }
  }
  if (!is.null(x$bayes)) {
    cat(sprintf(
       "  NUTS draws:             %d x %d chains\n",
      x$bayes$controls$samples, x$bayes$controls$chains
    ))
    cat(sprintf("  NUTS warmup:            %d per chain\n", x$bayes$controls$warmup))
    stored_metric <- x$bayes$controls$metric
    if (is.null(stored_metric)) stored_metric <- "diagonal"
    cat(sprintf("  NUTS metric:            %s\n", stored_metric))
    if (!is.null(x$bayes$controls$parallel_chains)) cat(sprintf(
      "  concurrent chain workers: %d (%d thread(s) each)\n",
      x$bayes$controls$parallel_chains,
      x$bayes$controls$threads_per_chain
    ))
    cat(sprintf("  divergences:            %d\n", x$bayes$divergences))
    cat(sprintf(
      "  warmup divergences:     %d\n",
      x$bayes$warmup_divergences
    ))
    cat(sprintf(
      "  maximum tree-depth hits: %d\n",
      x$bayes$maximum_tree_depth_hits
    ))
    if ("rhat" %in% names(x$bayes$summary)) {
      cat(sprintf("  maximum R-hat:          %.4f\n", max(x$bayes$summary$rhat)))
    }
    if ("ess_bulk" %in% names(x$bayes$summary)) {
      cat(sprintf(
        "  minimum bulk ESS:       %.1f\n",
        min(x$bayes$summary$ess_bulk)
      ))
    }
    cat(sprintf(
      "  Bayesian elapsed:       %.3f seconds\n",
      unname(x$bayes_timing[["elapsed"]])
    ))
  }
  if (!is.null(x$ctsem_summary)) {
    cat(sprintf(
      "  transformed draws:       %d\n",
      nrow(x$transformed_draws)
    ))
    cat("  ctsem summary:           available via summary()\n")
  }
  if (!is.null(x$cotima_summary)) {
    cat(sprintf(
      "  study-specific parameter sets: %d (coef(type=\"study\"))\n",
      length(x$cotima_summary$study_specific_parameters)
    ))
  }
  invisible(x)
}

#' Extract coefficients from a Julia CoTiMA fit
#'
#' @param object A fitted object returned by \code{ctmaJuliaFit()}.
#' @param type Coefficient representation to return.
#' @param ... Additional arguments reserved for method compatibility.
#' @return The requested coefficient vector, table, matrices, or study-specific
#'   parameter collection.
#' @method coef ctmaJuliaFit
#' @export
coef.ctmaJuliaFit <- function(
    object,
    type = c(
      "raw", "population", "tipred", "matrices", "study", "transformed"
    ),
    ...) {
  type <- match.arg(type)
  value <- object$coefficients[[type]]
  if (is.null(value) && type != "raw") {
    stop(
      paste0(
        "No transformed coefficients are stored; refit with ",
        "uncertainty=TRUE and transformed=TRUE."
      ),
      call. = FALSE
    )
  }
  value
}

#' Extract the covariance matrix from a Julia CoTiMA fit
#'
#' @param object A fitted object returned by \code{ctmaJuliaFit()}.
#' @param ... Additional arguments reserved for method compatibility.
#' @return A covariance matrix calculated from Bayesian draws or the observed
#'   information matrix.
#' @method vcov ctmaJuliaFit
#' @export
vcov.ctmaJuliaFit <- function(object, ...) {
  if (!is.null(object$bayes)) {
    matrix_draws <- do.call(
      rbind,
      lapply(
        seq_len(object$bayes$controls$chains),
        function(chain) object$bayes$draws[, , chain]
      )
    )
    return(stats::cov(matrix_draws))
  }
  if (is.null(object$information$covariance)) {
    stop(
      "No covariance matrix is stored; refit with uncertainty=TRUE.",
      call. = FALSE
    )
  }
  object$information$covariance
}

#' Add transformed ctsem samples to a Julia fit
#'
#' Mirrors the role of `ctsem::ctAddSamples()`, but starts from a completed
#' `ctmaJuliaFit`. For ML fits it draws from the stored Julia Hessian-normal
#' approximation; for Bayesian fits it uses the actual retained NUTS draws.
#' A lightweight `ctmaFit(fit=FALSE)` setup is sufficient because ctsem's
#' packaged precompiled model is used only for parameter constraints and
#' generated quantities, not for estimation.
#'
#' @param fit A completed `ctmaJuliaFit` object.
#' @param setup A matching fitted or `fit=FALSE` ctsem/CoTiMA object.
#' @param nsamples Number of transformed draws. For Bayesian input, at most
#'   this many retained NUTS draws are used.
#' @param cores Cores passed to ctsem's draw-transformation routine.
#' @param seed Optional reproducibility seed.
#' @param digits Number of decimal places passed to ctsem when constructing
#'   the transformed summary. The default is four.
#' @return The updated `ctmaJuliaFit` object.
#' @export
ctAddJuliaSamples <- function(
    fit, setup, nsamples = 1000L, cores = 2L, seed = NULL, digits = 4L) {
  .ctmaVersionNote()
  .ctma_timer <- .ctmaRegisterExecutionTimer("ctAddJuliaSamples")
  on.exit(.ctma_timer(), add = TRUE)
  pooled <- if (inherits(fit, "CoTiMAPooledJuliaFit")) fit else NULL
  if (!is.null(pooled)) fit <- .ctmaPooledJuliaNative(pooled)
  if (!inherits(fit, "ctmaJuliaFit")) {
    stop("'fit' must inherit from 'ctmaJuliaFit'.", call. = FALSE)
  }
  nsamples <- as.integer(nsamples)[1L]
  cores <- as.integer(cores)[1L]
  digits <- suppressWarnings(as.integer(digits)[1L])
  if (is.na(nsamples) || nsamples < 2L) {
    stop("'nsamples' must be an integer of at least 2.", call. = FALSE)
  }
  if (is.na(cores) || cores < 1L) {
    stop("'cores' must be a positive integer.", call. = FALSE)
  }
  if (is.na(digits) || digits < 1L) {
    stop("'digits' must be a positive integer.", call. = FALSE)
  }
  setup_source <- .ctmaJuliaReadObject(setup)
  bridge_fit <- if (.ctmaJuliaIsUnfittedCoTiMA(setup_source)) {
    .ctmaJuliaGenerateSetup(setup_source)
  } else if (is.list(setup_source) &&
             is.list(setup_source$fit) &&
             !is.null(setup_source$fit$standata) &&
             !is.null(setup_source$fit$setup)) {
    setup_source$fit
  } else if (is.list(setup_source) && !is.null(setup_source$standata) &&
             !is.null(setup_source$setup) &&
             is.null(setup_source$stanfit)) {
    setup_source
  } else {
    .ctmaJuliaExtractFit(setup_source)
  }
  bridge_fit <- .ctmaJuliaAttachConstraintModel(bridge_fit)
  estimates <- as.numeric(coef(fit, type = "raw"))
  if (length(bridge_fit$stanfit$rawest) != length(estimates)) {
    stop(
      sprintf(
        paste0(
          "The ctsem setup has %d raw parameters, but the Julia fit has %d; ",
          "they do not describe the same parameterization."
        ),
        length(bridge_fit$stanfit$rawest), length(estimates)
      ),
      call. = FALSE
    )
  }

  reuse_transformed <- !is.null(fit$ctsem_summary) &&
    !is.null(fit$transformed_draws) &&
    is.matrix(fit$transformed_draws) &&
    nrow(fit$transformed_draws) == nsamples &&
    !is.null(fit$coefficients$transformed) &&
    identical(attr(fit$ctsem_summary, "ctmaJulia_digits"), digits)
  if (reuse_transformed) {
    message("Reusing stored transformed Julia draws.")
    transformed <- list(
      draws = fit$transformed_draws,
      parameters = fit$coefficients$transformed,
      summary = fit$ctsem_summary,
      ctsem_fit = fit$ctsem_fit
    )
  } else if (is.null(fit$bayes)) {
    transformed <- .ctmaJuliaTransformSummary(
      fit = bridge_fit,
      estimates = estimates,
      covariance = vcov(fit),
      loglik = if (is.null(fit$optimization)) {
        fit$likelihood$julia
      } else {
        fit$optimization$loglik
      },
      draws = nsamples,
      cores = cores,
      seed = seed,
      raw_draws = if (!is.null(fit$transformed_draws) &&
          is.matrix(fit$transformed_draws) &&
          nrow(fit$transformed_draws) == nsamples) {
        fit$transformed_draws
      } else {
        NULL
      },
      digits = digits
    )
  } else {
    raw_draws <- do.call(
      rbind,
      lapply(
        seq_len(fit$bayes$controls$chains),
        function(chain) fit$bayes$draws[, , chain]
      )
    )
    if (nrow(raw_draws) > nsamples) {
      keep <- round(seq.int(1L, nrow(raw_draws), length.out = nsamples))
      raw_draws <- raw_draws[keep, , drop = FALSE]
    }
    transformed <- .ctmaJuliaTransformSummary(
      fit = bridge_fit,
      estimates = colMeans(raw_draws),
      covariance = stats::cov(raw_draws),
      loglik = fit$likelihood$julia,
      draws = nrow(raw_draws),
      cores = cores,
      seed = seed,
      raw_draws = raw_draws,
      digits = digits
    )
  }
  cotima <- .ctmaJuliaCoTiMAWEC(
    object = setup_source,
    fit = bridge_fit,
    raw_draws = transformed$draws
  )
  fit$transformed_draws <- transformed$draws
  fit$ctsem_summary <- transformed$summary
  fit$ctsem_fit <- transformed$ctsem_fit
  fit$cotima_summary <- cotima
  fit$coefficients$transformed <- transformed$parameters
  fit$coefficients$population <- transformed$summary$popmeans
  fit$coefficients$tipred <- transformed$summary$tipreds
  fit$coefficients$matrices <- transformed$summary$parmatrices
  fit$coefficients$study <- if (is.null(cotima)) {
    NULL
  } else {
    cotima$study_specific_parameters
  }
  if (is.null(pooled)) fit else .ctmaWrapPooledJuliaFit(fit, pooled)
}

#' @rdname ctAddJuliaSamples
#' @export
ctsmAddJuliaSamples <- ctAddJuliaSamples

#' Transform Julia point estimates without calculating uncertainty
#'
#' Applies ctsem's parameter constraints and generated-quantity transformations
#' to the raw maximum-likelihood estimates stored in a completed
#' `ctmaJuliaFit`.  Unlike [ctAddJuliaSamples()], this function does not require
#' a Hessian or covariance matrix, does not draw samples, and does not repeat
#' optimization.  It is intended for inspecting population, TI-predictor,
#' matrix, and study-specific point estimates after a fast fit performed with
#' `uncertainty = FALSE` and `transformed = FALSE`.
#'
#' The resulting transformed components are point estimates only.  They must
#' not be interpreted as standard errors, confidence intervals, or evidence
#' that the optimizer converged.
#'
#' @param fit A completed `ctmaJuliaFit` object.
#' @param setup A matching fitted, `fit=FALSE`, or prepared ctsem/CoTiMA setup
#'   object supplying ctsem's transformation model.
#' @param cores Number of cores passed to ctsem's transformation code.
#' @param digits Number of decimal places used internally by the ctsem summary.
#' @return The updated `ctmaJuliaFit` object. Its `coefficients$transformed`,
#'   `coefficients$population`, `coefficients$tipred`,
#'   `coefficients$matrices`, and, where available, `coefficients$study`
#'   components contain point transformations. No covariance matrix,
#'   transformed draws, or uncertainty summary is added.
#' @export
ctmaTransformJuliaFit <- function(fit, setup, cores = 1L, digits = 4L) {
  .ctmaVersionNote()
  .ctma_timer <- .ctmaRegisterExecutionTimer("ctmaTransformJuliaFit")
  on.exit(.ctma_timer(), add = TRUE)
  pooled <- if (inherits(fit, "CoTiMAPooledJuliaFit")) fit else NULL
  if (!is.null(pooled)) fit <- .ctmaPooledJuliaNative(pooled)
  if (!inherits(fit, "ctmaJuliaFit")) {
    stop("'fit' must inherit from 'ctmaJuliaFit'.", call. = FALSE)
  }
  cores <- suppressWarnings(as.integer(cores)[1L])
  digits <- suppressWarnings(as.integer(digits)[1L])
  if (is.na(cores) || cores < 1L) {
    stop("'cores' must be a positive integer.", call. = FALSE)
  }
  if (is.na(digits) || digits < 1L) {
    stop("'digits' must be a positive integer.", call. = FALSE)
  }

  setup_source <- .ctmaJuliaReadObject(setup)
  bridge_fit <- if (.ctmaJuliaIsUnfittedCoTiMA(setup_source)) {
    .ctmaJuliaGenerateSetup(setup_source)
  } else if (is.list(setup_source) &&
             is.list(setup_source$fit) &&
             !is.null(setup_source$fit$standata) &&
             !is.null(setup_source$fit$setup)) {
    setup_source$fit
  } else if (is.list(setup_source) && !is.null(setup_source$standata) &&
             !is.null(setup_source$setup) &&
             is.null(setup_source$stanfit)) {
    setup_source
  } else {
    .ctmaJuliaExtractFit(setup_source)
  }
  bridge_fit <- .ctmaJuliaAttachConstraintModel(bridge_fit)
  estimates <- as.numeric(coef(fit, type = "raw"))
  expected <- length(bridge_fit$stanfit$rawest)
  if (expected != length(estimates)) {
    stop(
      sprintf(
        paste0(
          "The ctsem setup has %d raw parameters, but the Julia fit has %d; ",
          "they do not describe the same parameterization."
        ),
        expected, length(estimates)
      ),
      call. = FALSE
    )
  }
  if (any(!is.finite(estimates))) {
    stop("The Julia fit contains non-finite raw estimates.", call. = FALSE)
  }

  # ctsem's transformation routine accepts posterior/sample rows. Repeating
  # the same ML point twice evaluates generated quantities deterministically
  # while deliberately introducing no sampling or Hessian approximation.
  point_rows <- matrix(
    rep(estimates, each = 2L), nrow = 2L, byrow = FALSE
  )
  loglik <- if (!is.null(fit$optimization)) {
    fit$optimization$loglik
  } else {
    fit$likelihood$julia
  }
  transformed <- .ctmaJuliaTransformSummary(
    fit = bridge_fit,
    estimates = estimates,
    covariance = matrix(0, nrow = length(estimates), ncol = length(estimates)),
    loglik = loglik,
    draws = 2L,
    cores = cores,
    raw_draws = point_rows,
    digits = digits
  )
  cotima <- .ctmaJuliaCoTiMAWEC(
    object = setup_source,
    fit = bridge_fit,
    raw_draws = point_rows
  )

  fit$coefficients$transformed <- transformed$parameters
  fit$coefficients$population <- transformed$summary$popmeans
  fit$coefficients$tipred <- transformed$summary$tipreds
  fit$coefficients$matrices <- transformed$summary$parmatrices
  fit$coefficients$study <- if (is.null(cotima)) {
    NULL
  } else {
    cotima$study_specific_parameters
  }
  fit$ctsem_fit <- transformed$ctsem_fit
  fit$point_transformation <- list(
    completed = TRUE,
    uncertainty = FALSE,
    draws = FALSE,
    setup_call = match.call()
  )
  if (is.null(pooled)) fit else .ctmaWrapPooledJuliaFit(fit, pooled)
}

#' Construct or return a ctsem summary for a Julia fit
#'
#' @param object A `ctmaJuliaFit` object.
#' @param ctsem_fit Optional fitted ctsem or CoTiMA object supplying ctsem's
#'   compiled transformation setup. This bridges results originally estimated
#'   from lightweight `ctmaFit(fit=FALSE)` input without repeating Julia
#'   optimization or Hessian calculation.
#' @param transformed_draws Number of normal draws from the stored Julia
#'   estimate and covariance.
#' @param transformed_cores Cores passed to ctsem's transformation code.
#' @param seed Optional random seed for reproducible transformed summaries.
#' @param digits Number of decimal places used when printing the regular
#'   and study-specific summary estimates.
#' @param study_specific Whether to include study-specific parameter summaries
#'   evaluated at each study's observed TI-predictor profile. This supports
#'   dummy-coded and weighted-effect-coded predictors.
#' @param pointEstimate NULL inherits the fit's reporting policy; auto uses optimized mode for ML when available, mean selects earlier draw means.
#' @param ... Reserved for method compatibility.
#' @export
summary.ctmaJuliaFit <- function(
    object, ctsem_fit = NULL, transformed_draws = 1000L,
    transformed_cores = 1L, seed = NULL, digits = 4L,
    study_specific = FALSE, pointEstimate = NULL, ...) {
  digits <- suppressWarnings(as.integer(digits)[1L])
  if (is.na(digits) || digits < 1L) {
    stop("'digits' must be a positive integer.", call. = FALSE)
  }
  if (!is.logical(study_specific) || length(study_specific) != 1L ||
      is.na(study_specific)) {
    stop("'study_specific' must be TRUE or FALSE.", call. = FALSE)
  }
  add_study_specific <- function(summary_object, fit_object) {
    if(is.null(pointEstimate)) pointEstimate <- fit_object$argumentList$pointEstimate
    if(is.null(pointEstimate)) pointEstimate <- "auto"
    summary_object <- .ctmaPointSummary(summary_object,fit_object,pointEstimate)
    if (isTRUE(study_specific)) {
      tables <- fit_object$coefficients$study
      if (is.null(tables)) {
        stop(
          paste0(
            "No study-specific parameters are stored. Supply the matching ",
            "fit=FALSE or fitted CoTiMA object through 'ctsem_fit'."
          ),
          call. = FALSE
        )
      }
      summary_object$study_specific <- tables
    }
    attr(summary_object, "ctmaJulia_digits") <- digits
    attr(summary_object, "cotima_fit_checks") <- .ctmaFitChecks(fit_object)
    class(summary_object) <- unique(c(
      "summary.ctmaJuliaFitStudy",
      class(summary_object)
    ))
    summary_object
  }
  if (is.null(ctsem_fit)) {
    if (!is.null(object$ctsem_summary)) {
      return(add_study_specific(object$ctsem_summary, object))
    }
    stop(
      paste0(
        "No ctsem summary is stored. Supply a fitted ctsem/CoTiMA object via ",
        "'ctsem_fit', or refit with uncertainty=TRUE and transformed=TRUE."
      ),
      call. = FALSE
    )
  }

  updated <- ctAddJuliaSamples(
    fit = object,
    setup = ctsem_fit,
    nsamples = transformed_draws,
    cores = transformed_cores,
    seed = seed,
    digits = digits
  )
  add_study_specific(updated$ctsem_summary, updated)
}

#' Print a ctsem summary followed by study-specific Julia estimates
#'
#' The ordinary ctsem summary is printed unchanged. Study-specific tables are
#' appended only when `summary(..., study_specific = TRUE)` created this
#' summary subclass.
#'
#' @param x A summary returned by `summary.ctmaJuliaFit()` with
#'   `study_specific = TRUE`.
#' @param width Output width passed to ctsem's summary printer.
#' @param row.names Whether to print row names in the study-specific tables.
#' @param digits Number of significant digits used for printed estimates.
#' @param ... Additional arguments passed to the ordinary summary printer.
#' @export
print.summary.ctmaJuliaFitStudy <- function(
    x, width = getOption("width"), row.names = TRUE,
    digits = attr(x, "ctmaJulia_digits"), ...) {
  if (is.null(digits)) digits <- 4L
  digits <- suppressWarnings(as.integer(digits)[1L])
  if (is.na(digits) || digits < 1L) {
    stop("'digits' must be a positive integer.", call. = FALSE)
  }

  format_estimates <- function(value, section = NULL) {
    if (is.matrix(value)) {
      value <- as.data.frame(value, check.names = FALSE)
    }
    if (is.data.frame(value)) {
      identifier_columns <- c("row", "col", "parameter_index")
      for (column in names(value)) {
        if (is.numeric(value[[column]]) &&
            !column %in% identifier_columns) {
          value[[column]] <- formatC(
            value[[column]],
            format = "f",
            digits = digits
          )
        }
      }
    }
    if (is.numeric(value) && length(value) == 1L &&
        !section %in% c("npars", "nsamples")) {
      value <- formatC(value, format = "f", digits = digits)
    }
    value
  }

  study_tables <- x$study_specific

  regular_summary <- x
  regular_summary$study_specific <- NULL
  regular_summary$pointEstimateNote <- NULL
  regular_summary$legacyPointEstimateNote <- NULL
  regular_summary$pointEstimate <- NULL
  class(regular_summary) <- setdiff(
    class(regular_summary),
    "summary.ctmaJuliaFitStudy"
  )
  for (section in names(regular_summary)) {
    regular_summary[[section]] <- format_estimates(
      regular_summary[[section]],
      section = section
    )
  }

  print(
    regular_summary,
    width = width,
    row.names = row.names,
    digits = digits,
    ...
  )

  if (!is.null(study_tables)) {
    cat("\nStudy-specific parameter estimates\n")
    cat("----------------------------------\n")
    for (study_name in names(study_tables)) {
      cat("\n", study_name, "\n", sep = "")
      print(
        format_estimates(study_tables[[study_name]]),
        row.names = row.names,
        right = TRUE
      )
    }
  }

  checks <- attr(x, "cotima_fit_checks")
  if (!is.null(checks)) .ctmaPrintFitChecks(checks)
  .ctmaPrintPointNote(x)
  invisible(x)
}

.ctmaJuliaIsCtsemFit <- function(x) {
  is.list(x) &&
    !is.null(x$standata) &&
    is.list(x$stanfit) &&
    !is.null(x$stanfit$rawest)
}

.ctmaJuliaReadObject <- function(object) {
  if (is.character(object) && length(object) == 1L) {
    if (!file.exists(object)) {
      stop("The supplied RDS file does not exist.", call. = FALSE)
    }
    return(readRDS(object))
  }
  object
}

.ctmaJuliaNamedCI <- function(x, name) {
  if (is.null(x) || is.null(names(x))) {
    return(NULL)
  }
  index <- match(tolower(name), tolower(names(x)))
  if (is.na(index)) NULL else x[[index]]
}

.ctmaJuliaIsUnfittedCoTiMA <- function(object) {
  inherits(object, "CoTiMAFit") &&
    is.list(object) &&
    !is.null(object$data) &&
    !is.null(object$ctModel) &&
    is.null(object$standata) &&
    is.null(object$stanmodeltext) &&
    (is.null(object$studyFitList) || !length(object$studyFitList))
}

.ctmaJuliaSuppressSetupMessage <- function(message) {
  message <- as.character(message)[1L]
  any(vapply(
    c(
      "Maximum likelihood estimation requested",
      "Binary indicators use a linearised (moment-matched Gaussian) measurement update on the stan backend"
    ),
    function(pattern) grepl(pattern, message, fixed = TRUE),
    logical(1)
  ))
}

.ctmaJuliaGenerateSetup <- function(object) {
  if (!requireNamespace("ctsem", quietly = TRUE)) {
    stop("Package 'ctsem' is required for unfitted CoTiMA input.", call. = FALSE)
  }
  stored_args <- .ctmaJuliaNamedCI(
    object$argumentList, "CoTiMAStanctArgs"
  )
  argument <- function(name, default) {
    value <- .ctmaJuliaNamedCI(stored_args, name)
    if (is.null(value)) default else value
  }
  same_initial <- .ctmaJuliaNamedCI(
    object$argumentList, "sameInitialTimes"
  )
  if (is.null(same_initial)) {
    same_initial <- FALSE
  }

  ct_fit <- ctsem::ctFit
  setup_args <- list(
    datalong = object$data,
    model = object$ctModel,
    fit = FALSE,
    optimize = TRUE,
    priors = FALSE,
    intoverstates = argument("intoverstates", TRUE),
    binomial = argument("binomial", FALSE),
    intoverpop = argument("intoverpop", "auto"),
    stationary = argument("stationary", FALSE),
    sameInitialTimes = same_initial,
    savesubjectmatrices = argument("savesubjectmatrices", FALSE),
    verbose = 0
  )
  # The bundled Julia engine reads the Stan parameter layout from this setup.
  if ("backend" %in% names(formals(ct_fit))) setup_args$backend <- "stan"

  message("Generating ctsem setup data.")
  withCallingHandlers(
    do.call(ct_fit, setup_args),
    message = function(condition) {
      if (.ctmaJuliaSuppressSetupMessage(conditionMessage(condition))) {
        invokeRestart("muffleMessage")
      }
    }
  )
}

.ctmaJuliaRawParameterCount <- function(fit) {
  parameter_names <- getFromNamespace("ctFitgetparnamesfromraw", "ctsem")(
    fit
  )
  sum(lengths(parameter_names))
}

.ctmaJuliaDefaultRawStart <- function(fit) {
  nraw <- .ctmaJuliaRawParameterCount(fit)
  start <- numeric(nraw)
  popsetup <- fit$setup$popsetup
  if (is.data.frame(popsetup) &&
      all(c("param", "transform") %in% names(popsetup))) {
    free <- !is.na(popsetup$param) & popsetup$param > 0L
    rows <- popsetup[free, c("param", "transform"), drop = FALSE]
    rows$matrix_name <- sub("[0-9]+$", "", rownames(rows))
    rows <- rows[!duplicated(rows$param), , drop = FALSE]
    # These are neutral raw-scale values under ctsem's standard transforms:
    # stable drift diagonals near -0.3 and covariance diagonals near one.
    # They keep all covariance matrices positive definite without putting the
    # softplus transforms in a saturated region. Unbounded effects and
    # correlation coordinates remain at zero.
    assign_start <- function(matrix_names, value) {
      index <- rows$param[
        rows$transform == 1L & rows$matrix_name %in% matrix_names
      ]
      index <- index[index <= nraw]
      start[index] <<- value
    }
    assign_start("DRIFT", 0.9)
    assign_start(c("DIFFUSION", "MANIFESTVAR"), -1.2)
    assign_start("T0VAR", -0.75)
  }
  start
}

.ctmaJuliaPrepareInput <- function(object, study = 1L, inits = NULL) {
  if (!.ctmaJuliaIsUnfittedCoTiMA(object)) {
    return(list(
      fit = .ctmaJuliaExtractFit(object, study = study),
      setup_only = FALSE,
      starting_values = "retained"
    ))
  }

  fit <- .ctmaJuliaGenerateSetup(object)
  nraw <- .ctmaJuliaRawParameterCount(fit)
  if (is.null(inits)) {
    raw_start <- .ctmaJuliaDefaultRawStart(fit)
    starting_values <- "model-derived"
  } else {
    raw_start <- as.numeric(inits)
    starting_values <- "user-supplied"
  }
  if (length(raw_start) != nraw) {
    stop(
      sprintf("'inits' has length %d; expected %d.", length(raw_start), nraw),
      call. = FALSE
    )
  }
  if (any(!is.finite(raw_start))) {
    stop("'inits' contains non-finite values.", call. = FALSE)
  }
  fit$stanfit <- list(
    rawest = raw_start,
    optimfit = list(value = NULL)
  )
  list(
    fit = fit,
    setup_only = TRUE,
    starting_values = starting_values
  )
}

.ctmaJuliaSessionInitialized <- function() {
  state <- get(".julia", envir = asNamespace("JuliaCall"), inherits = FALSE)
  isTRUE(state$initialized)
}

.ctmaJuliaAttachConstraintModel <- function(fit) {
  if (!requireNamespace("ctsem", quietly = TRUE) ||
      !requireNamespace("rstan", quietly = TRUE)) {
    stop("Packages 'ctsem' and 'rstan' are required to transform draws.",
         call. = FALSE)
  }
  compiled <- fit[["stanmodel", exact = TRUE]]
  if (!inherits(compiled, "stanmodel")) {
    packaged_models <- get("stanmodels", envir = asNamespace("ctsem"))
    generated_data <- !is.null(fit$standata$gendata) &&
      as.integer(fit$standata$gendata)[1L] == 1L
    compiled <- if (generated_data) {
      packaged_models$ctsmgen
    } else {
      packaged_models$ctsm
    }
    fit$stanmodel <- compiled
  }
  if (is.null(fit$stanfit)) {
    nraw <- .ctmaJuliaRawParameterCount(fit)
    fit$stanfit <- list(
      rawest = numeric(nraw),
      cov = matrix(NA_real_, nraw, nraw),
      rawposterior = matrix(numeric(), nrow = 0L, ncol = nraw),
      optimfit = list(value = NULL)
    )
  }
  if (is.null(fit$stanfit$stanfit)) {
    fit$stanfit$stanfit <- getFromNamespace("stan_reinitsf", "ctsem")(
      fit$stanmodel, fit$standata
    )
  }
  if (!inherits(fit, "ctStanFit")) {
    class(fit) <- "ctStanFit"
  }
  fit
}

.ctmaJuliaDimensionsFromFit <- function(fit) {
  scalar_integer <- function(name) {
    value <- fit$standata[[name]]
    if (is.null(value) || !length(value)) {
      stop(sprintf("ctsem standata has no '%s' dimension.", name), call. = FALSE)
    }
    as.integer(value)[1L]
  }
  list(
    ndatapoints = scalar_integer("ndatapoints"),
    nsubjects = scalar_integer("nsubjects"),
    nmanifest = scalar_integer("nmanifest"),
    nlatent = scalar_integer("nlatent"),
    nlatentpop = scalar_integer("nlatentpop"),
    nparams = scalar_integer("nparams"),
    ntipred = scalar_integer("ntipred"),
    ntipredeffects = scalar_integer("ntipredeffects"),
    rawest_length = length(fit$stanfit$rawest)
  )
}

.ctmaJuliaPerformanceAssessment <- function(
    dimensions, julia_threads, chunksize, carefulFit) {
  value <- function(name) as.numeric(dimensions[[name]])[1L]
  nsubjects <- value("nsubjects")
  ndatapoints <- value("ndatapoints")
  nbase <- value("nparams")
  nti <- value("ntipredeffects")
  nraw <- value("rawest_length")
  threads <- max(1, as.integer(julia_threads)[1L])
  chunk <- max(1, as.integer(chunksize)[1L])
  observations_per_subject <- ndatapoints / max(nsubjects, 1)
  local_ad_passes <- ceiling(nbase / chunk)
  subject_work_per_thread <- nsubjects * local_ad_passes / threads

  # This intentionally identifies only the clearest small/moderate cases. It
  # is an advisory, not a fitted runtime model: optimizer conditioning and the
  # number of required iterations can dominate static model dimensions.
  compact_model <- ndatapoints <= 10000 && nraw <= 60
  small_subject_model <- nsubjects <= 1000 && nraw <= 100
  low_parallel_work <- subject_work_per_thread <= 5000 && nraw <= 75
  likely_stan <- compact_model || small_subject_model || low_parallel_work

  classification <- if (likely_stan) {
    "stan_may_be_faster"
  } else {
    "no_clear_advantage"
  }
  message <- if (likely_stan) {
    sprintf(
      paste0(
        "Performance advisory (heuristic): this model has %d subjects, %d ",
        "data rows, %d base parameters, %d TI-effect parameters, and %d raw ",
        "parameters. With %d Julia threads, ctsem's compiled Stan backend may ",
        "complete maximum-likelihood estimation faster than the experimental ",
        "Julia backend. Runtime also depends strongly on conditioning and the ",
        "number of optimizer iterations; this is not a compatibility or ",
        "accuracy warning. To interrupt in RStudio, press Esc or click the ",
        "Stop button. An embedded Julia call may respond only when it reaches ",
        "a safe point; if it remains unresponsive, use Session > Interrupt R, ",
        "and restart the R session only as a last resort. Set ",
        "performance_warning=FALSE to suppress this advisory."
      ),
      as.integer(nsubjects), as.integer(ndatapoints), as.integer(nbase),
      as.integer(nti), as.integer(nraw), threads
    )
  } else {
    paste0(
      "Static model dimensions do not indicate a clear runtime advantage for ",
      "either backend; a short same-start benchmark is more informative."
    )
  }
  list(
    classification = classification,
    message = message,
    metrics = list(
      nsubjects = nsubjects,
      ndatapoints = ndatapoints,
      observations_per_subject = observations_per_subject,
      base_parameters = nbase,
      tipred_effect_parameters = nti,
      raw_parameters = nraw,
      julia_threads = threads,
      local_ad_passes = local_ad_passes,
      subject_work_per_thread = subject_work_per_thread,
      carefulFit = carefulFit
    )
  )
}

.ctmaJuliaTransformSummary <- function(
    fit, estimates, covariance, loglik, draws, cores, seed = NULL,
    raw_draws = NULL, digits = 4L) {
  if (!requireNamespace("ctsem", quietly = TRUE)) {
    stop("Package 'ctsem' is required for transformed summaries.", call. = FALSE)
  }

  old_seed_exists <- exists(
    ".Random.seed", envir = .GlobalEnv, inherits = FALSE
  )
  if (old_seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  if (!is.null(seed)) {
    seed <- as.integer(seed)[1L]
    if (is.na(seed)) {
      stop("'seed' must be coercible to one integer.", call. = FALSE)
    }
    set.seed(seed)
    on.exit({
      if (old_seed_exists) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
  }

  sampled <- !is.null(raw_draws)
  normal_draws <- getFromNamespace("ctOptimNormalDraws", "ctsem")
  update_transformed <- getFromNamespace("ctOptimUpdateTransformed", "ctsem")
  if (is.null(raw_draws)) {
    raw_draws <- normal_draws(estimates, covariance, draws)
  } else {
    raw_draws <- as.matrix(raw_draws)
  }

  transformed_fit <- fit
  transformed_fit$stanfit$rawest <- estimates
  transformed_fit$stanfit$cov <- covariance
  transformed_fit$stanfit$rawposterior <- raw_draws
  transformed_fit$stanfit$samples <- raw_draws
  transformed_fit$stanfit$transformedparsfull <- list(ll = loglik)
  if (!is.null(transformed_fit$stanfit$optimfit)) {
    transformed_fit$stanfit$optimfit$value <- loglik
  }
  transformed_fit <- update_transformed(
    transformed_fit,
    samples = raw_draws,
    cores = cores
  )

  transformed_summary <- summary(
    transformed_fit,
    digits = digits,
    residualcov = !is.null(transformed_fit$stanfit$kalman),
    priorcheck = FALSE
  )
  transformed_summary <- .ctmaPointSummary(transformed_summary,
    transformed_fit,"auto",optimize=!sampled,priors=FALSE)
  attr(transformed_summary, "ctmaJulia_digits") <- as.integer(digits)[1L]
  list(
    draws = raw_draws,
    parameters = transformed_fit$stanfit$transformedpars,
    summary = transformed_summary,
    ctsem_fit = transformed_fit
  )
}

.ctmaJuliaCoTiMAWEC <- function(object, fit, raw_draws) {
  if (!inherits(object, "CoTiMAFit")) {
    return(NULL)
  }

  popsetup <- fit$setup$popsetup
  popvalues <- fit$setup$popvalues
  tie_setup <- as.matrix(fit$standata$TIPREDEFFECTsetup)
  tie_setup[is.na(tie_setup)] <- 0
  tipreds <- as.matrix(fit$standata$tipredsdata)
  if (!ncol(tie_setup) || !ncol(tipreds)) {
    return(NULL)
  }
  nparams <- as.integer(fit$standata$nparams)[1L]
  if (nrow(tie_setup) != nparams || ncol(tie_setup) != ncol(tipreds)) {
    stop("Unexpected ctsem TI-effect setup dimensions.", call. = FALSE)
  }

  study_tipreds <- unique(tipreds)
  # With dummy coding, k studies are represented by k-1 predictors: the
  # reference study is the all-zero row. Weighted-effect coding likewise has
  # one unique predictor profile per study. Unique rows are therefore the
  # authoritative study count for this calculation.
  nstudies <- nrow(study_tipreds)
  reported_nstudies <- NA_integer_
  if (length(object$studyList)) {
    reported_nstudies <- length(object$studyList)
  } else if (!is.null(object$n.studies) && length(object$n.studies)) {
    reported_nstudies <- suppressWarnings(as.integer(object$n.studies)[1L])
  }
  if (!is.na(reported_nstudies) && reported_nstudies > 0L &&
      reported_nstudies != nstudies) {
    stop(
      sprintf(
        "Found %d unique TI-predictor patterns but CoTiMA reports %d studies.",
        nstudies, reported_nstudies
      ),
      call. = FALSE
    )
  }

  # ctsem 4 renamed the shared matrix registry; retain older ctsem support.
  registry_name <- if (exists(".ctMatricesList", asNamespace("ctsem"),
                              inherits = FALSE)) ".ctMatricesList" else "ctStanMatricesList"
  drift_code <- unname(getFromNamespace(registry_name, "ctsem")()$all["DRIFT"])
  if (length(drift_code) != 1L || is.na(drift_code)) {
    stop("ctsem matrix registry does not identify DRIFT.", call. = FALSE)
  }
  drift_rows <- which(
    !is.na(popsetup$matrix) & !is.na(popsetup$copyrow) &
      !is.na(popsetup$param) &
      popsetup$matrix == drift_code & popsetup$copyrow == 0 &
      popsetup$param > 0
  )
  drift_setup <- popsetup[drift_rows, , drop = FALSE]
  original_nlatent <- suppressWarnings(as.integer(object$n.latent)[1L])
  if (!length(original_nlatent) || is.na(original_nlatent) ||
      original_nlatent < 1L) {
    model_nlatent <- bridge_nlatent <- NULL
    if (!is.null(object$ctModel$n.latent)) {
      model_nlatent <- suppressWarnings(as.integer(object$ctModel$n.latent)[1L])
    }
    if (!is.null(fit$standata$nlatent)) {
      bridge_nlatent <- suppressWarnings(as.integer(fit$standata$nlatent)[1L])
    }
    original_nlatent <- if (!is.null(model_nlatent) &&
                            !is.na(model_nlatent) && model_nlatent > 0L) {
      model_nlatent
    } else {
      bridge_nlatent
    }
  }
  original_drift_rows <- which(
    !is.na(drift_setup$row) & !is.na(drift_setup$col) &
      drift_setup$row <= original_nlatent &
      drift_setup$col <= original_nlatent
  )
  drift_setup <- drift_setup[original_drift_rows, , drop = FALSE]
  if (!nrow(drift_setup)) {
    stop(
      paste0(
        "TI-predictor profiles were found, but no valid DRIFT parameter rows ",
        "could be identified in the regenerated ctsem setup."
      ),
      call. = FALSE
    )
  }

  scale_time <- object$argumentList$scaleTime
  if (is.null(scale_time)) {
    scale_time <- 1
  }
  scale_time <- as.numeric(scale_time)[1L]
  tform <- getFromNamespace("tform", "ctsem")
  out <- vector("list", nstudies)
  nonfinite_draws_removed <- 0L

  study_numbers <- seq_len(nstudies)
  if (length(object$studyList) == nstudies) {
    candidate <- vapply(
      object$studyList,
      function(x) {
        value <- x$originalStudyNo
        if (is.null(value) || !length(value)) {
          NA_character_
        } else {
          as.character(value[[1L]])
        }
      },
      character(1L)
    )
    if (!all(is.na(candidate))) {
      study_numbers <- candidate
    }
  }

  for (study_index in seq_len(nstudies)) {
    values <- matrix(
      NA_real_, nrow = nrow(raw_draws), ncol = nrow(drift_setup)
    )
    for (drift_index in seq_len(nrow(drift_setup))) {
      param_index <- as.integer(drift_setup$param[drift_index])
      local_raw <- raw_draws[, param_index]
      for (tipred_index in seq_len(ncol(tie_setup))) {
        effect_index <- tie_setup[param_index, tipred_index]
        if (!is.na(effect_index) && effect_index > 0) {
          local_raw <- local_raw +
            raw_draws[, nparams + effect_index] *
            study_tipreds[study_index, tipred_index]
        }
      }
      setup_row <- match(
        rownames(drift_setup)[drift_index], rownames(popsetup)
      )
      values[, drift_index] <- tform(
        local_raw,
        popsetup$transform[setup_row],
        popvalues$multiplier[setup_row],
        popvalues$meanscale[setup_row],
        popvalues$offset[setup_row],
        popvalues$inneroffset[setup_row],
        fit$setup$extratforms
      ) * scale_time
    }

    statistics <- t(vapply(
      seq_len(ncol(values)),
      function(column) {
        x <- values[, column]
        finite <- is.finite(x)
        if (!any(finite)) {
          stop(
            sprintf(
              paste0(
                "All transformed draws are non-finite for study %d, ",
                "DRIFT parameter '%s'. Inspect its TI-predictor profile and ",
                "raw covariance."
              ),
              study_index, drift_setup$parname[column]
            ),
            call. = FALSE
          )
        }
        nonfinite_draws_removed <<-
          nonfinite_draws_removed + sum(!finite)
        x <- x[finite]
        c(
          Mean = mean(x),
          sd = stats::sd(x),
          `2.5%` = stats::quantile(x, 0.025, names = FALSE),
          `50%` = stats::quantile(x, 0.5, names = FALSE),
          `97.5%` = stats::quantile(x, 0.975, names = FALSE),
          Tvalues = mean(x) / stats::sd(x)
        )
      },
      numeric(6L)
    ))
    table <- data.frame(
      matrix = rep("DRIFT", nrow(drift_setup)),
      row = drift_setup$row,
      col = drift_setup$col,
      statistics,
      check.names = FALSE
    )
    rownames(table) <- make.unique(
      paste("DRIFT", drift_setup$parname), sep = "__cell"
    )
    out[[study_index]] <- table
  }
  names(out) <- paste0("Study ", study_numbers)

  if (nonfinite_draws_removed > 0L) {
    warning(
      sprintf(
        paste0(
          "Removed %d non-finite transformed study-parameter draw(s) before ",
          "calculating summaries. This can occur in the tails of a ",
          "Hessian-normal approximation."
        ),
        nonfinite_draws_removed
      ),
      call. = FALSE
    )
  }

  list(
    study_specific_parameters = out,
    WEC_estimates_original_time_scale = if (
      isTRUE(object$argumentList$WEC)
    ) out else NULL,
    study_tipreds = study_tipreds,
    scaleTime = scale_time,
    nonfinite_draws_removed = nonfinite_draws_removed
  )
}

.ctmaJuliaExtractFit <- function(object, study = 1L) {
  if (.ctmaJuliaIsCtsemFit(object)) {
    return(object)
  }

  study_fit_list <- if (is.list(object)) object$studyFitList else NULL
  if (.ctmaJuliaIsCtsemFit(study_fit_list)) {
    return(study_fit_list)
  }

  study <- as.integer(study)[1L]
  if (is.na(study) || study < 1L) {
    stop("'study' must be a positive integer.", call. = FALSE)
  }
  if (is.list(study_fit_list) && length(study_fit_list) >= study &&
      .ctmaJuliaIsCtsemFit(study_fit_list[[study]])) {
    return(study_fit_list[[study]])
  }

  stop(
    paste0(
      "Could not locate a ctsem fit. Supply a direct ctsem fit or a CoTiMA ",
      "object whose studyFitList is a ctsem fit or contains one."
    ),
    call. = FALSE
  )
}

# Reject future/unknown covariance meanings before loading an estimation engine.
.ctmaJuliaCheckCovarianceCode <- function(standata) {
  code <- standata$choleskymats
  if (!is.numeric(code) || length(code) != 1L || !is.finite(code) ||
      !code %in% c(-1, 0, 1, 2)) stop(
    "Unsupported ctsem covariance transformation code. CoTiMA2's bundled Julia ",
    "engine supports -1/0 (raw correlation), 1 (Cholesky), and 2 (ctsem 4 z). ",
    "Update CoTiMA2 or use backend = 'ctsemjulia'; the model was not fitted.",
    call. = FALSE)
  invisible(code)
}
