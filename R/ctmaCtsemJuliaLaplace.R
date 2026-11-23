# Adapter for ctsem's native Julia engine (public backend="ctsemjulia").
# Feature-detected support for the ctsem 3.12 and 4.x development APIs.
# Keep placement optimization separate from posterior sampling diagnostics.

.ctmaCtsemJuliaCorrectionMode <- function(namespace) {
  if (is.null(namespace)) return("missing")
  if (exists("ctLaplaceCorrect", envir = namespace, inherits = FALSE)) {
    return("separate")
  }
  if (!exists("ctLaplaceCheck", envir = namespace, inherits = FALSE)) {
    return("missing")
  }
  check <- get("ctLaplaceCheck", envir = namespace, inherits = FALSE)
  arguments <- names(formals(check))
  if (all(c("correction", "refine") %in% arguments)) {
    return("integrated")
  }
  "missing"
}

.ctmaCtsemJuliaLaplaceStatus <- function() {
  ctsem_available <- requireNamespace("ctsem", quietly = TRUE)
  version <- if (ctsem_available) {
    as.character(utils::packageVersion("ctsem"))
  } else {
    NA_character_
  }
  namespace <- if (ctsem_available) asNamespace("ctsem") else NULL
  correction_mode <- .ctmaCtsemJuliaCorrectionMode(namespace)
  api <- c(
    ctLaplaceCheck = !is.null(namespace) &&
      exists("ctLaplaceCheck", envir = namespace, inherits = FALSE),
    # ctsem <= 3.12 exported ctLaplaceCorrect(). In ctsem 4.0 its correction
    # and optional refinement are integrated into ctLaplaceCheck().
    ctLaplaceCorrect = correction_mode != "missing",
    ctJuliaStatus = !is.null(namespace) &&
      exists("ctJuliaStatus", envir = namespace, inherits = FALSE),
    ctJuliaSetup = !is.null(namespace) &&
      exists("ctJuliaSetup", envir = namespace, inherits = FALSE),
    outerLevelOnly = !is.null(namespace) &&
      exists(".ctAnyVarying", envir = namespace, inherits = FALSE) &&
      exists(".ctOuterVaryingColumns", envir = namespace, inherits = FALSE)
  )
  connector_available <- requireNamespace("JuliaConnectoR", quietly = TRUE)
  available <- ctsem_available && all(api) && connector_available
  reason <- if (!ctsem_available) {
    "Package 'ctsem' is not installed."
  } else if (!all(api)) {
    paste0(
      "The installed ctsem ", version,
      " does not provide the complete native Julia/Laplace API (missing: ",
      paste(names(api)[!api], collapse = ", "), ")."
    )
  } else if (!connector_available) {
    "Package 'JuliaConnectoR' is not installed."
  } else {
    NULL
  }
  structure(list(
    available = available,
    ctsem_version = version,
    JuliaConnectoR = connector_available,
    api = api,
    laplace_correction = correction_mode,
    reason = reason
  ), class = "ctmaCtsemJuliaLaplaceStatus")
}

.ctmaCtsemJuliaInstallCode <- function() {
  paste(
    'ctsemJuliaLib <- path.expand("~/R/ctsem-juliaFit-lib")',
    'dir.create(ctsemJuliaLib, recursive = TRUE, showWarnings = FALSE)',
    '.libPaths(c(ctsemJuliaLib, .libPaths()))',
    'if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")',
    'install.packages("JuliaConnectoR", lib = ctsemJuliaLib)',
    'remotes::install_github(',
    '  "cdriveraus/ctsem", ref = "juliaFit",',
    '  lib = ctsemJuliaLib, dependencies = TRUE, upgrade = "never"',
    ')',
    '',
    '# Restart R, then select the isolated library before loading CoTiMA2:',
    'ctsemJuliaLib <- path.expand("~/R/ctsem-juliaFit-lib")',
    '.libPaths(c(ctsemJuliaLib, .libPaths()))',
    'library(CoTiMA2)',
    sep = "\n"
  )
}

.ctmaRequireCtsemJuliaLaplace <- function() {
  status <- .ctmaCtsemJuliaLaplaceStatus()
  if (!isTRUE(status$available)) {
    stop(
      "The experimental ctsem Julia/Laplace adapter is unavailable. ",
      status$reason,
      "\n\nAt your own risk, install Charles Driver's experimental ctsem ",
      "juliaFit branch into an isolated R library. Do not replace a working ",
      "production ctsem installation:\n\n",
      .ctmaCtsemJuliaInstallCode(),
      call. = FALSE
    )
  }
  invisible(status)
}

.ctmaCheckCtsemJuliaRequest <- function(backend, workflow,
                                        supported = FALSE) {
  requested <- length(backend) && !is.na(backend[1L]) &&
    identical(tolower(as.character(backend[1L])), "ctsemjulia")
  if (!requested) return(invisible(NULL))

  status <- .ctmaRequireCtsemJuliaLaplace()
  if (!isTRUE(supported)) stop(
    "A compatible ctsem Julia/Laplace API is installed, but ", workflow,
    "() does not yet support backend = 'ctsemjulia'. Use ctmaFit() for ",
    "the experimental native ctsem Julia backend.",
    call. = FALSE
  )
  invisible(status)
}

.ctmaCtsemJuliaRandomEffectsRequested <- function(setup, backend_args) {
  selectors <- c(
    "study_random_effects", "subject_random_effects",
    "study_parameters", "subject_parameters"
  )
  explicit <- any(vapply(selectors, function(name) {
    value <- backend_args[[name]]
    !is.null(value) && length(value) && any(nzchar(as.character(value)))
  }, logical(1L)))
  masks <- c("reDRIFT", "reDIFFUSION", "reT0VAR")
  explicit <- explicit || any(vapply(masks, function(name) {
    value <- backend_args[[name]]
    !is.null(value) && length(value) && any(value %in% TRUE, na.rm = TRUE)
  }, logical(1L)))
  if (explicit) return(TRUE)

  pars <- setup$ctModel$pars
  if (!is.data.frame(pars)) return(FALSE)
  varying <- grep(
    "^indvarying($|_)", names(pars), value = TRUE
  )
  any(vapply(varying, function(name) {
    any(pars[[name]] %in% TRUE, na.rm = TRUE)
  }, logical(1L)))
}

.ctmaAnnotateCtsemJuliaResult <- function(result, prepared = NULL,
                                          use_laplace = FALSE) {
  optim <- result$optim
  estimate <- result$estimate
  sampling <- is.list(result$sample)
  flag <- function(x) if (is.logical(x) && length(x) == 1L && !is.na(x)) x else NA
  optimizer_converged <- flag(optim$converged)
  posterior_converged <- flag(result$sample$converged)
  first_non_null <- function(...) {
    values <- list(...)
    for (value in values) if (!is.null(value)) return(value)
    NULL
  }

  # ctsem 3.12.0 separates optimizer diagnostics (`optim`) from estimates and
  # uncertainty (`estimate`).  Retain both native components, but also expose a
  # stable, compact CoTiMA2 view so downstream diagnostics need not depend on
  # the exact ctsem development-branch layout.
  inner_trace <- if (!is.null(optim$trace) &&
      "inner_converged" %in% names(optim$trace)) {
    as.integer(optim$trace$inner_converged)
  } else integer()
  native_inner_flag <- result$laplace$inner_converged
  result$ctsemjulia_diagnostics <- list(
    converged = if (sampling) posterior_converged else
      flag(first_non_null(optim$converged, estimate$converged)),
    convergence_source = if (sampling) "posterior sampling" else "optimization",
    optimizer_converged = optimizer_converged,
    posterior_converged = if (sampling) posterior_converged else NA,
    sampler = if (sampling) result$sample$sampler else NULL,
    gradient_norm = as.numeric(
      first_non_null(optim$gradient_norm, estimate$gradient_norm, NA_real_)
    )[1L],
    iterations = as.integer(
      first_non_null(optim$iterations, estimate$iterations, NA_integer_)
    )[1L],
    optimizer_resolution = list(
      predicted_gain = optim$predicted_gain,
      convergence_tolerance = optim$convergence_tolerance,
      resolution = optim$resolution,
      convergence_floor = optim$convergence_floor
    ),
    likelihood = .ctmaStoredLikelihood(result),
    laplace_accuracy = .ctmaStoredLaplaceAccuracy(result),
    inner_modes = list(
      expected = if (isTRUE(use_laplace) && length(inner_trace)) {
        max(inner_trace, na.rm = TRUE)
      } else NA_integer_,
      converged = if (isTRUE(use_laplace) && length(inner_trace)) {
        as.integer(utils::tail(inner_trace, 1L))
      } else NA_integer_,
      all_converged = if (isTRUE(use_laplace)) {
        if (!is.null(native_inner_flag)) {
          isTRUE(native_inner_flag)
        } else {
          length(inner_trace) > 0L &&
            all(inner_trace == max(inner_trace, na.rm = TRUE))
        }
      } else NA
    )
  )

  subject_loglik <- estimate$subject_loglik
  study_only <- isTRUE(use_laplace) &&
    !is.null(prepared) &&
    !length(prepared$selected_subject_parameters) &&
    length(prepared$selected_parameters)
  repeated <- length(subject_loglik) > 1L &&
    length(unique(as.numeric(subject_loglik))) == 1L
  usable <- length(subject_loglik) > 0L && !(study_only || repeated)
  result$ctsemjulia_diagnostics$casewise_loglik <- list(
    available = length(subject_loglik) > 0L,
    usable_for_loo = usable,
    scope = if (!length(subject_loglik)) {
      "unavailable"
    } else if (study_only) {
      "study-marginal, not personwise"
    } else "personwise",
    repeated_constant = repeated,
    reason = if (study_only) {
      paste0(
        "The model integrates study-level random effects. The native ",
        "subject_loglik vector is not a valid personwise decomposition of ",
        "the marginal likelihood and must not be used for LOO or WAIC."
      )
    } else if (repeated) {
      paste0(
        "All native subject_loglik entries are identical. Treat them as a ",
        "backend reporting placeholder until independently verified; do not ",
        "use them for LOO or WAIC."
      )
    } else NULL
  )
  if (!isTRUE(usable) && length(subject_loglik)) message(
    "Note: ctsemjulia returned subject_loglik values that are not a valid ",
    "personwise decomposition for this fit. Do not use them for LOO or WAIC; ",
    "see fit$ctsemjulia_diagnostics$casewise_loglik."
  )
  if (inherits(result, "ctJuliaFit")) {
    matrices <- if (is.list(result$model_spec) && is.numeric(estimate$raw)) {
      tryCatch(suppressMessages(
        get("ctBackendParMatrices", envir = asNamespace("ctsem"))(result)
      ), error = function(e) NULL)
    } else NULL
    result$ctsemjulia_diagnostics$population_matrices <- matrices
    class(result) <- unique(c("ctmaCtsemJuliaFit", class(result)))
  }
  result
}

.ctmaLaplaceStudyMembership <- function(setup, study_id = NULL,
    study_id_ti = NULL, study_dummy_ti = NULL, verbose = TRUE) {
  if (!is.null(study_id)) {
    values <- as.character(study_id)
    if (length(values) != nrow(setup$data)) stop(
      "study_id must have one value per row of the long data set.",
      call. = FALSE
    )
    subject <- setup$data[[setup$ctModel$subjectIDname]]
    split_study <- split(values, subject)
    if (any(lengths(lapply(split_study, unique)) != 1L)) stop(
      "Each subject must belong to exactly one study.", call. = FALSE
    )
    labels <- unique(values)
    dummy_names <- if (is.null(study_dummy_ti)) character() else {
      as.character(study_dummy_ti)
    }
    missing_dummy <- setdiff(dummy_names, setup$ctModel$TIpredNames)
    if (length(missing_dummy)) stop(
      "Unknown study dummy TI column(s): ",
      paste(missing_dummy, collapse = ", "), ".", call. = FALSE
    )
    return(list(
      row_id = match(values, labels), labels = labels,
      dummy_names = dummy_names
    ))
  }

  model <- setup$ctModel
  data <- setup$data
  ti_names <- model$TIpredNames
  subject <- data[[model$subjectIDname]]
  if (is.null(subject) || anyNA(subject) || !all(ti_names %in% names(data)))
    stop("Native study membership requires subject IDs and the model's TI columns.",
      call. = FALSE)
  first <- !duplicated(subject)
  subject_id <- match(subject, subject[first])
  ti <- data[first, ti_names, drop = FALSE]
  for (name in ti_names) {
    if (!identical(unname(data[[name]]), unname(ti[[name]][subject_id])))
      stop("TI predictor varies within a subject: ", name, call. = FALSE)
  }
  # Reuse the study-contrast inference with subject-level predictors from
  # the native setup rather than relying on ctsem's legacy Stan standata.
  generated <- list(ctstanmodelbase = model,
    standata = list(tipredsdata = as.matrix(ti), nsubjects = nrow(ti)))
  if (!is.null(study_id_ti)) {
    sid <- .ctma_study_id(generated, study_id = NULL,
      study_id_ti = study_id_ti)
    return(list(
      row_id = sid$id[subject_id], labels = sid$labels,
      dummy_names = character()
    ))
  }

  inferred <- .ctma_infer_study_id(
    generated, study_dummy_ti = study_dummy_ti, verbose = verbose
  )
  labels <- unique(as.character(inferred$id))
  list(
    row_id = inferred$id[subject_id], labels = labels,
    dummy_names = inferred$dummy_names
  )
}

.ctmaPrepareCtsemJuliaLaplace <- function(
    object,
    study_id = NULL,
    study_id_ti = NULL,
    study_dummy_ti = NULL,
    reDRIFT = NULL,
    reDIFFUSION = NULL,
    reT0VAR = NULL,
    subject_random_effects = NULL,
    study_random_effects = NULL,
    subject_parameters = NULL,
    study_parameters = NULL,
    remove_study_fixed_effects = TRUE,
    sdscale_study = 1,
    verbose = TRUE) {
  status <- .ctmaRequireCtsemJuliaLaplace()
  if (!.ctma_is_unfitted_cotima(object)) stop(
    "The prototype adapter currently requires the object returned by ",
    "ctmaFit(..., fit = FALSE).", call. = FALSE
  )
  if (!is.numeric(sdscale_study) || length(sdscale_study) != 1L ||
      !is.finite(sdscale_study) || sdscale_study <= 0) stop(
    "sdscale_study must be one finite positive number.", call. = FALSE
  )

  datalong <- object$data
  membership <- list(
    row_id = NULL, labels = character(), dummy_names = character()
  )

  model <- object$ctModel
  subject_name <- model$subjectIDname
  if (is.null(subject_name) || !nzchar(subject_name)) subject_name <- "id"
  model$subjectIDname <- subject_name
  model$groupIDnames <- character()

  if (!"indvarying" %in% names(model$pars)) model$pars$indvarying <- FALSE
  if (!"sdscale" %in% names(model$pars)) model$pars$sdscale <- 1
  model$pars$indvarying_study <- FALSE
  model$pars$sdscale_study <- model$pars$sdscale
  adjustable <- !is.na(model$pars$param)
  model$pars$sdscale_study[adjustable] <- sdscale_study

  determined <- is.list(model$determinedDiffusions)
  fixed_variance <- !is.null(model$determinedDiffusions$stationaryVariance)
  fixed_drift <- model$determinedDiffusions$drift_parameters
  if (fixed_variance && (length(subject_parameters) || length(subject_random_effects)))
    stop("stationaryVariance currently supports study-level DRIFT effects only, not person-level effects.", call. = FALSE)

  parameters_in_matrices <- function(matrices, level) {
    if (is.null(matrices)) return(character())
    matrices <- unique(toupper(as.character(matrices)))
    if (determined && identical(level, "study")) {
      if (any(!matrices %in% "DRIFT")) stop(
        "With determinedDiffusions, only DRIFT may have study random effects; DIFFUSION and T0VAR are derived.", call. = FALSE)
      return(fixed_drift)
    }

    if (determined && identical(level, "subject") && any(matrices %in% c("DRIFT", "DIFFUSION", "T0VAR", "PARS")))
      stop("Determined-diffusion process random effects are supported at study level only.", call. = FALSE)
    matrices <- matrices[nzchar(matrices) & !is.na(matrices)]
    available <- unique(toupper(as.character(stats::na.omit(
      model$pars$matrix
    ))))
    unknown_matrices <- setdiff(matrices, available)
    if (length(unknown_matrices)) stop(
      "Unknown ", level, "-level random-effect matrix name(s): ",
      paste(unknown_matrices, collapse = ", "), ". Available matrix names: ",
      paste(sort(available), collapse = ", "), ".", call. = FALSE
    )
    rows <- !is.na(model$pars$param) &
      toupper(as.character(model$pars$matrix)) %in% matrices
    parameters <- unique(as.character(model$pars$param[rows]))
    if (!length(parameters)) stop(
      "No free raw parameters were found in the requested ", level,
      "-level random-effect matrices.", call. = FALSE
    )
    parameters
  }

  # Reuse CoTiMA2's established matrix-selection semantics. This also validates
  # matrix dimensions and records exactly which raw parameter rows were chosen.
  generated <- suppressMessages(
    .ctma_generate_ctsem_setup(object, verbose = FALSE)
  )
  selected <- character()
  if (determined) {
    if (any(as.logical(c(reDIFFUSION, reT0VAR)), na.rm = TRUE)) stop(
      "With determinedDiffusions, DIFFUSION and T0VAR cannot have separate random effects.", call. = FALSE)
    if (!is.null(reDRIFT)) {
      if (length(reDRIFT) != 1L || is.na(reDRIFT) || !is.logical(reDRIFT) && !is.numeric(reDRIFT))
        stop("For determined-diffusion models use reDRIFT=TRUE or study_random_effects='DRIFT'.", call. = FALSE)
      if (isTRUE(as.logical(reDRIFT))) selected <- fixed_drift
    }
    spec <- NULL
  } else if (!is.null(reDRIFT) || !is.null(reDIFFUSION) || !is.null(reT0VAR)) {
    spec <- .ctma_re_spec(generated, reDRIFT, reDIFFUSION, reT0VAR)
    selected <- spec$parameter[spec$type == 1L]
    if (any(spec$type != 1L)) stop(
      "The first native-Laplace prototype supports raw parameter random ",
      "effects only; covariance-cell T0VAR effects are not yet translated.",
      call. = FALSE
    )
  } else {
    spec <- NULL
  }
  if (!is.null(study_parameters)) {
    selected <- unique(c(selected, as.character(study_parameters)))
  }
  selected <- unique(c(
    selected,
    parameters_in_matrices(study_random_effects, "study")
  ))
  selected <- unique(selected[nzchar(selected) & !is.na(selected)])
  unknown <- setdiff(selected, unique(stats::na.omit(model$pars$param)))
  if (length(unknown)) stop(
    "Unknown study-level raw parameter(s): ",
    paste(unknown, collapse = ", "), ".", call. = FALSE
  )
  if (determined && any(!selected %in% fixed_drift))
    stop("Only original DRIFT coordinates may have study random effects with determinedDiffusions.", call. = FALSE)
  selected_rows <- !is.na(model$pars$param) & model$pars$param %in% selected
  model$pars$indvarying_study[selected_rows] <- TRUE

  subject_parameters <- unique(c(
    as.character(subject_parameters),
    parameters_in_matrices(subject_random_effects, "subject")
  ))
  subject_parameters <- subject_parameters[
    nzchar(subject_parameters) & !is.na(subject_parameters)
  ]
  if (determined && any(subject_parameters %in% model$determinedDiffusions$auxiliary))
    stop("Determined-diffusion auxiliary coordinates cannot have person-level random effects.", call. = FALSE)
  if (length(subject_parameters)) {
    unknown_subject <- setdiff(
      subject_parameters, unique(stats::na.omit(model$pars$param))
    )
    if (length(unknown_subject)) stop(
      "Unknown subject-level raw parameter(s): ",
      paste(unknown_subject, collapse = ", "), ".", call. = FALSE
    )
    subject_rows <- !is.na(model$pars$param) &
      model$pars$param %in% subject_parameters
    model$pars$indvarying[subject_rows] <- TRUE
  }

  has_subject_effects <- any(model$pars$indvarying %in% TRUE, na.rm = TRUE)
  has_study_effects <- any(
    model$pars$indvarying_study %in% TRUE, na.rm = TRUE
  )
  if (!has_subject_effects && !has_study_effects) stop(
    "Select at least one subject- or study-level random effect.",
    call. = FALSE
  )
  selected_subject_rows <- model$pars$indvarying %in% TRUE &
    !is.na(model$pars$param)
  selected_subject_parameters <- unique(as.character(
    model$pars$param[selected_subject_rows]
  ))

  if (has_study_effects) {
    membership <- .ctmaLaplaceStudyMembership(
      object, study_id = study_id, study_id_ti = study_id_ti,
      study_dummy_ti = study_dummy_ti, verbose = verbose
    )
    datalong$study <- membership$row_id
    model$groupIDnames <- "study"
  }

  if (isTRUE(remove_study_fixed_effects) && length(membership$dummy_names)) {
    columns <- paste0(membership$dummy_names, "_effect")
    columns <- intersect(columns, names(model$pars))
    if (length(columns)) model$pars[selected_rows, columns] <- FALSE
  }

  structure(list(
    data = datalong,
    model = model,
    selected_parameters = selected,
    selected_subject_parameters = selected_subject_parameters,
    study_labels = membership$labels,
    study_dummy_names = membership$dummy_names,
    ctsem_status = status,
    fit_call = quote(ctsem::ctFit(
      datalong = data, model = model, backend = "julia",
      intoverpop = "laplace", cores = cores
    ))
  ), class = "ctmaCtsemJuliaLaplacePrep")
}

.ctmaFitCtsemJuliaBackend <- function(
    setup, backend_args, cores, cotima_call = NULL,
    inits = NULL, priors = FALSE, verbose = 0, estimation = "ML", compute = NULL) {
  prepare_names <- c(
    "study_id", "study_id_ti", "study_dummy_ti", "reDRIFT",
    "reDIFFUSION", "reT0VAR", "subject_random_effects",
    "study_random_effects", "subject_parameters", "study_parameters",
    "remove_study_fixed_effects", "sdscale_study"
  )
  fit_args <- list()
  fit_args$inits <- if (!is.null(backend_args$inits)) {
    backend_args$inits
  } else inits
  fit_args$priors <- if (!is.null(backend_args$priors)) {
    backend_args$priors
  } else priors
  fit_args$verbose <- if (!is.null(backend_args$verbose)) {
    backend_args$verbose
  } else verbose
  fit_args$stationary <- backend_args$stationary
  fit_args$optimize <- identical(tolower(estimation), "ml")
  if (!is.null(backend_args$optimcontrol)) {
    if (!is.list(backend_args$optimcontrol)) stop(
      "backendArgs$optimcontrol must be NULL or a named list.",
      call. = FALSE
    )
    fit_args$optimcontrol <- backend_args$optimcontrol
  }
  if (!is.null(backend_args$sampleControl)) {
    if (!is.list(backend_args$sampleControl)) stop(
      "backendArgs$sampleControl must be NULL or a named list.",
      call. = FALSE
    )
    fit_args$sampleControl <- backend_args$sampleControl
  }
  fit_args <- fit_args[!vapply(fit_args, is.null, logical(1L))]

  use_laplace <- .ctmaCtsemJuliaRandomEffectsRequested(
    setup, backend_args
  )
  fit_args$intoverpop <- .ctmaCtsemJuliaRoute(
    backend_args, estimation, use_laplace
  )
  if (identical(compute$tasks, "likelihood")) fit_args$compute <- compute
  if (use_laplace) {
    prepare_args <- backend_args[
      intersect(prepare_names, names(backend_args))
    ]
    prepared <- do.call(
      .ctmaPrepareCtsemJuliaLaplace,
      c(list(object = setup), prepare_args)
    )
    result <- do.call(
      .ctmaFitCtsemJuliaLaplace,
      c(list(prepared = prepared, cores = cores), fit_args)
    )
  } else {
    prepared <- NULL
    result <- do.call(
      .ctmaFitCtsemJuliaDirect,
      c(list(setup = setup, cores = cores), fit_args)
    )
  }
  result$source_object <- setup
  result$ctsemjulia_prepared <- prepared
  result$ctsemjulia_mode <- .ctmaCtsemJuliaMode(
    result, fit_args$intoverpop, !fit_args$optimize, use_laplace
  )
  result$cotima_call <- cotima_call
  result$backend <- "ctsemjulia"
  result$argumentList <- setup$argumentList
  result$argumentList$backend <- "ctsemjulia"
  result$argumentList$backendArgs <- backend_args
  result$argumentList$fit <- TRUE
  .ctmaAnnotateCtsemJuliaResult(
    result, prepared = prepared,
    use_laplace = result$ctsemjulia_mode %in% c("laplace", "quadrature")
  )
}

.ctmaFitCtsemJuliaDirect <- function(
    setup, cores = 1L, intoverpop = "auto", compute = NULL, ...) {
  .ctmaRequireCtsemJuliaLaplace()
  if (!.ctma_is_unfitted_cotima(setup)) stop(
    "setup must be returned by ctmaFit(..., fit = FALSE).",
    call. = FALSE
  )
  cores <- as.integer(cores)[1L]
  if (is.na(cores) || cores < 1L) stop(
    "cores must be one positive integer.", call. = FALSE
  )
  args <- c(list(
    datalong = setup$data, model = setup$ctModel, backend = "julia",
    intoverpop = intoverpop, cores = cores), list(...))
  if (!is.null(compute)) return(.ctmaNativeLikelihoodFit(args, compute))
  .ctmaNativeCtsemJuliaFit(args)
}

.ctmaFitCtsemJuliaLaplace <- function(
    prepared, cores = 1L, intoverpop = "laplace", compute = NULL, ...) {
  .ctmaRequireCtsemJuliaLaplace()
  if (!inherits(prepared, "ctmaCtsemJuliaLaplacePrep")) stop(
    "prepared must be returned by .ctmaPrepareCtsemJuliaLaplace().",
    call. = FALSE
  )
  cores <- as.integer(cores)[1L]
  if (is.na(cores) || cores < 1L) stop(
    "cores must be one positive integer.", call. = FALSE
  )
  args <- c(list(
    datalong = prepared$data, model = prepared$model, backend = "julia",
    intoverpop = intoverpop, cores = cores), list(...))
  if (!is.null(compute)) return(.ctmaNativeLikelihoodFit(args, compute))
  .ctmaNativeCtsemJuliaFit(args)
}
