# Internal prototype adapter for ctsem's native Julia/Laplace backend.
#
# This file deliberately exports no user-facing function and does not alter
# CoTiMA2's existing backend dispatch. It isolates the compatibility checks and
# model translation needed before a future backend = "ctsemjulia" can be made
# public. The implementation follows Charles Driver's ctsem 3.12.0 juliaFit
# API; feature detection is used because that API is not yet on CRAN.

.ctmaCtsemJuliaLaplaceStatus <- function() {
  ctsem_available <- requireNamespace("ctsem", quietly = TRUE)
  version <- if (ctsem_available) {
    as.character(utils::packageVersion("ctsem"))
  } else {
    NA_character_
  }
  namespace <- if (ctsem_available) asNamespace("ctsem") else NULL
  api <- c(
    ctLaplaceCheck = !is.null(namespace) &&
      exists("ctLaplaceCheck", envir = namespace, inherits = FALSE),
    ctLaplaceCorrect = !is.null(namespace) &&
      exists("ctLaplaceCorrect", envir = namespace, inherits = FALSE),
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

  generated <- suppressMessages(
    .ctma_generate_ctsem_setup(setup, verbose = FALSE)
  )
  if (!is.null(study_id_ti)) {
    sid <- .ctma_study_id(generated, study_id = NULL,
      study_id_ti = study_id_ti)
    subject_id <- generated$standata$subject
    return(list(
      row_id = sid$id[subject_id], labels = sid$labels,
      dummy_names = character()
    ))
  }

  inferred <- .ctma_infer_study_id(
    generated, study_dummy_ti = study_dummy_ti, verbose = verbose
  )
  subject_id <- generated$standata$subject
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

  membership <- .ctmaLaplaceStudyMembership(
    object, study_id = study_id, study_id_ti = study_id_ti,
    study_dummy_ti = study_dummy_ti, verbose = verbose
  )
  datalong <- object$data
  datalong$study <- membership$row_id

  model <- object$ctModel
  subject_name <- model$subjectIDname
  if (is.null(subject_name) || !nzchar(subject_name)) subject_name <- "id"
  model$subjectIDname <- subject_name
  model$groupIDnames <- "study"

  if (!"indvarying" %in% names(model$pars)) model$pars$indvarying <- FALSE
  if (!"sdscale" %in% names(model$pars)) model$pars$sdscale <- 1
  model$pars$indvarying_study <- FALSE
  model$pars$sdscale_study <- model$pars$sdscale
  adjustable <- !is.na(model$pars$param)
  model$pars$sdscale_study[adjustable] <- sdscale_study

  parameters_in_matrices <- function(matrices, level) {
    if (is.null(matrices)) return(character())
    matrices <- unique(toupper(as.character(matrices)))
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
  if (!is.null(reDRIFT) || !is.null(reDIFFUSION) || !is.null(reT0VAR)) {
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
  selected_rows <- !is.na(model$pars$param) & model$pars$param %in% selected
  model$pars$indvarying_study[selected_rows] <- TRUE

  subject_parameters <- unique(c(
    as.character(subject_parameters),
    parameters_in_matrices(subject_random_effects, "subject")
  ))
  subject_parameters <- subject_parameters[
    nzchar(subject_parameters) & !is.na(subject_parameters)
  ]
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
    inits = NULL, priors = FALSE, verbose = 0) {
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
  if (!is.null(backend_args$optimcontrol)) {
    if (!is.list(backend_args$optimcontrol)) stop(
      "backendArgs$optimcontrol must be NULL or a named list.",
      call. = FALSE
    )
    fit_args$optimcontrol <- backend_args$optimcontrol
  }
  fit_args <- fit_args[!vapply(fit_args, is.null, logical(1L))]

  use_laplace <- .ctmaCtsemJuliaRandomEffectsRequested(
    setup, backend_args
  )
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
  result$ctsemjulia_mode <- if (use_laplace) "laplace" else "direct"
  result$cotima_call <- cotima_call
  result$backend <- "ctsemjulia"
  result$argumentList <- setup$argumentList
  result$argumentList$backend <- "ctsemjulia"
  result$argumentList$backendArgs <- backend_args
  result$argumentList$fit <- TRUE
  result
}

.ctmaFitCtsemJuliaDirect <- function(setup, cores = 1L, ...) {
  .ctmaRequireCtsemJuliaLaplace()
  if (!.ctma_is_unfitted_cotima(setup)) stop(
    "setup must be returned by ctmaFit(..., fit = FALSE).",
    call. = FALSE
  )
  cores <- as.integer(cores)[1L]
  if (is.na(cores) || cores < 1L) stop(
    "cores must be one positive integer.", call. = FALSE
  )
  ctsem::ctFit(
    datalong = setup$data, model = setup$ctModel, backend = "julia",
    intoverpop = "auto", cores = cores, ...
  )
}

.ctmaFitCtsemJuliaLaplace <- function(prepared, cores = 1L, ...) {
  .ctmaRequireCtsemJuliaLaplace()
  if (!inherits(prepared, "ctmaCtsemJuliaLaplacePrep")) stop(
    "prepared must be returned by .ctmaPrepareCtsemJuliaLaplace().",
    call. = FALSE
  )
  cores <- as.integer(cores)[1L]
  if (is.na(cores) || cores < 1L) stop(
    "cores must be one positive integer.", call. = FALSE
  )
  ctsem::ctFit(
    datalong = prepared$data,
    model = prepared$model,
    backend = "julia",
    intoverpop = "laplace",
    cores = cores,
    ...
  )
}
