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

.ctmaRequireCtsemJuliaLaplace <- function() {
  status <- .ctmaCtsemJuliaLaplaceStatus()
  if (!isTRUE(status$available)) {
    stop(
      "The experimental ctsem Julia/Laplace adapter is unavailable. ",
      status$reason,
      " Use an isolated ctsem juliaFit installation for prototype testing; ",
      "do not replace a working production ctsem installation.",
      call. = FALSE
    )
  }
  invisible(status)
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

  generated <- .ctma_generate_ctsem_setup(setup, verbose = FALSE)
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

  # Reuse CoTiMA2's established matrix-selection semantics. This also validates
  # matrix dimensions and records exactly which raw parameter rows were chosen.
  generated <- .ctma_generate_ctsem_setup(object, verbose = FALSE)
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
  selected <- unique(selected[nzchar(selected) & !is.na(selected)])
  if (!length(selected)) stop(
    "Select at least one study-level random effect through reDRIFT, ",
    "reDIFFUSION, reT0VAR, or study_parameters.", call. = FALSE
  )
  unknown <- setdiff(selected, unique(stats::na.omit(model$pars$param)))
  if (length(unknown)) stop(
    "Unknown study-level raw parameter(s): ",
    paste(unknown, collapse = ", "), ".", call. = FALSE
  )
  selected_rows <- !is.na(model$pars$param) & model$pars$param %in% selected
  model$pars$indvarying_study[selected_rows] <- TRUE

  if (!is.null(subject_parameters)) {
    subject_parameters <- unique(as.character(subject_parameters))
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

  if (isTRUE(remove_study_fixed_effects) && length(membership$dummy_names)) {
    columns <- paste0(membership$dummy_names, "_effect")
    columns <- intersect(columns, names(model$pars))
    if (length(columns)) model$pars[selected_rows, columns] <- FALSE
  }

  structure(list(
    data = datalong,
    model = model,
    selected_parameters = selected,
    study_labels = membership$labels,
    study_dummy_names = membership$dummy_names,
    ctsem_status = status,
    fit_call = quote(ctsem::ctFit(
      datalong = data, model = model, backend = "julia",
      intoverpop = "laplace", cores = cores
    ))
  ), class = "ctmaCtsemJuliaLaplacePrep")
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
