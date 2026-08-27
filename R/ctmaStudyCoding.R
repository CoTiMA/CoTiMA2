.ctma_resolve_study_coding <- function(
    studyCoding = c("dummy", "weighted"), studyVarying = NULL,
    WEC = FALSE) {
  if (!is.logical(WEC) || length(WEC) != 1L || is.na(WEC)) {
    stop("'WEC' must be TRUE or FALSE.", call. = FALSE)
  }
  coding <- tolower(as.character(studyCoding)[1L])
  if (is.na(coding) || !coding %in% c("dummy", "weighted")) {
    stop(
      "'studyCoding' must be either 'dummy' or 'weighted'.",
      call. = FALSE
    )
  }

  if (isTRUE(WEC)) {
    warning(
      paste0(
        "'WEC=TRUE' is deprecated. It is now interpreted as ",
        "studyCoding='weighted' with studyVarying=c('CINT', ",
        "'DIFFUSION', 'T0VAR'). Use the new arguments directly; ",
        "the legacy WEC argument will be removed in a future version."
      ),
      call. = FALSE
    )
    coding <- "weighted"
    if (is.null(studyVarying)) {
      studyVarying <- c("CINT", "DIFFUSION", "T0VAR")
    }
  } else if (identical(coding, "weighted") && is.null(studyVarying)) {
    studyVarying <- c("CINT", "DIFFUSION", "T0VAR")
  }

  if (!is.null(studyVarying)) {
    studyVarying <- unique(toupper(as.character(studyVarying)))
    allowed <- c(
      "CINT", "DIFFUSION", "T0VAR", "DRIFT", "T0MEANS",
      "MANIFESTMEANS", "MANIFESTVAR", "LAMBDA"
    )
    invalid <- setdiff(studyVarying, allowed)
    if (length(invalid)) {
      stop(
        sprintf(
          "Unsupported 'studyVarying' value(s): %s. Allowed values are: %s.",
          paste(invalid, collapse = ", "), paste(allowed, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  list(coding = coding, varying = studyVarying)
}

.ctma_apply_study_varying <- function(
    model, study_effect_columns, studyVarying = NULL,
    phantom_intercepts = FALSE) {
  if (is.null(studyVarying) || !length(study_effect_columns)) {
    return(model)
  }
  columns <- intersect(study_effect_columns, colnames(model$pars))
  if (!length(columns)) return(model)

  model$pars[, columns] <- FALSE
  ordinary <- setdiff(studyVarying, "CINT")
  if (length(ordinary)) {
    rows <- model$pars$matrix %in% ordinary
    model$pars[rows, columns] <- TRUE
  }
  if ("CINT" %in% studyVarying) {
    cint_rows <- model$pars$matrix %in% "CINT"
    if (isTRUE(phantom_intercepts)) {
      parameters <- as.character(model$pars$param)
      parameters[is.na(parameters)] <- ""
      cint_rows <- cint_rows |
        (model$pars$matrix %in% "T0MEANS" &
           grepl("_cint($|_)", parameters, ignore.case = TRUE))
    }
    model$pars[cint_rows, columns] <- TRUE
  }

  # Fixed or structurally absent ctsem entries must never receive TI effects.
  missing_parameters <- is.na(model$pars$param)
  model$pars[missing_parameters, columns] <- FALSE
  model
}
