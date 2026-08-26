# Internal planning and validation for Julia's subject-subset warm start.

.ctmaWarmStartMode <- function(value) {
  if (is.logical(value) && length(value) == 1L && !is.na(value)) {
    return(if (isTRUE(value)) "always" else "never")
  }
  if (is.character(value) && length(value) == 1L &&
      tolower(value) %in% c("auto", "always", "never")) {
    return(tolower(value))
  }
  stop("'warm_start' must be 'auto', TRUE, or FALSE.", call. = FALSE)
}

.ctmaWarmStartStudy <- function(standata, source_object, nsubjects) {
  nstudies <- if (is.list(source_object)) {
    length(source_object$studyList)
  } else 0L
  if (!is.finite(nstudies) || nstudies < 2L) return(rep.int(1L, nsubjects))
  tipreds <- as.matrix(standata$tipredsdata)
  if (nrow(tipreds) != nsubjects || ncol(tipreds) < nstudies - 1L) {
    return(rep.int(1L, nsubjects))
  }
  study_design <- tipreds[, seq_len(nstudies - 1L), drop = FALSE]
  keys <- apply(round(study_design, 12L), 1L, paste, collapse = "\r")
  match(keys, unique(keys))
}

.ctmaWarmStartPlan <- function(
    fit, source_object, mode = "auto", fraction = 0.10,
    subjects_per_parameter = 20L, minimum_per_study = 25L,
    maximum_subjects = 5000L, binary_minimum = 20L, seed = NULL) {
  mode <- .ctmaWarmStartMode(mode)
  nsubjects <- as.integer(fit$standata$nsubjects)[1L]
  nbase <- as.integer(fit$standata$nparams)[1L]
  controls <- c(
    fraction = fraction,
    subjects_per_parameter = subjects_per_parameter,
    minimum_per_study = minimum_per_study,
    maximum_subjects = maximum_subjects,
    binary_minimum = binary_minimum
  )
  if (any(!is.finite(controls)) || fraction <= 0 || fraction > 1 ||
      any(controls[-1L] < 1)) {
    stop(
      "Warm-start fraction must be in (0, 1], and its count controls must be ",
      "positive.", call. = FALSE
    )
  }
  counts <- as.integer(round(controls[-1L]))
  names(counts) <- names(controls)[-1L]
  if (mode == "never") return(list(
    use = FALSE, reason = "disabled", subjects = integer(), nsubjects = nsubjects
  ))
  if (!is.finite(nsubjects) || nsubjects < 2L || !is.finite(nbase) ||
      nbase < 1L) return(list(
    use = FALSE, reason = "invalid model dimensions", subjects = integer(),
    nsubjects = nsubjects
  ))

  subject_rows <- tabulate(
    as.integer(fit$standata$subject), nbins = nsubjects
  )
  argument_list <- if (is.list(source_object)) {
    source_object$argumentList
  } else NULL
  pooled_ri <- identical(argument_list$randomInterceptStructure, "pooled")
  minimum_waves <- if (pooled_ri) 3L else 2L
  informative <- which(subject_rows >= minimum_waves)
  required_informative <- min(
    nsubjects, counts[["subjects_per_parameter"]] * nbase
  )
  if (length(informative) < min(required_informative, nsubjects)) {
    return(list(
      use = FALSE,
      reason = sprintf(
        "only %d subject(s) have at least %d observations; %d are required",
        length(informative), minimum_waves,
        min(required_informative, nsubjects)
      ),
      subjects = integer(), nsubjects = nsubjects
    ))
  }

  strata <- .ctmaWarmStartStudy(fit$standata, source_object, nsubjects)
  groups <- split(seq_len(nsubjects), strata)
  per_study_required <- sum(vapply(
    groups, function(index) min(length(index), counts[["minimum_per_study"]]),
    integer(1L)
  ))
  target <- max(
    ceiling(fraction * nsubjects),
    required_informative,
    per_study_required
  )
  # maximum_subjects is a soft computational ceiling: structural minima win.
  target <- max(
    required_informative, per_study_required,
    min(target, counts[["maximum_subjects"]])
  )
  target <- min(target, nsubjects)
  if (mode == "auto" && target >= ceiling(0.8 * nsubjects)) {
    return(list(
      use = FALSE,
      reason = "the defensible pilot would use at least 80% of all subjects",
      subjects = integer(), nsubjects = nsubjects, target = target
    ))
  }

  if (!is.null(seed)) {
    seed <- suppressWarnings(as.integer(seed)[1L])
    if (is.na(seed)) stop("'seed' must be coercible to an integer.", call. = FALSE)
  } else seed <- 1L
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  safe_sample <- function(x, size) {
    if (size <= 0L) return(x[0L])
    if (size >= length(x)) return(x)
    sample(x, size)
  }

  selected <- unique(unlist(lapply(groups, function(index) {
    preferred <- intersect(index, informative)
    pool <- c(preferred, setdiff(index, preferred))
    take <- min(length(pool), counts[["minimum_per_study"]])
    if (take) safe_sample(pool, take) else integer()
  }), use.names = FALSE))
  remaining <- setdiff(seq_len(nsubjects), selected)
  preferred_remaining <- c(
    intersect(remaining, informative), setdiff(remaining, informative)
  )
  if (length(selected) < target) {
    selected <- c(
      selected,
      safe_sample(preferred_remaining, target - length(selected))
    )
  }

  # Ensure adequate representation of both categories for every binary
  # manifest. Add complete subjects; never sample individual rows.
  manifesttype <- as.integer(fit$standata$manifesttype)
  binary <- which(manifesttype == 1L)
  Y <- as.matrix(fit$standata$Y)
  row_subject <- as.integer(fit$standata$subject)
  binary_requirements <- list()
  for (manifest in binary) for (value in 0:1) {
    supporting <- unique(row_subject[
      is.finite(Y[, manifest]) & Y[, manifest] == value
    ])
    required <- min(length(supporting), counts[["binary_minimum"]])
    present <- intersect(selected, supporting)
    if (length(present) < required) {
      add <- setdiff(supporting, selected)
      selected <- c(selected, safe_sample(add, required - length(present)))
    }
    binary_requirements[[paste(manifest, value, sep = ":")]] <- c(
      available = length(supporting), selected = sum(selected %in% supporting),
      required = required
    )
  }
  selected <- sort(unique(as.integer(selected)))
  selected_informative <- sum(selected %in% informative)
  if (selected_informative < required_informative) {
    add <- setdiff(informative, selected)
    selected <- sort(c(
      selected, safe_sample(add, required_informative - selected_informative)
    ))
  }
  if (length(selected) >= nsubjects) return(list(
    use = FALSE, reason = "structural checks require the full dataset",
    subjects = integer(), nsubjects = nsubjects, target = length(selected)
  ))
  list(
    use = TRUE, reason = "eligible", subjects = selected,
    nsubjects = nsubjects, selected = length(selected),
    fraction = length(selected) / nsubjects, requested_fraction = fraction,
    nbase = nbase, minimum_waves = minimum_waves,
    informative_selected = sum(selected %in% informative),
    required_informative = required_informative,
    studies = length(groups), study_counts = table(strata[selected]),
    binary_requirements = binary_requirements, seed = seed
  )
}
