.ctmaRedHetWEC <- function(fit, label) {
  means <- fit$modelResults$DRIFTCoeffViaWECMean
  standard_errors <- fit$modelResults$DRIFTCoeffViaWECSD
  if (is.null(means) || is.null(standard_errors)) {
    stop(label, " does not contain reconstructed WEC study effects. ",
      "Refit it with WEC=TRUE and retain transformed model results.",
      call. = FALSE)
  }
  n_studies <- length(fit$studyList)
  n_latent <- as.integer(fit$n.latent)[1L]
  expected <- n_studies * n_latent^2
  means <- as.numeric(unlist(means, use.names = FALSE))
  standard_errors <- as.numeric(unlist(standard_errors, use.names = FALSE))
  if (length(means) != expected || length(standard_errors) != expected) {
    stop(label, " has an incompatible number of reconstructed WEC estimates.",
      call. = FALSE)
  }
  means <- matrix(means, nrow = n_studies, byrow = TRUE)
  standard_errors <- matrix(standard_errors, nrow = n_studies, byrow = TRUE)
  parameter_names <- fit$parameterNames$DRIFT
  if (is.null(parameter_names) || length(parameter_names) != ncol(means)) {
    parameter_names <- names(fit$modelResults$DRIFT)
  }
  if (is.null(parameter_names) || length(parameter_names) != ncol(means)) {
    parameter_names <- paste0("DRIFT_", seq_len(ncol(means)))
  }
  parameter_names <- sub("^DRIFT\\s+", "", as.character(parameter_names))
  colnames(means) <- colnames(standard_errors) <- parameter_names
  list(effects = means, standard_errors = standard_errors,
    parameter_names = parameter_names)
}

.ctmaRedHetReduction <- function(unconditional, moderated) {
  ratio_reduction <- function(before, after) {
    ifelse(is.finite(before) & before > 0 & is.finite(after),
      1 - after / before, NA_real_)
  }
  descriptive_before <- apply(unconditional$effects, 2L, stats::var,
    na.rm = TRUE)
  descriptive_after <- apply(moderated$effects, 2L, stats::var,
    na.rm = TRUE)
  h0 <- unconditional$analysis$heterogeneity
  h1 <- moderated$analysis$heterogeneity
  data.frame(
    parameter = unconditional$parameter_names,
    scale = "transformed DRIFT; two-stage residual heterogeneity",
    unconditional_variance = as.numeric(h0["tau2Drift", ]),
    moderated_variance = as.numeric(h1["tau2Drift", ]),
    proportional_reduction = ratio_reduction(
      as.numeric(h0["tau2Drift", ]), as.numeric(h1["tau2Drift", ])),
    unconditional_descriptive_variance = descriptive_before,
    moderated_descriptive_variance = descriptive_after,
    descriptive_proportional_reduction = ratio_reduction(
      descriptive_before, descriptive_after),
    unconditional_Q = as.numeric(h0["Q_Drift", ]),
    moderated_Q = as.numeric(h1["Q_Drift", ]),
    Q_proportional_reduction = ratio_reduction(
      as.numeric(h0["Q_Drift", ]), as.numeric(h1["Q_Drift", ])),
    unconditional_I2 = as.numeric(h0["I2_Drift", ]),
    moderated_I2 = as.numeric(h1["I2_Drift", ]),
    I2_proportional_reduction = ratio_reduction(
      as.numeric(h0["I2_Drift", ]), as.numeric(h1["I2_Drift", ])),
    row.names = NULL,
    check.names = FALSE
  )
}

.ctmaRedHetFitObject <- function(fit) {
  out <- fit$studyFitList
  if (inherits(out, "ctmaJuliaFit")) out <- out$ctsem_fit
  if (is.null(out$stanfit$rawposterior) &&
      inherits(fit$engine_fit, "ctmaJuliaFit")) {
    out <- fit$engine_fit$ctsem_fit
  }
  out
}

.ctmaRedHetTransformRaw <- function(values, transform) {
  if (is.na(transform) || !nzchar(transform)) return(values)
  vapply(values, function(value) {
    eval(parse(text = transform), envir = list(param = value))
  }, numeric(1L))
}

.ctmaRedHetDriftDraws <- function(fit, undoTimeScaling = TRUE) {
  ctfit <- .ctmaRedHetFitObject(fit)
  posterior <- as.matrix(ctfit$stanfit$rawposterior)
  if (!nrow(posterior)) stop(
    "Discrete-time heterogeneity requires raw draws in both fits. For a ",
    "Julia fit, retain or add transformed samples before calling ctmaRedHet().",
    call. = FALSE)
  map <- ctmaParMap(ctfit)
  drift <- map$rawpopmeans[map$rawpopmeans$matrix == "DRIFT", , drop = FALSE]
  n_latent <- as.integer(fit$n.latent)[1L]
  n_studies <- length(fit$studyList)
  if (nrow(drift) != n_latent^2) stop(
    "The fitted DRIFT parameter map is incompatible with fit$n.latent.",
    call. = FALSE)

  ti <- map$tipredeffectparams
  ti <- ti[ti$rawpopmeans_col %in% drift$rawpopmeans_col, , drop = FALSE]
  ti_names <- ctfit$ctstanmodelbase$TIpredNames
  profiles <- as.data.frame(ctfit$data$tipredsdata)
  if (!length(ti_names) || !nrow(ti)) {
    profiles <- matrix(numeric(), nrow = n_studies, ncol = 0L)
  } else {
    if (is.null(colnames(profiles))) colnames(profiles) <- ti_names
    profiles <- unique(profiles[, ti_names, drop = FALSE])
    if (nrow(profiles) != n_studies) stop(
      "Could not identify exactly one TI-predictor profile per study.",
      call. = FALSE)
  }
  study_dummy_names <- head(ti_names, max(0L, n_studies - 1L))
  scale_time <- if (isTRUE(undoTimeScaling)) {
    as.numeric(fit$argumentList$scaleTime)[1L]
  } else 1
  if (!is.finite(scale_time)) scale_time <- 1

  base <- posterior[, drift$rawposterior_col, drop = FALSE]
  transformed <- function(raw) {
    out <- raw
    for (j in seq_len(ncol(raw))) {
      out[, j] <- .ctmaRedHetTransformRaw(raw[, j], drift$transform[j]) *
        scale_time
    }
    out
  }
  add_ti <- function(raw, study, include) {
    if (!nrow(ti)) return(raw)
    selected <- ti$TI_predictor %in% include
    for (j in which(selected)) {
      profile_value <- profiles[study, ti$TI_predictor[j]]
      if (is.finite(profile_value) && profile_value != 0) {
        target <- match(ti$rawpopmeans_col[j], drift$rawpopmeans_col)
        raw[, target] <- raw[, target] +
          profile_value * posterior[, ti$rawposterior_col[j]]
      }
    }
    raw
  }

  total <- predicted <- vector("list", n_studies)
  moderator_names <- setdiff(ti_names, study_dummy_names)
  for (study in seq_len(n_studies)) {
    total[[study]] <- transformed(add_ti(base, study, ti_names))
    predicted[[study]] <- transformed(add_ti(base, study, moderator_names))
  }
  list(total = total, predicted = predicted, n_latent = n_latent)
}

.ctmaRedHetDiscreteEffects <- function(fit, dt, undoTimeScaling = TRUE) {
  draws <- .ctmaRedHetDriftDraws(fit, undoTimeScaling)
  n_studies <- length(draws$total)
  n_parameters <- draws$n_latent^2
  effect <- standard_error <- matrix(NA_real_, n_studies, n_parameters)
  for (study in seq_len(n_studies)) {
    ndraw <- nrow(draws$total[[study]])
    residual <- matrix(NA_real_, ndraw, n_parameters)
    for (sample in seq_len(ndraw)) {
      total_matrix <- matrix(draws$total[[study]][sample, ],
        draws$n_latent, byrow = TRUE)
      predicted_matrix <- matrix(draws$predicted[[study]][sample, ],
        draws$n_latent, byrow = TRUE)
      residual[sample, ] <- as.numeric(t(
        OpenMx::expm(dt * total_matrix) -
          OpenMx::expm(dt * predicted_matrix)
      ))
    }
    effect[study, ] <- colMeans(residual)
    standard_error[study, ] <- apply(residual, 2L, stats::sd)
  }
  parameter_names <- paste0("dt", dt, "_", fit$parameterNames$DRIFT)
  colnames(effect) <- colnames(standard_error) <- parameter_names
  list(effects = effect, standard_errors = standard_error,
    parameter_names = parameter_names)
}

.ctmaRedHetOneComparison <- function(null, moderated) {
  null$analysis <- .ctmaBiGMetaAnalysis(null$effects, null$standard_errors)
  moderated$analysis <- .ctmaBiGMetaAnalysis(
    moderated$effects, moderated$standard_errors)
  list(null = null, moderated = moderated,
    reduction = .ctmaRedHetReduction(null, moderated))
}

.ctmaRedHetContinuous <- function(ctmaFitObject, ctmaFitObjectMod,
    activeDirectory, digits, call, dt = NULL, undoTimeScaling = TRUE) {
  null <- .ctmaRedHetWEC(ctmaFitObject, "'ctmaFitObject'")
  moderated <- .ctmaRedHetWEC(ctmaFitObjectMod, "'ctmaFitObjectMod'")
  if (!identical(dim(null$effects), dim(moderated$effects)) ||
      !identical(null$parameter_names, moderated$parameter_names)) {
    stop("The two fits must contain the same studies and DRIFT parameters ",
      "in the same order.", call. = FALSE)
  }
  comparison <- .ctmaRedHetOneComparison(null, moderated)
  null <- comparison$null
  moderated <- comparison$moderated
  reduction <- comparison$reduction

  combined <- function(component) {
    out <- cbind(null$analysis[[component]], moderated$analysis[[component]])
    colnames(out) <- c(
      paste0("w/o Mod ", null$parameter_names),
      paste0("with Mod ", null$parameter_names)
    )
    out
  }
  fixed <- combined("common")
  random <- combined("random")
  heterogeneity <- combined("heterogeneity")
  study_null <- cbind(null$effects, null$standard_errors)
  study_moderated <- cbind(moderated$effects, moderated$standard_errors)
  colnames(study_null) <- colnames(study_moderated) <- c(
    null$parameter_names, paste0(null$parameter_names, " SE"))

  legacy_reduction <- rbind(
    tau2Drift = 100 * reduction$proportional_reduction,
    Q_Drift = 100 * reduction$Q_proportional_reduction,
    I2_Drift = 100 * reduction$I2_proportional_reduction
  )
  colnames(legacy_reduction) <- null$parameter_names
  discrete_results <- NULL
  discrete_legacy <- "No dt effects estimated."
  if (!is.null(dt)) {
    discrete_results <- lapply(dt, function(interval) {
      discrete_null <- .ctmaRedHetDiscreteEffects(
        ctmaFitObject, interval, undoTimeScaling)
      discrete_moderated <- .ctmaRedHetDiscreteEffects(
        ctmaFitObjectMod, interval, undoTimeScaling)
      .ctmaRedHetOneComparison(discrete_null, discrete_moderated)
    })
    names(discrete_results) <- paste0("Time Interval = ", dt)
    discrete_legacy <- lapply(discrete_results, function(x) {
      list(
        heterogeneity_reduction = x$reduction,
        unconditional = x$null$analysis,
        moderated = x$moderated$analysis
      )
    })
  }
  legacy_summary <- list(
    continuousTime = list(
      `fixedEffects (w/o Mod | with Mod)` = round(fixed, digits),
      `randomEffects (w/o Mod | with Mod)` = round(random, digits),
      fitObj_DriftAndSE = round(study_null, digits),
      fitObjMod_DriftAndSE = round(study_moderated, digits),
      `heterogeneity (w/o Mod | with Mod)` = round(heterogeneity, digits),
      HeterogeneityReduction = round(legacy_reduction, digits),
      Note = paste(
        "Primary reduction uses sampling-error-adjusted two-stage residual",
        "tau-squared. Negative values denote increased residual heterogeneity."
      )
    ),
    discreteTime = discrete_legacy
  )

  structure(list(
    activeDirectory = activeDirectory,
    plot.type = NULL,
    model.type = "RedHet",
    n.studies = nrow(null$effects),
    n.latent = ctmaFitObject$n.latent,
    studyList = list(ctmaFitObject, ctmaFitObjectMod),
    heterogeneity_reduction = reduction,
    unconditional = null$analysis,
    moderated = moderated$analysis,
    reconstructed_study_effects = list(
      unconditional = null[c("effects", "standard_errors")],
      moderated_residual = moderated[c("effects", "standard_errors")]
    ),
    discrete_time = discrete_results,
    method = paste(
      "Two-stage residual heterogeneity from WEC-reconstructed study DRIFT",
      "effects; the moderated effects exclude substantive moderator",
      "contributions and therefore represent deviations from model-predicted",
      "study means."
    ),
    summary = legacy_summary,
    call = call
  ), class = c("ctmaRedHet", "CoTiMAFit"))
}

#' @export
print.ctmaRedHet <- function(x, digits = 4L, ...) {
  cat("ctmaRedHet: reduction in residual DRIFT heterogeneity\n\n")
  print(x$heterogeneity_reduction, digits = digits, row.names = FALSE)
  cat("\nMethod: ", x$method, "\n", sep = "")
  invisible(x)
}

#' @export
summary.ctmaRedHet <- function(object, ...) object
