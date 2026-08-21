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

.ctmaRedHetContinuous <- function(ctmaFitObject, ctmaFitObjectMod,
    activeDirectory, digits, call) {
  null <- .ctmaRedHetWEC(ctmaFitObject, "'ctmaFitObject'")
  moderated <- .ctmaRedHetWEC(ctmaFitObjectMod, "'ctmaFitObjectMod'")
  if (!identical(dim(null$effects), dim(moderated$effects)) ||
      !identical(null$parameter_names, moderated$parameter_names)) {
    stop("The two fits must contain the same studies and DRIFT parameters ",
      "in the same order.", call. = FALSE)
  }
  null$analysis <- .ctmaBiGMetaAnalysis(null$effects, null$standard_errors)
  moderated$analysis <- .ctmaBiGMetaAnalysis(
    moderated$effects, moderated$standard_errors)
  reduction <- .ctmaRedHetReduction(null, moderated)

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
    discreteTime = "No dt effects estimated."
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
