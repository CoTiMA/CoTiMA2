.ctmaStoredLikelihood <- function(object) {
  # Preserve source definitions: checked likelihood is not the optimized
  # Laplace objective, and a missing method must not be guessed.
  estimate <- object$estimate
  scalar <- function(x) if (is.numeric(x) && length(x) == 1L &&
    is.finite(x)) as.numeric(x) else NA_real_
  method <- estimate$loglik_method
  if (!is.character(method) || length(method) != 1L || is.na(method))
    method <- "not recorded"
  list(method = method, loglik = scalar(estimate$loglik),
    loglik_laplace = scalar(estimate$loglik_laplace),
    loglik_se = scalar(estimate$loglik_se))
}

.ctmaStoredLaplaceAccuracy <- function(object) {
  # Read cached metadata only: never call ctLaplaceCheck or initialize Julia.
  laplace <- object$laplace
  applicable <- is.list(laplace) && length(laplace) > 0L
  if (!is.list(laplace)) laplace <- list()
  scalar <- function(x) if (is.numeric(x) && length(x) == 1L &&
    is.finite(x)) as.numeric(x) else NA_real_
  text <- function(x) if (is.character(x) && length(x) == 1L &&
    !is.na(x) && nzchar(x)) x else NULL
  correction <- if (is.list(laplace$correction)) laplace$correction else list()
  check <- if (is.list(correction$check)) correction$check else list()
  status <- text(check$status)
  # "checked" records an operation, not a verdict about its accuracy.
  if (is.null(status)) status <- text(correction$status)
  if (is.null(status)) status <- "not assessed in stored diagnostics"
  applied <- if (is.logical(correction$applied) &&
    length(correction$applied) == 1L && !is.na(correction$applied))
      correction$applied else NA
  list(applicable = applicable, status = status,
    check_status = text(check$status), correction_status = text(correction$status),
    correction_applied = applied, gap = scalar(correction$gap),
    maximum_shift = scalar(check$max_shift),
    shift_noise = scalar(check$noise), shift_threshold = scalar(check$bar))
}

.ctmaLaplaceAccuracyNote <- function(accuracy) {
  if (!isTRUE(accuracy$applicable)) return(character())
  note <- paste0("Laplace approximation accuracy: ", accuracy$status,
    " (stored metadata).")
  if (identical(tolower(accuracy$status), "checked")) note <- paste(note,
    "A completed check is not an accuracy verdict; no verdict is stored.")
  if (is.finite(accuracy$gap)) note <- paste0(note,
    " Reported integration gap = ", format(signif(accuracy$gap, 4)), ".")
  if (is.finite(accuracy$maximum_shift)) note <- paste0(note,
    " Reported maximum standardized shift = ",
    format(signif(accuracy$maximum_shift, 4)), ".")
  if (!is.na(accuracy$correction_applied)) note <- paste0(note,
    " Correction applied: ", if (accuracy$correction_applied) "yes" else "no", ".")
  paste(note, "Optimizer convergence does not certify integration accuracy.")
}

.ctmaIsNativeFit <- function(object) {
  is.list(object) && inherits(object, c("ctJuliaFit", "ctmaJuliaFit", "ctStanFit", "ctFit"))
}

.ctmaHasJuliaFit <- function(object) {
  if (inherits(object, c("ctJuliaFit", "ctmaJuliaFit"))) return(TRUE)
  fits <- object$studyFitList
  if (.ctmaIsNativeFit(fits)) return(inherits(fits, c("ctJuliaFit", "ctmaJuliaFit")))
  is.list(fits) && any(vapply(fits, function(x)
    inherits(x, c("ctJuliaFit", "ctmaJuliaFit")), logical(1)))
}

#' Separate convergence, identification and scale diagnostics
#'
#' @description Reports stored fit diagnostics without re-fitting or starting
#' Julia. Optimizer convergence is distinguished from identification and from
#' posterior sampling and Laplace integration-accuracy diagnostics. Missing
#' diagnostics are reported as unknown.
#' @param object A CoTiMA2 Julia fit or a conventional CoTiMAFit object.
#' @param print Whether to print the diagnostic report.
#' @details New native ctsemjulia fits cache population matrices at fitting
#' time. A binary linear-predictor SD above 10 at Time 0 or at the population
#' stationary covariance is flagged as a scale heuristic. This is not an
#' identification test, a validated numerical cutoff, or an instruction to fix
#' diffusion. Older fits without cached matrices remain assessable for their
#' stored optimizer and identification diagnostics. A clean report does not
#' certify likelihood accuracy, standard errors or posterior convergence.
#' Stored raw-coordinate and transformed-summary diagnostics are labelled
#' separately. SMC's undefined R-hat is not a chain-convergence warning. Pooled
#' native storage and primary-study rosters are recognized independently.
#' Cached Laplace check status is reported independently of optimizer convergence.
#' A stored status of checked does not imply an accuracy pass; unresolved checks
#' remain unresolved. Reading this information never reruns the check.
#' Source-reported checked and Laplace log likelihoods, log-likelihood SE,
#' and optimizer predicted gain, resolution and convergence floor are displayed
#' separately when stored. Missing likelihood methods remain unknown.
#' @return Invisibly, a diagnostic list, or study-indexed lists for ctmaInit fits.
#' @export
ctmaFitDiagnostics <- function(object, print = TRUE) {
  if (!is.list(object)) stop("object must be a fitted model list.", call. = FALSE)
  if (.ctmaIsNativeFit(object$studyFitList)) {
    checks <- .ctmaFitChecks(object$studyFitList)
    if (isTRUE(print)) .ctmaPrintFitChecks(checks, "Pooled model")
    return(invisible(checks))
  }
  if (is.list(object$studyFitList) && length(object$studyFitList)) {
    studies <- lapply(object$studyFitList, .ctmaFitChecks)
    numbers <- object$primaryStudyList$studyNumbers
    if (length(numbers) != length(studies)) numbers <- object$statisticsList$originalStudyNumbers
    if (length(numbers) != length(studies)) numbers <- seq_along(studies)
    names(studies) <- paste("Study", unlist(numbers))
    if (isTRUE(print)) for (label in names(studies)) {
      .ctmaPrintFitChecks(studies[[label]], label)
    }
    return(invisible(studies))
  }
  checks <- .ctmaFitChecks(object)
  if (isTRUE(print)) .ctmaPrintFitChecks(checks)
  invisible(checks)
}

.ctmaStoredPosteriorDiagnostics <- function(object) {
  frames <- list()
  for (slot in c("resultsSummary", "ctsem_summary", "summary")) {
    summary <- object[[slot]]
    if (!is.list(summary)) next
    for (section in names(summary)) {
      tab <- summary[[section]]
      if (is.data.frame(tab) && any(c("Rhat", "n_eff") %in% names(tab))) {
        frames[[paste(slot, section, sep = "$")]] <- tab
      }
    }
  }
  values <- function(statistic) {
    out <- numeric()
    for (section in names(frames)) {
      tab <- frames[[section]]
      if (!is.numeric(tab[[statistic]])) next
      labels <- rownames(tab)
      if (is.null(labels)) labels <- as.character(seq_len(nrow(tab)))
      out <- c(out, stats::setNames(tab[[statistic]], paste(section, labels, sep = ": ")))
    }
    out
  }
  list(rhat = values("Rhat"), ess = values("n_eff"))
}

.ctmaFitChecks <- function(object) {
  first <- function(...) {
    for (x in list(...)) if (!is.null(x) && length(x)) return(x)
    NULL
  }
  flag <- function(x) if (is.logical(x) && length(x) == 1L && !is.na(x)) x else NA
  if (!is.list(object)) stop("A stored fit must be a list.", call. = FALSE)
  optim <- first(object$optim, object$optimization, list())
  identification <- object$identifiability
  sample <- first(object$sample, object$bayes$diagnostics, list())
  # The bundled NUTS result stores posterior::summarise_draws output directly,
  # rather than the native engine's sample$rhat/sample$ess schema.
  if (is.null(object$sample) && is.list(object$bayes) &&
      is.data.frame(object$bayes$summary)) {
    tab <- object$bayes$summary
    labels <- if (is.character(tab$variable)) tab$variable else rownames(tab)
    sample <- list(sampler = "nuts", divergent = object$bayes$divergences,
      rhat = if (is.numeric(tab$rhat)) stats::setNames(tab$rhat, labels),
      ess = if (is.numeric(tab$ess_bulk)) stats::setNames(tab$ess_bulk, labels))
    # There is no stored global convergence verdict: leave it unknown.
  }
  sampling <- !is.null(object$sample) || identical(tolower(object$estimation), "bayes")
  sampler <- first(sample$sampler, "unknown")
  rhat_applicable <- !identical(tolower(sampler), "smc")
  stored <- .ctmaStoredPosteriorDiagnostics(object)
  raw_rhat <- if (is.numeric(sample$rhat)) sample$rhat else numeric()
  all_rhat <- c(raw_rhat, stored$rhat)
  finite_max <- function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
  finite_min <- function(x) if (is.numeric(x) && any(is.finite(x))) min(x[is.finite(x)]) else NA_real_
  worst <- if (rhat_applicable && any(is.finite(all_rhat)))
    which.max(replace(all_rhat, !is.finite(all_rhat), -Inf)) else integer()
  weak <- first(identification$nweak, NA_integer_)
  saturated <- isTRUE(optim$saturated) || length(optim$saturated_parameters) > 0L
  raw <- first(object$estimate$raw, object$optimization$estimates)
  m <- object$ctsemjulia_diagnostics$population_matrices
  model <- first(object$modelbase, object$ctstanmodelbase, object$model)
  binary <- if (is.list(model)) which(model$manifesttype == 1) else integer()
  scale <- list(assessed = FALSE, flagged = character(), cutoff_sd = 10)
  if (length(binary) && is.list(m) && is.matrix(m$LAMBDA)) {
    for (label in c("T0cov", "asymDIFFUSIONcov")) {
      v <- m[[label]]; lambda <- m$LAMBDA
      if (!is.matrix(v) || ncol(lambda) != nrow(v) || nrow(v) != ncol(v)) next
      variances <- diag(lambda %*% v %*% t(lambda))
      valid <- binary[is.finite(variances[binary]) & variances[binary] >= 0]
      if (length(valid)) scale$assessed <- TRUE
      large <- valid[sqrt(variances[valid]) > scale$cutoff_sd]
      if (length(large)) scale$flagged <- c(scale$flagged,
        paste0("manifest ", large, " (", label, ", SD=",
          format(signif(sqrt(variances[large]), 4)), ")"))
    }
  }
  list(sampling = sampling,
    optimizer_converged = flag(first(optim$converged, optim$optimizer_converged)),
    gradient_norm = first(optim$gradient_norm, optim$maximum_absolute_gradient, NA_real_),
    likelihood = .ctmaStoredLikelihood(object),
    optimizer_resolution = list(predicted_gain = optim$predicted_gain,
      convergence_tolerance = optim$convergence_tolerance,
      resolution = optim$resolution, convergence_floor = optim$convergence_floor),
    laplace_accuracy = .ctmaStoredLaplaceAccuracy(object),
    unidentified_directions = weak,
    unidentified_parameters = first(identification$parameters, character()),
    negative_directions = first(identification$negative, NA_integer_),
    saturated = saturated, saturated_parameters = optim$saturated_parameters,
    nonfinite_estimates = if (is.numeric(raw)) any(!is.finite(raw)) else NA,
    scale = scale,
    posterior = list(sampler = sampler,
      reported_converged = flag(sample$converged),
      divergences = first(sample$divergent, NA_integer_),
      rhat_applicable = rhat_applicable,
      nonfinite_rhat = rhat_applicable && length(all_rhat) > 0L && any(!is.finite(all_rhat)),
      rhat_raw = if (rhat_applicable) finite_max(raw_rhat) else NA_real_,
      rhat_transformed = if (rhat_applicable) finite_max(stored$rhat) else NA_real_,
      rhat = if (rhat_applicable) finite_max(all_rhat) else NA_real_,
      worst_rhat_parameter = if (length(worst) && length(names(all_rhat))) names(all_rhat)[worst] else NA_character_,
      ess_raw = finite_min(sample$ess),
      ess_transformed = finite_min(stored$ess),
      ess_transformed_definition = "source-reported n_eff",
      rhat_target = first(sample$rhat_target, 1.01)))
}

.ctmaPrintFitChecks <- function(x, label = NULL) {
  cat("\nCoTiMA2 fit checks", if (!is.null(label)) paste0(" -- ", label), "\n", sep = "")
  verdict <- function(x) if (is.na(x)) "unknown" else if (x) "yes" else "NO"
  cat(if (x$sampling) "  Placement optimizer converged: " else "  Optimizer converged: ",
    verdict(x$optimizer_converged), "; gradient diagnostic: ",
    format(signif(x$gradient_norm, 4)), "\n", sep = "")
  likelihood <- x$likelihood
  if (is.finite(likelihood$loglik)) {
    cat("  Reported log likelihood (", likelihood$method, "): ",
      format(likelihood$loglik, digits = 10), "\n", sep = "")
    if (is.finite(likelihood$loglik_se)) cat("    Reported log-likelihood SE: ",
      format(signif(likelihood$loglik_se, 4)), "\n", sep = "")
    if (is.finite(likelihood$loglik_laplace)) cat("    Laplace log likelihood: ",
      format(likelihood$loglik_laplace, digits = 10),
      "; keep separate from the reported likelihood.\n", sep = "")
  }
  for (field in names(x$optimizer_resolution)) {
    value <- x$optimizer_resolution[[field]]
    if (is.numeric(value) && length(value) == 1L && is.finite(value)) cat(
      "    Optimizer ", gsub("_", " ", field), ": ",
      format(signif(value, 5)), "\n", sep = "")
  }
  if (isTRUE(x$laplace_accuracy$applicable)) cat("  ",
    .ctmaLaplaceAccuracyNote(x$laplace_accuracy), "\n", sep = "")
  cat("  Identification: ", if (is.na(x$unidentified_directions)) {
    "not assessed in stored diagnostics"
  } else paste(x$unidentified_directions, "unidentified direction(s) reported"), "\n", sep = "")
  if (length(x$unidentified_parameters)) cat("    Parameters: ",
    paste(x$unidentified_parameters, collapse = ", "), "\n", sep = "")
  if (is.finite(x$negative_directions) && x$negative_directions > 0) cat(
    "    Negative-curvature directions reported: ", x$negative_directions, "\n", sep = "")
  if (x$saturated) cat("  Parameter-transform saturation reported",
    if (length(x$saturated_parameters)) paste0(": ", paste(x$saturated_parameters, collapse = ", ")),
    ". Small gradients may be misleading.\n", sep = "")
  if (isTRUE(x$nonfinite_estimates)) cat("  Non-finite estimates present.\n")
  if (length(x$scale$flagged)) cat("  Binary scale flag (heuristic, predictor SD > 10): ",
    paste(x$scale$flagged, collapse = "; "),
    ". Check likelihood accuracy and sensitivity; do not automatically fix diffusion.\n", sep = "")
  if (!x$scale$assessed) cat("  Binary scale: not assessed (no applicable cached covariance).\n")
  if (x$scale$assessed && !length(x$scale$flagged)) cat(
    "  Binary scale: no predictor SD > 10 flag in cached matrices (heuristic).\n")
  if (x$sampling) cat("  Posterior diagnostics (", x$posterior$sampler,
    "): reported convergence=", verdict(x$posterior$reported_converged),
    "; divergences=", x$posterior$divergences,
    "; maximum finite R-hat=", if (isFALSE(x$posterior$rhat_applicable))
      "not applicable for SMC" else format(signif(x$posterior$rhat, 4)),
    ". Missing fields require sampler-specific checks.\n", sep = "")
  if (x$sampling && isTRUE(x$posterior$rhat_applicable)) {
    cat("    Raw-coordinate maximum R-hat: ", format(signif(x$posterior$rhat_raw, 4)),
      "; stored transformed-summary maximum: ", format(signif(x$posterior$rhat_transformed, 4)), "\n", sep = "")
    cat("    Minimum raw ESS: ", format(signif(x$posterior$ess_raw, 4)),
      "; minimum stored transformed n_eff: ", format(signif(x$posterior$ess_transformed, 4)), "\n", sep = "")
    if (is.finite(x$posterior$rhat) && x$posterior$rhat > x$posterior$rhat_target) cat(
      "    R-hat exceeds target ", x$posterior$rhat_target, ": ",
      x$posterior$worst_rhat_parameter, ". Inspect the affected posterior quantities.\n", sep = "")
  }
  if (x$sampling && isFALSE(x$posterior$rhat_applicable)) cat("    Minimum sample-reported ESS: ",
    format(signif(x$posterior$ess_raw, 4)), "; R-hat is intentionally undefined for SMC.\n", sep = "")
  if (x$sampling && x$posterior$nonfinite_rhat) cat("    Non-finite R-hat values are present.\n")
  cat("  Optimizer convergence does not establish identification, reliable SEs or likelihood accuracy.\n")
  invisible(x)
}

#' @export
summary.ctmaCtsemJuliaFit <- function(object, ..., pointEstimate = NULL) {
  answer <- NextMethod("summary")
  if(is.null(pointEstimate)) pointEstimate <- object$argumentList$pointEstimate
  if(is.null(pointEstimate)) pointEstimate <- "auto"
  answer <- .ctmaPointSummary(answer,object,pointEstimate)
  attr(answer, "cotima_fit_checks") <- .ctmaFitChecks(object)
  class(answer) <- unique(c("summary.ctmaCtsemJuliaFit", class(answer)))
  answer
}

#' @export
print.summary.ctmaCtsemJuliaFit <- function(x, ...) {
  NextMethod("print")
  .ctmaPrintPointNote(x)
  .ctmaPrintFitChecks(attr(x, "cotima_fit_checks"))
  invisible(x)
}
