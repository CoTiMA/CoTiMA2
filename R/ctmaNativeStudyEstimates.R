.ctmaNativeStudyFit <- function(fit) {
  if (inherits(fit, "ctJuliaFit")) return(fit)
  if (inherits(fit$studyFitList, "ctJuliaFit")) return(fit$studyFitList)
  NULL
}

.ctmaNativeStudyProfiles <- function(fit, study_design = NULL) {
  if (!is.null(fit$model_spec$laplace$nrandom) && fit$model_spec$laplace$nrandom > 0)
    stop("Study reconstruction requires a fit with study contrasts and no hierarchical/person-level random effects.", call. = FALSE)
  pars <- fit$modelbase$pars
  if (is.null(pars)) pars <- fit$model$pars
  if (is.data.frame(pars)) {
    flags <- grep("^indvarying($|_)", names(pars), value = TRUE)
    if (any(vapply(flags, function(v) any(pars[[v]] %in% TRUE), logical(1))))
      stop("Study reconstruction requires no hierarchical/person-level random effects, including augmented-state effects.", call. = FALSE)
  }
  ti <- fit$model_spec$tipred_data
  if (!is.matrix(ti) || !ncol(ti) || any(!is.finite(ti)))
    stop("The native fit must contain finite study contrast predictors.", call. = FALSE)
  first <- which(!duplicated(as.data.frame(ti)))
  design <- ti[first, , drop = FALSE]
  if (nrow(design) != ncol(design) + 1L || qr(cbind(1, design))$rank != nrow(design))
    stop("Expected a saturated full-rank study contrast design, with K profiles and K-1 predictors. Additional moderators are not supported by this native reconstruction path.", call. = FALSE)
  if (!is.null(study_design) && !isTRUE(all.equal(unname(as.matrix(study_design)),
      unname(design), check.attributes = FALSE, tolerance = 1e-12)))
    stop("For native fits, study_design must match the fitted profile order; use the returned studyDesign to identify profiles.", call. = FALSE)
  rownames(design) <- paste0("profile_", seq_len(nrow(design)))
  list(first = first, design = design)
}

.ctmaNativeStudyExtract <- function(fit, profiles, cores, raw = NULL) {
  if (!is.null(raw)) {
    fit$estimate$rawposterior <- as.matrix(raw)
    colnames(fit$estimate$rawposterior) <- colnames(fit$estimate$cov)
  }
  out <- ctsem::ctExtract(fit, subjectMatrices = TRUE,
    subjects = profiles$first, realid = FALSE, nsamples = "all", cores = cores)
  if (!is.null(raw) && !isTRUE(all.equal(unname(out$rawpars), unname(as.matrix(raw)),
      tolerance = 0, check.attributes = FALSE)))
    stop("The native extractor did not retain all requested parameter evaluations.", call. = FALSE)
  a <- out$subj_DRIFT
  if (length(dim(a)) != 4L || dim(a)[2L] != nrow(profiles$design) || any(!is.finite(a)))
    stop("Native extraction did not return finite study DRIFT matrices.", call. = FALSE)
  out
}

#' Reconstruct native Julia study DRIFT estimates and joint covariance
#'
#' @param fit A native ctsem Julia fit or pooled CoTiMA2 fit with a saturated
#'   study-indicator design and no hierarchical or person-level random effects.
#' @param cores Positive integer number of extraction cores.
#' @param derivativeStep Positive finite relative finite-difference step.
#' @return A list containing studyEstimates (ML modes and delta-method SEs),
#'   covariance (including cross-study covariance), jacobian, studyDesign,
#'   representativeSubjects and scaleTime. DRIFT is on the original time scale.
#' @details Uses the optimizer estimate and its full raw covariance, not means
#'   or SDs of approximate draws. The native extractor evaluates the model at
#'   the estimate and central perturbations, preserving nonlinear transforms.
#'   Profiles are ordered by first occurrence in the fitted data; studyDesign
#'   identifies them. Bayesian fits and additional moderators are not supported.
#' @export
ctmaStudyEstimates <- function(fit, cores = 1L, derivativeStep = 1e-5) {
  native <- .ctmaNativeStudyFit(fit)
  if (is.null(native)) stop("ctmaStudyEstimates currently requires a native Julia fit.", call. = FALSE)
  optimize <- native$args$resolved$optimize
  if (is.null(optimize)) optimize <- native$args$input$optimize
  if (is.null(optimize)) optimize <- native$args$optimize
  if (!isTRUE(optimize)) stop("ctmaStudyEstimates requires an ML optimization fit.", call. = FALSE)
  if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) || cores < 1 || cores != as.integer(cores))
    stop("cores must be one positive integer.", call. = FALSE)
  if (!is.numeric(derivativeStep) || length(derivativeStep) != 1L ||
      !is.finite(derivativeStep) || derivativeStep <= 0)
    stop("derivativeStep must be finite and positive.", call. = FALSE)
  profiles <- .ctmaNativeStudyProfiles(native)
  z <- as.numeric(native$estimate$raw)
  v <- native$estimate$cov
  p <- length(z)
  if (!p || any(!is.finite(z)) || !is.matrix(v) || !identical(dim(v), c(p, p)) ||
      any(!is.finite(v)) || !isSymmetric(unname(v)))
    stop("A finite estimate and matching symmetric raw covariance are required.", call. = FALSE)
  h <- derivativeStep * pmax(1, abs(z))
  points <- rbind(z, sweep(diag(h, p), 2L, z, "+"), sweep(-diag(h, p), 2L, z, "+"))
  extracted <- .ctmaNativeStudyExtract(native, profiles, cores, points)
  a <- extracted$subj_DRIFT
  values <- matrix(a, nrow = nrow(points))
  jac <- t(sweep(values[1L + seq_len(p), , drop = FALSE] -
    values[1L + p + seq_len(p), , drop = FALSE], 1L, 2 * h, "/"))
  scale <- .ctmaSummaryTimeFactor(fit)
  values <- values * scale
  jac <- jac * scale
  covariance <- jac %*% v %*% t(jac)
  # Restore exact symmetry after floating-point matrix multiplication.
  covariance <- (covariance + t(covariance)) / 2
  index <- expand.grid(study = rownames(profiles$design), row = seq_len(dim(a)[3L]),
    col = seq_len(dim(a)[4L]), stringsAsFactors = FALSE)
  index$parameter <- paste0("DRIFT[", index$row, ",", index$col, "]")
  index$estimate <- as.numeric(values[1L, ])
  index$standard_error <- sqrt(pmax(0, diag(covariance)))
  labels <- paste(index$study, index$parameter, sep = ":")
  dimnames(covariance) <- list(labels, labels)
  rownames(jac) <- labels
  list(studyEstimates = index, covariance = covariance, jacobian = jac,
    studyDesign = profiles$design, representativeSubjects = profiles$first,
    scaleTime = scale, estimator = "optimizer mode; delta-method covariance from full raw covariance")
}

.ctmaReMomentsNative <- function(fit, native, study_design, probs, weighting,
    SEsource, cores, include_matrices, include_individual) {
  profiles <- .ctmaNativeStudyProfiles(native, study_design)
  raw <- native$estimate$rawposterior
  if (!is.matrix(raw) || nrow(raw) < 2L)
    stop("Native moment summaries require at least two joint uncertainty draws.", call. = FALSE)
  source <- .ctmaReMomentsSESource(native, SEsource)
  ex <- .ctmaNativeStudyExtract(native, profiles, cores)
  a <- ex$subj_DRIFT
  n <- dim(a)[3L]
  labels <- as.vector(outer(seq_len(n), seq_len(n), function(r, c) paste0("DRIFT[", r, ",", c, "]")))
  transformed <- array(aperm(a, c(1L, 3L, 4L, 2L)),
    dim = c(dim(a)[1L], n * n, dim(a)[2L]),
    dimnames = list(NULL, labels, rownames(profiles$design)))
  moments <- .ctmaReMomentsDistribution(transformed, probs, weighting)
  matrices <- NULL
  if (isTRUE(include_matrices)) {
    matrices <- list(DRIFT = aperm(a, c(1L, 3L, 4L, 2L)))
    for (m in c("DIFFUSIONcov", "T0cov")) if (!is.null(ex[[paste0("subj_", m)]]))
      matrices[[m]] <- aperm(ex[[paste0("subj_", m)]], c(1L, 3L, 4L, 2L))
  }
  if (length(matrices)) for (m in names(matrices)) dimnames(matrices[[m]]) <-
    list(NULL, seq_len(dim(matrices[[m]])[2L]), seq_len(dim(matrices[[m]])[3L]), rownames(profiles$design))
  map <- expand.grid(row = seq_len(n), col = seq_len(n))
  map$parameter <- labels
  map$matrix <- "DRIFT"
  out <- structure(list(
    method = "post-hoc empirical moments of native study DRIFT draws",
    weighting = weighting, covarianceWeighting = "equal study weights", SEsource = source,
    pooledFit = native, scaleTime = .ctmaSummaryTimeFactor(fit),
    parameterMap = list(rawpopmeans = map, note = "Evaluated DRIFT matrix entries, not a raw-coordinate mapping."),
    studyDesign = profiles$design, representativeSubjects = profiles$first,
    rawStudyDraws = NULL, transformedStudyDraws = transformed,
    drawsRequested = nrow(raw), drawsRetained = dim(a)[1L], drawsDropped = nrow(raw) - dim(a)[1L],
    studyParameters = moments$study, distributionParameters = moments$distribution,
    studySE = moments$study_se, studyWeights = moments$weights,
    effectiveNumberOfStudies = moments$effective_studies,
    covarianceDraws = moments$covariance, correlationDraws = moments$correlation,
    studyMatrices = matrices, matrixSummary = if (is.null(matrices)) NULL else .ctmaReMomentsMatrixSummary(matrices, probs),
    individualRandomEffects = NULL, probs = probs), class = c("CoTiMAFit", "CoTiMAReMoments"))
  out$summary <- .ctmaReMomentsSummary(out, probs, print_note = FALSE)
  .ctma_print_message("Native Julia study DRIFT moments use transformed joint draws. Their means and SDs are not ML modes and Hessian SEs; use ctmaStudyEstimates() for those. The WEC raw intercept is not generally an arithmetic mean on the DRIFT scale.")
  out
}
