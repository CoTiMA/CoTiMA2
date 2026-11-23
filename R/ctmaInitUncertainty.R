# Keep native no-width diagnostics separate from numerical covariance projections.
.ctmaInitNoWidth <- function(fit) {
  # This is ctsem's explicit reporting decision, not a cutoff on eigenvalues
  # or a guess based on every parameter mentioned in a weak direction.
  if (!is.null(fit$sample$chains)) return(character())
  unique(as.character(fit$uncertainty$intervalcheck$unidentified))
}

.ctmaInitMaskNativeSummary <- function(fit, result = fit$resultsSummary) {
  bad <- .ctmaInitNoWidth(fit)
  if (!length(bad) || !is.list(result)) return(result)
  p <- fit$model$pars
  if (is.null(p)) p <- fit$ctstanmodelbase$pars
  tab <- result$parmatrices
  if (!is.data.frame(p) || !is.data.frame(tab) ||
      !all(c("matrix", "row", "col", "param") %in% names(p)) ||
      !all(c("matrix", "row", "col") %in% names(tab))) return(result)
  # Own parameter cells only. Derived covariance families are withheld
  # conservatively if their source parameters have no estimable width.
  flagged <- rep(FALSE, nrow(tab))
  for (j in which(!is.na(p$param) & p$param %in% bad)) {
    flagged <- flagged | (tab$matrix == p$matrix[j] &
      tab$row == p$row[j] & tab$col == p$col[j])
  }
  source_bad <- unique(p$matrix[!is.na(p$param) & p$param %in% bad])
  derived <- c(DIFFUSION = "DIFFUSIONcov", T0VAR = "T0cov",
    MANIFESTVAR = "MANIFESTcov")
  flagged <- flagged | tab$matrix %in% unname(derived[intersect(names(derived), source_bad)])
  # Stationary T0 covariance depends on both A and Q; DT effects depend on A.
  stationary <- isTRUE(as.logical(fit$model$stationary)) || isTRUE(fit$args$stationary)
  if ("DRIFT" %in% source_bad) flagged <- flagged | tab$matrix %in% c("dtDRIFT", "DRIFTdiscrete")
  if (stationary && any(c("DRIFT", "DIFFUSION") %in% source_bad))
    flagged <- flagged | tab$matrix %in% c("T0cov", "ASYMDIFFUSIONcov")
  uncertainty <- tolower(names(tab)) %in% c("sd", "se", "2.5%", "50%", "97.5%", "z")
  for (column in names(tab)[uncertainty]) tab[[column]][flagged] <- NA_real_
  result$parmatrices <- tab
  result$uncertainty_note <- paste0("Uncertainty is unavailable for ctsem-flagged coordinates: ",
    paste(bad, collapse = ", "), ". Derived covariance uncertainty is withheld when its source has no estimable width; point estimates are retained.")
  result
}

.ctmaInitWaveDiagnostic <- function(data, model, random_intercepts = FALSE) {
  if (!is.data.frame(data) || !all(c("id", "time") %in% names(data))) return(NULL)
  manifests <- model$manifestNames
  keep <- is.finite(data$time) & !is.na(data$id)
  available <- intersect(manifests, names(data))
  if (length(available)) keep <- keep & rowSums(!is.na(data[, available, drop = FALSE])) > 0L
  d <- data[keep, , drop = FALSE]
  counts <- vapply(split(d$time, d$id), function(t) length(unique(t)), integer(1))
  p <- model$pars
  has_random <- is.character(random_intercepts) && any(toupper(random_intercepts) %in% c("MANIFEST", "CINT"))
  if (is.data.frame(p)) has_random <- has_random || any(p$matrix %in% c("MANIFESTMEANS", "CINT") & p$indvarying %in% TRUE)
  if (!length(counts) || !has_random || max(counts) > 2L) return(NULL)
  paste0("At most two observed waves per person in a separately fitted primary study with its own random-intercept covariance: an unrestricted within-person dynamic process and random-intercept covariance cannot generally be separated without additional identifying constraints. More persons do not resolve this design limitation. Unequal observation intervals or other model constraints can change identification; inspect ctsem's identification diagnostics. The model is retained unchanged.")
}

.ctmaRefreshInitUncertainty <- function(object) {
  fits <- object$studyFitList
  if (!is.list(fits) || length(object$n.studies) != 1L || length(fits) != object$n.studies || !is.list(object$summary)) return(object)
  notes <- character()
  mask_family <- function(table, row, family, cells, ci = FALSE) {
    if (is.null(table) || !is.matrix(table) || nrow(table) < row) return(table)
    if (ci) {
      cols <- as.vector(rbind(2L * cells - 1L, 2L * cells))
    } else cols <- 1L + 2L * cells
    cols <- cols[cols <= ncol(table)]
    table[row, cols] <- NA
    table
  }
  for (i in seq_along(fits)) {
    f <- fits[[i]]
    s <- .ctmaInitMaskNativeSummary(f)
    object$studyFitList[[i]]$resultsSummary <- s
    label <- sub("^Study No ", "", rownames(object$summary$estimates)[i])
    bad <- .ctmaInitNoWidth(f)
    if (length(bad)) notes <- c(notes, paste0("Study ", label,
      ": unavailable uncertainty for ", paste(bad, collapse = ", ")))
    p <- s$parmatrices
    if (!is.data.frame(p)) next
    for (family in c("DRIFT", "DIFFUSIONcov", "T0cov")) {
      q <- p[p$matrix == family, , drop = FALSE]
      if (!nrow(q) || !"sd" %in% names(q)) next
      q <- q[order(q$col, q$row), , drop = FALSE]
      cells <- which(is.na(q$sd))
      if (!length(cells)) next
      if (family == "DRIFT") {
        for (field in c("estimates", "drift_estimates_original_time_scale"))
          object$summary[[field]] <- mask_family(object$summary[[field]], i, family, cells)
        for (field in c("drift_CI_original_time_scale"))
          object$summary[[field]] <- mask_family(object$summary[[field]], i, family, cells, TRUE)
        # General CI table: DRIFT, DIFFUSION, T0VAR, each with lower/upper pairs.
        offset <- 0L
      } else if (family == "DIFFUSIONcov") {
        object$summary$diff_estimates_original_time_scale <- mask_family(
          object$summary$diff_estimates_original_time_scale, i, family, cells)
        object$summary$diff_CI_original_time_scale <- mask_family(
          object$summary$diff_CI_original_time_scale, i, family, cells, TRUE)
        offset <- 2L * object$n.latent^2
      } else offset <- 4L * object$n.latent^2
      cols <- offset + as.vector(rbind(2L * cells - 1L, 2L * cells))
      if (is.matrix(object$summary$confidenceIntervals))
        object$summary$confidenceIntervals[i, cols] <- NA_real_
    }
    if (length(f$uncertainty$intervalcheck$nflagged) && f$uncertainty$intervalcheck$nflagged > 0L)
      notes <- c(notes, paste0("Study ", label, ": unstable curvature-based interval widths for ",
        paste(f$uncertainty$intervalcheck$parameters, collapse = ", "), "."))
    if (length(f$identifiability$nweak) && f$identifiability$nweak > 0L)
      notes <- c(notes, paste0("Study ", label, ": ", f$identifiability$nweak,
        " weak parameter direction(s); optimizer convergence does not establish identification."))
    if (identical(f$uncertainty$certification$status, "saturated"))
      notes <- c(notes, paste0("Study ", label, ": saturated parameter transformation."))
    if (is.data.frame(f$corbound) && any(f$corbound$at_bound %in% TRUE))
      notes <- c(notes, paste0("Study ", label, ": random-effect correlation at its boundary."))
    wave <- .ctmaInitWaveDiagnosticFromFit(f)
    if (!is.null(wave)) notes <- c(notes, paste0("Study ", label, ": ", wave))
  }
  object$summary$identification_notes <- if (length(notes))
    structure(unique(notes), class = c("ctmaIdentificationNotes", "character")) else character()
  object
}

.ctmaInitWaveDiagnosticFromFit <- function(fit) {
  d <- fit$data
  if (!is.list(d) || is.null(d$subject) || is.null(d$time)) return(NULL)
  dat <- data.frame(id = d$subject, time = as.numeric(d$time))
  if (is.matrix(d$Y) && nrow(d$Y) == nrow(dat)) {
    keep <- rowSums(!is.na(d$Y)) > 0L
    dat <- dat[keep, , drop = FALSE]
  }
  .ctmaInitWaveDiagnostic(dat, fit$model)
}

# Undefined transformed draws cannot supply an interval for that coordinate.
# Preserve finite-coordinate behavior and keep unavailable quantiles as NA.
.ctmaInitDrawQuantile <- function(x, probability) {
  if (!length(x) || any(!is.finite(x))) return(NA_real_)
  stats::quantile(x, probability)
}
