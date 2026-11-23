# Observed study design and ctmaPrep-compatible, non-destructive partitions.
.ctmaObservedStudyDesign <- function(data, manifest_names = NULL) {
  empty <- list(nWaves = NA_integer_, minWaves = NA_integer_, maxWaves = NA_integer_,
    nSubjects = 0L, nSubjectsTwoWaves = 0L, nSubjectsThreeOrMoreWaves = 0L)
  if (!is.data.frame(data) || !all(c("id", "time") %in% names(data))) return(empty)
  keep <- !is.na(data$id) & is.finite(data$time)
  manifest_names <- intersect(manifest_names, names(data))
  if (length(manifest_names))
    keep <- keep & rowSums(!is.na(data[, manifest_names, drop = FALSE])) > 0L
  d <- data[keep, , drop = FALSE]
  counts <- vapply(split(d$time, as.character(d$id), drop = TRUE),
    function(t) length(unique(t)), integer(1))
  if (!length(counts)) return(empty)
  list(nWaves = max(counts), minWaves = min(counts), maxWaves = max(counts),
    nSubjects = length(counts), nSubjectsTwoWaves = sum(counts == 2L),
    nSubjectsThreeOrMoreWaves = sum(counts >= 3L))
}

.ctmaPrepStudySubset <- function(studies, indices) {
  if (!length(indices)) return(NULL)
  n <- studies$n.studies
  if (length(n) != 1L || !is.finite(n) || n < 1L ||
      anyNA(indices) || any(indices < 1L | indices > n) || anyDuplicated(indices))
    stop("Cannot partition an inconsistent primary-study list.", call. = FALSE)
  out <- studies
  # Global metadata can coincidentally have the same length as the study list.
  globals <- c("n.studies", "moderatorLabels", "moderatorValues", "activeDirectory",
    "plot.type", "excelSheets", "reUseEmprawData", "summary")
  for (name in setdiff(names(studies), globals)) {
    value <- studies[[name]]
    if (is.list(value) && !is.data.frame(value) && length(value) == n)
      out[[name]] <- value[indices]
  }
  if (is.data.frame(studies$summary) && nrow(studies$summary) == n)
    out$summary <- studies$summary[indices, , drop = FALSE]
  out$n.studies <- length(indices)
  if ("excelSheets" %in% names(out)) {
    # Never retain a workbook containing the unselected studies or mutate its
    # reference object in place. Rebuild the available subset summary only.
    out["excelSheets"] <- list(NULL)
    if (is.data.frame(out$summary)) {
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "summary")
      openxlsx::writeData(wb, "summary", out$summary)
      out$excelSheets <- wb
    }
  }
  out
}

.ctmaStudyDesignTable <- function(studies, long_data, model = NULL,
    random_intercepts = FALSE) {
  n <- studies$n.studies
  if (length(n) != 1L || length(long_data) != n) return(NULL)
  p <- model$pars
  varying <- is.character(random_intercepts) &&
    any(toupper(random_intercepts) %in% c("CINT", "MANIFEST"))
  if (is.data.frame(p)) varying <- varying ||
    any(p$matrix %in% c("CINT", "MANIFESTMEANS") & p$indvarying %in% TRUE)
  rows <- lapply(seq_len(n), function(i) {
    d <- .ctmaObservedStudyDesign(long_data[[i]], model$manifestNames)
    number <- studies$studyNumbers[[i]]
    if (!length(number)) number <- i
    group <- if (is.na(d$nWaves)) "unknown" else if (d$nWaves == 2L)
      "twoWaves" else if (d$nWaves >= 3L) "threeOrMoreWaves" else "fewerThanTwoWaves"
    data.frame(study = as.character(number), as.data.frame(d), waveGroup = group,
      randomInterceptIdentificationRisk = if (is.na(d$nWaves)) NA else
        varying && d$nWaves <= 2L, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

.ctmaAttachStudyDesign <- function(object, studies, long_data, model = NULL,
    random_intercepts = FALSE) {
  design <- .ctmaStudyDesignTable(studies, long_data, model, random_intercepts)
  if (is.null(design)) return(object)
  object$studyDesign <- design
  object$studyLists <- list(
    twoWaves = .ctmaPrepStudySubset(studies, which(design$waveGroup == "twoWaves")),
    threeOrMoreWaves = .ctmaPrepStudySubset(studies, which(design$waveGroup == "threeOrMoreWaves")))
  object
}

.ctmaInitObservedLong <- function(object) {
  n <- object$n.studies
  fits <- object$studyFitList
  if (length(n) != 1L || !is.list(fits) || length(fits) != n) return(NULL)
  lapply(seq_len(n), function(i) {
    f <- fits[[i]]
    d <- f$data
    if (is.list(d) && !is.null(d$subject) && !is.null(d$time)) {
      out <- data.frame(id = d$subject, time = as.numeric(d$time))
      if (is.matrix(d$Y) && nrow(d$Y) == nrow(out))
        out <- out[rowSums(!is.na(d$Y)) > 0L, , drop = FALSE]
      return(out)
    }
    # Legacy primary fits retain wide data. Select outcome columns by their
    # numeric occasion suffix, not lexicographic order or a count of rows.
    wide <- f$empraw
    if (is.null(wide)) wide <- object$emprawList[[i]]
    if (!is.data.frame(wide) && !is.matrix(wide)) return(NULL)
    columns <- grep("_T[0-9]+$", colnames(wide), value = TRUE)
    if (!length(columns)) return(NULL)
    waves <- unique(sub("^.*_T", "", columns))
    pieces <- lapply(waves, function(w) {
      c <- columns[sub("^.*_T", "", columns) == w]
      seen <- rowSums(!is.na(wide[, c, drop = FALSE])) > 0L
      data.frame(id = which(seen), time = as.numeric(w))
    })
    do.call(rbind, pieces)
  })
}

.ctmaWarnStudySpecificTwoWaves <- function(init, structure) {
  # Native individual random effects shared across studies are NOT subjected
  # to a separate-primary-study identification warning.
  if (!identical(structure, "study-specific")) return(invisible(NULL))
  design <- init$studyDesign
  if (!is.data.frame(design)) {
    long <- .ctmaInitObservedLong(init)
    if (is.null(long)) return(invisible(NULL))
    studies <- init$primaryStudyList
    if (is.null(studies)) return(invisible(NULL))
    design <- .ctmaStudyDesignTable(studies, long)
  }
  if (is.null(design)) return(invisible(NULL))
  ids <- design$study[!is.na(design$maxWaves) & design$maxWaves <= 2L]
  if (length(ids)) warning(paste0(
    "ctmaFit: Study-specific random-intercept covariances requested for studies ",
    paste(ids, collapse = ", "), " with at most two observed waves per person. ",
    "These covariances and unrestricted within-person dynamics may require additional identifying constraints. ",
    "A common covariance requested through indVarying is a different pooled model; identification must be assessed jointly. ",
    "All studies are retained; no model constraints were changed."), call. = FALSE)
  invisible(NULL)
}
