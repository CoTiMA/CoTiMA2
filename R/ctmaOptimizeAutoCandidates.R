# Descriptive nominations from stored evidence only: no refitting, Hessian
# computation, likelihood evaluation, or inference that a study is defective.
.ctmaOptimizeIsAuto <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) &&
    identical(tolower(trimws(x)), "auto")
}

.ctmaOptimizeAutoCandidates <- function(context, candidates, references) {
  auto_candidates <- .ctmaOptimizeIsAuto(candidates)
  auto_references <- .ctmaOptimizeIsAuto(references)
  if (!auto_candidates && !auto_references) return(NULL)
  ids <- context$ids
  n <- length(ids)
  visits <- vapply(seq_len(n), function(i) {
    rows <- context$data[context$group == i, , drop = FALSE]
    counts <- vapply(split(rows$time, rows$id), function(t) {
      length(unique(t[is.finite(t)]))
    }, integer(1))
    stats::median(counts)
  }, numeric(1))
  tab <- data.frame(study = ids, n_subjects = context$counts,
    median_distinct_times = visits, localized_weak_directions = integer(n),
    primary_nonconverged = rep(FALSE, n), primary_weak_directions = integer(n),
    stringsAsFactors = FALSE)
  reasons <- replicate(n, character(), simplify = FALSE)
  notes <- character()
  evidence <- FALSE

  # Resolve reported diagnostic labels through the actual free TI layout.
  # Do not infer scientific identities from raw covariance-matrix row names.
  ti <- context$layout$ti
  ti <- ti[ti$predictor %in% context$columns, , drop = FALSE]
  label <- paste0("rawtipredeffect_", ti$base, "_", ti$predictor)
  key <- paste0("ti:", ti$base, "|", ti$predictor)
  direction_list <- context$native$identifiability$directions
  if (length(direction_list)) evidence <- TRUE
  for (d in direction_list) {
    p <- d$parameters
    w <- d$loadings
    if (!is.character(p) || !length(p) || length(p) > 4L ||
        !is.numeric(w) || length(w) != length(p) ||
        any(!is.finite(w)) || sum(w^2) == 0) next
    at <- match(p, label)
    missing <- is.na(at)
    at[missing] <- match(p[missing], key)
    # A shared/unmapped coefficient prevents study-local attribution. These
    # are reported loadings only; their localization is not causal evidence.
    if (anyNA(at) || anyDuplicated(label)) next
    study <- match(unique(ti$predictor[at]), context$columns)
    for (i in study) {
      tab$localized_weak_directions[i] <- tab$localized_weak_directions[i] + 1L
      reasons[[i]] <- c(reasons[[i]], "reported weak direction localized to study contrasts")
    }
  }

  # Optional retained primary fits can provide additional nominations, but
  # only when original IDs verify the alignment (never by position alone).
  initial <- context$init
  if (is.null(initial)) initial <- context$fit$ctmaInitFit
  if (is.list(initial) && is.list(initial$studyList) &&
      is.list(initial$studyFitList) &&
      !inherits(initial$studyFitList, c("ctJuliaFit", "ctStanFit", "ctmaJuliaFit"))) {
    primary_ids <- vapply(initial$studyList, function(s) {
      if (length(s$originalStudyNo) == 1L) as.character(s$originalStudyNo) else NA_character_
    }, character(1))
    if (!anyNA(primary_ids) && !anyDuplicated(primary_ids) &&
        length(primary_ids) == length(initial$studyFitList)) {
      for (i in seq_len(n)) {
        j <- match(ids[i], primary_ids)
        if (is.na(j)) next
        f <- initial$studyFitList[[j]]
        if (!is.list(f)) next
        if (inherits(f$studyFitList, c("ctJuliaFit", "ctStanFit", "ctmaJuliaFit"))) f <- f$studyFitList
        if (!is.list(f)) next
        op <- if (is.list(f$optim)) f$optim else list()
        opt <- if (is.list(f$optimization)) f$optimization else list()
        stan <- if (is.list(f$stanfit) && is.list(f$stanfit$optimfit)) f$stanfit$optimfit else list()
        converged <- if (!is.null(op$converged)) op$converged else
          if (!is.null(opt$converged)) opt$converged else stan$converged
        weak <- f$identifiability$nweak
        if (is.logical(converged) && length(converged) == 1L && !is.na(converged)) {
          evidence <- TRUE
          if (!converged) {
            tab$primary_nonconverged[i] <- TRUE
            reasons[[i]] <- c(reasons[[i]], "retained primary fit explicitly not converged")
          }
        }
        if (is.numeric(weak) && length(weak) == 1L && is.finite(weak) && weak > 0) {
          evidence <- TRUE
          tab$primary_weak_directions[i] <- as.integer(weak)
          reasons[[i]] <- c(reasons[[i]], "retained primary fit reports weak directions")
        }
      }
    } else notes <- c(notes, "Primary-fit diagnostics ignored: original study IDs could not verify their alignment.")
  }

  tab$flagged <- tab$localized_weak_directions > 0L |
    tab$primary_nonconverged | tab$primary_weak_directions > 0L
  tab$reason <- vapply(reasons, function(x) paste(unique(x), collapse = "; "), character(1))
  tab$exclusion_rank <- tab$reference_rank <- rep(NA_integer_, n)
  tab$nominated_exclusion <- tab$nominated_reference <- rep(FALSE, n)
  selected_candidates <- selected_references <- character()
  if (auto_candidates) {
    ranked <- which(tab$flagged)
    ranked <- ranked[order(!tab$primary_nonconverged[ranked],
      -tab$localized_weak_directions[ranked], -tab$primary_weak_directions[ranked],
      tab$n_subjects[ranked], ranked)]
    if (length(ranked)) tab$exclusion_rank[ranked] <- seq_along(ranked)
    chosen <- head(ranked, min(10L, max(0L, n - 2L)))
    tab$nominated_exclusion[chosen] <- TRUE
    selected_candidates <- ids[chosen]
    if (!length(chosen)) notes <- c(notes, if (n < 3L)
      "No automatic exclusions: at least two studies must remain." else if (!evidence)
      "No study-local diagnostic evidence was retained; no exclusion candidates nominated. Small N or few occasions alone do not identify a problem study." else
      "No study-local diagnostic flags met the nomination rule; no exclusion candidates nominated.")
    if (length(ranked) > length(chosen)) notes <- c(notes,
      "Automatic exclusion nominations were capped at ten and at nstudies minus two; all flags remain in studyDiagnostics.")
  }
  if (auto_references) {
    flags <- do.call(cbind, lapply(paste0(context$columns, "_effect"),
      function(nm) as.logical(context$model$pars[[nm]])))
    flags[is.na(flags)] <- FALSE
    partial <- any(rowSums(flags) > 0 & rowSums(flags) < ncol(flags))
    if (partial) notes <- c(notes,
      "No automatic alternative references: partial study-effect constraints require a model-specific reparameterization.") else {
      ranked <- which(!tab$flagged & seq_len(n) != n)
      ranked <- ranked[order(-tab$n_subjects[ranked], -visits[ranked], ranked)]
      if (length(ranked)) tab$reference_rank[ranked] <- seq_along(ranked)
      chosen <- head(ranked, 3L)
      selected_references <- ids[chosen]
      tab$nominated_reference[chosen] <- TRUE
    }
  }
  notes <- c(notes, "Nominations are descriptive sensitivity checks, not proof of study-specific nonidentification or grounds for deleting a study. The reference study is not symmetrically diagnosed by its own dummy coefficient.")
  list(candidateStudies = if (auto_candidates) selected_candidates else candidates,
    referenceCandidates = if (auto_references) selected_references else references,
    autoCandidates = auto_candidates, autoReferences = auto_references,
    studyDiagnostics = tab, notes = notes)
}
