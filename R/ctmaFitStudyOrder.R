# Select and order initialization records atomically by their original study IDs.
# Raw TI coefficients are tied to this order; selecting the same members alone
# is not sufficient to reuse an optimizer's raw vector.
.ctmaFitStudyOrder <- function(init, requested) {
  studies <- init$studyList
  ids <- vapply(studies, function(s) as.character(s$originalStudyNo), character(1))
  target <- as.character(unlist(requested, use.names = FALSE))
  if (!length(ids) || anyNA(ids) || anyDuplicated(ids) ||
      !length(target) || anyNA(target) || anyDuplicated(target) ||
      any(!nzchar(target)) || !all(target %in% ids))
    stop("primaryStudyList must contain unique, nonmissing original study IDs present in ctmaInitFit.", call. = FALSE)
  ix <- match(target, ids); n <- length(ids)
  take <- function(x, name) {
    if (!length(x)) return(x)
    if (length(x) != n) stop("Initialization field '", name,
      "' is not aligned with studyList; regenerate the initialization object.", call. = FALSE)
    x[ix]
  }
  init$studyList <- studies[ix]
  for (i in seq_along(ix)) init$studyList[[i]]$studyNumber <- init$studyList[[i]]$originalStudyNo
  for (name in intersect(c("studyFitList", "emprawList", "ind.mod.List"), names(init)))
    init[name] <- list(take(init[[name]], name))
  for (name in names(init$modelResults)) {
    x <- init$modelResults[[name]]
    if (is.list(x) && !is.data.frame(x))
      init$modelResults[name] <- list(take(x, paste0("modelResults$", name)))
  }
  # A historical initialization may retain a larger prep than studyList.
  # Match its independent original IDs instead of indexing it positionally.
  prep <- init$primaryStudyList
  if (!is.null(prep)) {
    prepids <- as.character(unlist(prep$studyNumbers, use.names = FALSE))
    if (anyNA(prepids) || anyDuplicated(prepids) || !all(target %in% prepids))
      stop("Initialization preparation has inconsistent original study IDs.", call. = FALSE)
    j <- match(target, prepids)
    indexed <- c("studyNumbers", "empcovs", "deltas", "sampleSizes", "moderators",
      "source", "pairwiseNs", "rawData", "startValues", "emprawList", "emprawLongList",
      "empMeans", "empVars", "ageM", "malePercent", "occupation", "country", "alphas",
      "targetVariables", "recodeVariables", "combineVariables", "combineVariablesNames",
      "missingVariables", "inits")
    for (name in intersect(indexed, names(prep))) {
      x <- prep[[name]]
      if (!length(x)) next
      if (length(x) != length(prepids)) stop("Initialization prep field '", name,
        "' is not aligned with its studyNumbers.", call. = FALSE)
      prep[name] <- list(x[j])
    }
    prep$n.studies <- length(target)
    # A workbook is a mutable reference object. Do not reorder it in place.
    prep$excelSheets <- NULL
    if (is.data.frame(prep$summary) && nrow(prep$summary) == length(prepids))
      prep$summary <- prep$summary[j, , drop = FALSE]
    init$primaryStudyList <- prep
  }
  # Summaries and wave diagnostics have their own row identities.
  for (name in c("estimates", "confidenceIntervals")) {
    x <- init$summary[[name]]
    if (!is.null(dim(x))) {
      rowids <- rownames(x)
      j <- if (is.null(rowids) && nrow(x) == n) ix else
        match(target, sub("^Study No ", "", rowids))
      if (anyNA(j)) stop("Initialization summary rows cannot be matched to original study IDs.", call. = FALSE)
      init$summary[name] <- list(x[j, , drop = FALSE])
    }
  }
  if (is.data.frame(init$studyDesign) && "study" %in% names(init$studyDesign))
    init$studyDesign <- init$studyDesign[match(target, as.character(init$studyDesign$study)), , drop = FALSE]
  stats <- init$statisticsList
  if (!is.null(stats)) {
    stats$originalStudyNumbers <- unlist(lapply(init$studyList, `[[`, "originalStudyNo"), use.names = FALSE)
    stats$allSampleSizes <- unlist(lapply(init$studyList, `[[`, "sampleSize"), use.names = FALSE)
    stats$allTpoints <- unlist(lapply(init$studyList, function(s)
      if (!is.null(s$timePoints)) s$timePoints else length(s$delta_t) + 1L), use.names = FALSE)
    stats$allDeltas <- unlist(lapply(init$studyList, `[[`, "delta_t"), use.names = FALSE)
    for (pair in list(c("SampleSize", "allSampleSizes"), c("Tpoints", "allTpoints"), c("Delta", "allDeltas"))) {
      x <- stats[[pair[2]]]
      if (length(x) && any(is.finite(x))) {
        for (fun in c("mean", "min", "max")) stats[[paste0(fun, pair[1])]] <- match.fun(fun)(x, na.rm = TRUE)
        if (pair[1] != "Delta") stats[[paste0("overall", pair[1])]] <- sum(x, na.rm = TRUE)
      }
    }
    stats$n.studies <- length(target)
    init$statisticsList <- stats
  }
  init$n.studies <- length(target)
  init
}

# ctmaPrep historically appends one empty NA-study record. Retain that format:
# several legacy consumers explicitly expect it. Normalize only at this boundary,
# and only when the study count and all study-indexed payloads confirm emptiness.
.ctmaPrepSelectionIDs <- function(prep) {
  if (!is.list(prep) || is.null(prep$studyNumbers))
    stop("primaryStudyList must contain studyNumbers.",call.=FALSE)
  ids <- unlist(prep$studyNumbers,use.names=FALSE)
  if (!anyNA(ids)) return(ids)
  last <- length(ids)
  if (last<2L || !is.na(ids[last]) || anyNA(ids[-last]) ||
      length(prep$n.studies)!=1L || is.na(prep$n.studies) ||
      prep$n.studies!=last-1L)
    stop("primaryStudyList contains invalid missing study IDs, not a confirmed empty closing row.",call.=FALSE)
  empty <- function(x) {
    if(is.null(x)||!length(x))return(TRUE)
    if(is.list(x))return(all(vapply(x,empty,logical(1))))
    if(is.matrix(x)||is.data.frame(x))return(FALSE)
    is.atomic(x)&&all(is.na(x))
  }
  metadata <- c('n.studies','summary','excelSheets','activeDirectory','plot.type',
    'moderatorLabels','moderatorValues','reUseEmprawData')
  for(name in setdiff(names(prep),metadata)) {
    x <- prep[[name]]
    if(!is.list(x)||is.data.frame(x))next
    if(length(x)!=last) {
      # Known indexed fields cannot have inconsistent lengths; global metadata
      # are listed above and are never inspected as study rows.
      if(name%in%c('studyNumbers','empcovs','deltas','sampleSizes','moderators',
          'source','pairwiseNs','rawData','startValues','emprawList','emprawLongList'))
        stop("primaryStudyList field '",name,"' is misaligned.",call.=FALSE)
      next
    }
    value<-x[[last]]
    if(name=='rawData'&&is.list(value)) {
      if(!empty(value$fileName)||!empty(value$studyNumbers))
        stop("The missing-ID closing row contains raw-data information.",call.=FALSE)
    } else if(!empty(value)) stop("The missing-ID closing row contains data in '",name,"'.",call.=FALSE)
  }
  ids[-last]
}
