#' ctmaFitToPrep
#'
#' @description Extracts information from fitted CoTiMA objects to (re-)create list of primary studies originally created with \code{\link{ctmaPrep}}
#'
#' @param ctmaFitObject  ctmaFitObject
#' @param reUseEmprawData Whether cached wide data should be reused in subsequent fits.
#' @param emprawLong Logical; additionally return cached empirical data in
#'   \code{emprawLongList}, one long-format data frame per study, with
#'   \code{id}, \code{time}, and manifest columns. Wide data remain in
#'   \code{emprawList}. Times use the original cached intervals, without applying
#'   \code{scaleTime}. Does not enable reuse unless \code{reUseEmprawData=TRUE}.
#'
#' @examples
#' ctmaLoadExamples("CoTiMAInitFit_3")
#' newStudyList <- ctmaFitToPrep(CoTiMAInitFit_3)
#' ctmaUnloadExamples("CoTiMAInitFit_3")
#'
#' @export ctmaFitToPrep
#'
#' @return list that could be used for fitting new CoTiMA models with \code{\link{ctmaInit}} or \code{\link{ctmaFit}}.
#'
ctmaFitToPrep <- function(ctmaFitObject=NULL, reUseEmprawData=FALSE, emprawLong=FALSE) {
  if (!is.list(ctmaFitObject) || !length(ctmaFitObject)) stop(
    "Supply a CoTiMA fit or a nonempty list of study definitions.", call. = FALSE)
  if (!is.logical(reUseEmprawData) || length(reUseEmprawData) != 1L ||
      is.na(reUseEmprawData)) stop("reUseEmprawData must be TRUE or FALSE.", call. = FALSE)
  if (!is.logical(emprawLong) || length(emprawLong) != 1L || is.na(emprawLong))
    stop("emprawLong must be TRUE or FALSE.", call. = FALSE)
  context <- ctmaFitObject
  # Historical direct pooled Julia results retain the CoTiMA setup here.
  if (is.null(context$studyList) && is.null(context$cotima_study) &&
      is.list(context$source_object) && is.list(context$source_object$studyList)) context <- context$source_object
  is_study <- function(x) is.list(x) && !is.null(x[["originalStudyNo", exact=TRUE]])
  entries <- if (is_study(context$cotima_study)) {
    list(context$cotima_study)
  } else if (is_study(context)) {
    list(context)
  } else if (is.list(context$studyList) && length(context$studyList)) {
    context$studyList
  } else if (all(vapply(context, is_study, logical(1)))) {
    context
  } else stop("No primary-study definitions were found. Supply the complete ",
    "ctmaInit object; a bare native fit without study metadata is insufficient.", call. = FALSE)
  ids <- lapply(entries, function(x) x[["originalStudyNo", exact=TRUE]])
  if (!all(vapply(ids, function(x) length(x)==1L && !is.na(x), logical(1))) ||
      anyDuplicated(as.character(unlist(ids)))) stop(
    "Study definitions must have unique, nonmissing originalStudyNo values.", call. = FALSE)
  mapping <- c(studyNumbers="originalStudyNo", empcovs="empcov", deltas="delta_t",
    sampleSizes="sampleSize", moderators="moderators", source="source",
    pairwiseNs="pairwiseN", rawData="rawData", startValues="startValues")
  result <- lapply(mapping, function(field) lapply(entries, function(x) x[[field, exact=TRUE]]))
  names(result) <- names(mapping)
  result$n.studies <- length(entries)
  result$summary <- "Summary is only available from ctmaPrep and not from ctmaFitToPrep."
  if (reUseEmprawData || emprawLong) {
    cached <- context$emprawList
    if (is.null(cached) && length(entries)==1L && !is.null(context$empraw))
      cached <- list(context$empraw)
    fits <- context$studyFitList
    if (is.null(cached) && is.list(fits) && !.ctmaIsNativeFit(fits) &&
        length(fits)==length(entries)) cached <- lapply(fits, function(x) x$empraw)
    if (is.null(cached) && all(vapply(entries,function(x)!is.null(x$empraw),logical(1))))
      cached <- lapply(entries, function(x) x$empraw)
    if (!is.list(cached) || length(cached)!=length(entries) ||
        !all(vapply(cached,function(x)(is.data.frame(x)||is.matrix(x)) &&
          nrow(x)>0L && ncol(x)>0L,logical(1)))) stop(
      "Cached wide-format data are missing or do not match the study definitions. ",
      "Supply the complete ctmaInit object, or use reUseEmprawData=FALSE.", call. = FALSE)
    for(i in seq_along(entries)) .ctmaValidateCachedStudyWide(cached[[i]],entries[[i]])
    result$emprawList <- cached
    # Preserve covariance metadata while explicitly choosing cached data over
    # re-simulation from those covariances in the next ctmaInit call.
    result$reUseEmprawData <- reUseEmprawData
    if (emprawLong) result$emprawLongList <- lapply(cached, .ctmaEmprawToLong)
    if (reUseEmprawData)
    result$rawData <- lapply(seq_along(entries), function(i) {
      spec <- result$rawData[[i]]
      if (!is.list(spec)) spec <- list()
      spec$studyNumbers <- ids[[i]]
      spec
    })
  }
  directory <- context$activeDirectory
  if (is.null(directory)) directory <- context$argumentList$activeDirectory
  if (!is.null(directory)) result$activeDirectory <- directory
  class(result) <- "CoTiMAFit"
  invisible(result)
}

# ctWideToLong selects manifest columns by name but interval columns by position.
# Explicitly order both to avoid lexical ordering errors at T10 and beyond.
.ctmaEmprawToLong <- function(wide) {
  columns <- colnames(wide)
  manifest_columns <- columns[grepl("_T[0-9]+$", columns)]
  manifests <- unique(sub("_T[0-9]+$", "", manifest_columns))
  occasions <- sort(unique(as.integer(sub(".*_T", "", manifest_columns))))
  if (!length(manifests) || length(occasions) < 2L ||
      !identical(occasions, seq.int(0L, max(occasions))))
    stop("Cached wide data require manifest columns ending in _T0, _T1, ... .", call. = FALSE)
  required <- unlist(lapply(occasions, function(t) paste0(manifests, "_T", t)), use.names=FALSE)
  intervals <- paste0("dT", seq_len(length(occasions)-1L))
  if (!all(c(required, intervals) %in% columns))
    stop("Cached wide data lack a complete manifest/time grid or dT interval columns.", call. = FALSE)
  long <- ctsem::ctWideToLong(wide[, c(required, intervals), drop=FALSE],
    Tpoints=length(occasions), n.manifest=length(manifests), manifestNames=manifests)
  as.data.frame(suppressMessages(ctsem::ctDeintervalise(long, id="id", dT="dT")))
}


.ctmaValidateCachedStudyWide <- function(data,study,nvar=NULL) {
  waves <- study$timePoints
  if(is.null(waves) && length(study$delta_t))waves<-length(study$delta_t)+1L
  if(length(waves)!=1L || !is.finite(waves) || waves<2L)return(invisible(TRUE))
  if(is.null(nvar) && is.matrix(study$empcov) && nrow(study$empcov)>0)
    nvar <- nrow(study$empcov)/waves
  if(length(nvar)!=1L || !is.finite(nvar) || nvar<1L || nvar!=as.integer(nvar))return(invisible(TRUE))
  expected <- as.integer(nvar*waves+waves-1L)
  if(!(is.data.frame(data)||is.matrix(data)) || ncol(data)!=expected)stop(
    "Cached wide data for Study ",study$originalStudyNo," do not match its ",waves,
    " waves (expected ",expected," columns). Recreate the prep object from the complete initialization and subset by original study IDs.",call.=FALSE)
  invisible(TRUE)
}
