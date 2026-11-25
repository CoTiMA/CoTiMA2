#' Convert a CoTiMA fit into a reusable model
#'
#' @param ctmaFitObject A complete CoTiMA initialization or pooled fit, or a ctmaPrep study list.
#' @param model NULL to retain the fitted specification, or a case-insensitive
#'   choice of "init", "ctmaInit", "fit", "ctmaFit", or "prep".
#' @param reUseEmprawData Reuse empirical wide data. Missing caches produce an
#'   error; FALSE permits regeneration from the study covariance matrices.
#' @param emprawLong Additionally include original-time long empirical data.
#' @return A CoTiMAModel for ctmaInit or ctmaFit; for model="prep", the
#'   result of ctmaFitToPrep. Explicit fitting arguments override inherited ones.
#' @details Stored long-data times are returned to original units. The saved
#'   scaleTime is applied once at refitting. Standard-model conversion resets
#'   scientific model settings. Two-wave studies disable inherited person-level
#'   random intercepts in standard initialization models; explicit arguments
#'   can restore them. Original-model conversion retains parameter constraints.
#' @export
ctmaFitToModel <- function(ctmaFitObject, model = NULL,
    reUseEmprawData = TRUE, emprawLong = FALSE) {
  for (flag in list(reUseEmprawData,emprawLong))
    if (!is.logical(flag) || length(flag)!=1L || is.na(flag))
      stop("reUseEmprawData and emprawLong must be TRUE or FALSE.",call.=FALSE)
  if (!is.list(ctmaFitObject)) stop("Supply a complete CoTiMA fit.", call. = FALSE)
  target <- if (is.null(model)) "original" else {
    if (!is.character(model) || length(model) != 1L || is.na(model))
      stop("model must be NULL, 'init', 'fit', or 'prep'.", call. = FALSE)
    sub("^ctma", "", tolower(model))
  }
  if (!target %in% c("original", "init", "fit", "prep"))
    stop("Unknown model target: ", target, call. = FALSE)
  if (target == "prep" && !is.null(ctmaFitObject$empcovs)) return(ctmaFitObject)
  if (target == "prep") return(ctmaFitToPrep(ctmaFitObject,
    reUseEmprawData = reUseEmprawData, emprawLong = emprawLong))
  if (is.list(ctmaFitObject$empcovs) && is.null(ctmaFitObject$studyList)) {
    prep <- ctmaFitObject
    if (target == "original") target <- "fit"
    if (target == "prep") return(prep)
    counts <- vapply(seq_len(prep$n.studies), function(i) {
      ncol(prep$empcovs[[i]]) / (length(prep$deltas[[i]]) + 1L)
    }, numeric(1))
    if (any(!is.finite(counts)) || any(counts < 1) ||
        any(counts != as.integer(counts)) || length(unique(counts)) != 1L)
      stop("Cannot infer n.latent from this prepared study list.",call.=FALSE)
    studyList <- lapply(seq_len(prep$n.studies), function(i) list(
      studyNumber=i, originalStudyNo=prep$studyNumbers[[i]],
      empcov=prep$empcovs[[i]], delta_t=prep$deltas[[i]],
      sampleSize=prep$sampleSizes[[i]], timePoints=length(prep$deltas[[i]])+1L,
      moderators=prep$moderators[[i]], pairwiseN=prep$pairwiseNs[[i]],
      rawData=prep$rawData[[i]], source=prep$source[[i]]))
    if (is.null(prep$emprawList)) {
      .ctma_print_message("Preparing empirical data with ctmaInit(fit=FALSE); primary-study models are not estimated. Set the R seed before conversion for reproducible synthesized data.")
      preparation <- ctmaInit(primaryStudies=prep,n.latent=as.integer(counts[1]),
        activeDirectory=if(is.null(prep$activeDirectory)) tempdir() else prep$activeDirectory,
        fit=FALSE,backend="stan",scaleTime=1,silentOverwrite=TRUE)
      prep$emprawList <- preparation$data
      prep$reUseEmprawData <- TRUE
    }
    long <- if (emprawLong) lapply(prep$emprawList,.ctmaEmprawToLong) else NULL
    return(structure(list(target=target, preserve=FALSE, primaryStudies=prep,
      ctModels=NULL, data=NULL, studyList=studyList,
      emprawList=prep$emprawList,emprawLongList=long,
      argumentList=list(n.latent=as.integer(counts[1]),scaleTime=1,
        activeDirectory=prep$activeDirectory),
      timeScale=list(dataScale="original",scaleTime=1,sourceDataScale="original"),
      adjustments="Prepared-study input: standard models will be generated without fitting primaries.",
      sourceFunction="ctmaPrep", generateEmpirical=is.null(prep$emprawList),
      requestEmprawLong=emprawLong),class="CoTiMAModel"))
  }
  source <- ctmaFitObject
  if (is.null(source$studyList) && !is.null(source$source_object$studyList))
    source <- source$source_object
  prep <- ctmaFitToPrep(source, reUseEmprawData, emprawLong)
  scale <- .ctmaSummaryTimeFactor(ctmaFitObject)
  fits <- ctmaFitObject$studyFitList
  native <- inherits(fits, c("ctJuliaFit", "ctStanFit", "ctmaJuliaFit"))
  is_init <- is.list(fits) && !native && length(fits) == prep$n.studies
  kind <- if (target == "original") if (is_init) "init" else "fit" else target
  get_model <- function(f) {
    m <- f$ctstanmodelbase
    if (is.null(m)) m <- f$modelbase
    if (is.null(m)) m <- f$ctModel
    if (is.null(m) && is.list(f$model) && is.data.frame(f$model$pars)) m <- f$model
    if (is.null(m) || !is.data.frame(m$pars))
      stop("The original model parameter table is unavailable.", call. = FALSE)
    m
  }
  models <- data <- NULL
  if (target == "original") {
    if (is_init) {
      models <- lapply(fits, get_model)
    } else {
      engine <- if (native) fits else ctmaFitObject
      models <- list(get_model(engine))
      data <- ctmaFitObject$data
      if (!is.data.frame(data)) data <- engine$model_spec$data
      if (!is.data.frame(data)) data <- source$data
      if (!is.data.frame(data)) stop("Original pooled long data are unavailable.",call.=FALSE)
      time <- models[[1L]]$timeName
      if (is.null(time)) time <- "time"
      if (!time %in% names(data)) stop("The model time column is missing.",call.=FALSE)
      data[[time]] <- data[[time]] / scale
    }
  }
  inherited <- ctmaFitObject$argumentList
  if (is.null(inherited)) inherited <- source$argumentList
  if (is.null(inherited)) inherited <- list()
  inherited$scaleTime <- scale
  if (is.null(inherited$n.latent)) inherited$n.latent <- ctmaFitObject$n.latent
  if (is.null(inherited$n.latent)) inherited$n.latent <- source$n.latent
  if (is.null(inherited$n.latent)) stop("The substantive n.latent is unavailable.",call.=FALSE)
  adjustments <- character()
  if (target != "original") {
    keep <- c("n.latent", "n.manifest", "scaleTime", "binaries", "activeDirectory")
    inherited <- inherited[intersect(names(inherited), keep)]
    adjustments <- "A standard model will be constructed; original constraints, moderator effects and study-level random effects are not inherited."
    waves <- vapply(prep$deltas, function(x) length(x) + 1L, integer(1))
    if (kind == "init" && any(waves <= 2L)) {
      inherited$indVarying <- FALSE
      adjustments <- c(adjustments, paste0(
        "Person-level random intercepts are disabled in the standard initialization model. Two-wave study IDs: ",
        paste(unlist(prep$studyNumbers)[waves <= 2L], collapse = ", "),
        ". Supply indVarying explicitly to override this setting; identification warnings may apply."))
    }
  }
  for (note in adjustments) .ctma_print_message(note)
  structure(list(target = kind, preserve = target == "original",
    primaryStudies = prep, ctModels = models, data = data,
    argumentList = inherited, studyList = source$studyList,
    emprawList = prep$emprawList, emprawLongList = prep$emprawLongList,
    timeScale = list(dataScale = "original", scaleTime = scale,
      sourceDataScale = "fitted"), adjustments = adjustments,
    sourceFunction = if (is_init) "ctmaInit" else "ctmaFit"), class = "CoTiMAModel")
}

.ctmaModelCallArguments <- function(model, call, frame, fun) {
  explicit <- as.list(call)[-1L]
  explicit$model <- NULL
  values <- lapply(explicit, eval, envir = frame)
  inherited <- model$argumentList
  inherited <- inherited[intersect(names(inherited), names(formals(fun)))]
  inherited[c("model", "ctmaInitFit", "primaryStudies", "primaryStudyList",
    "inits", "loadSingleStudyModelFit", "saveSingleStudyModelFit")] <- NULL
  if ("backend" %in% names(values) && !"backendArgs" %in% names(values) &&
      !identical(as.character(values$backend)[1L],as.character(inherited$backend)[1L])) {
    inherited$backendArgs <- list()
    inherited$CoTiMAStanctArgs <- NULL
    .ctma_print_message("Backend changed; engine-specific controls were reset. Supply backendArgs explicitly to set new controls.")
  }
  # Assignment preserves explicit NULL values, unlike modifyList.
  if ("backendArgs" %in% names(inherited) && !"CoTiMAStanctArgs" %in% names(values))
    inherited$CoTiMAStanctArgs <- NULL
  inherited[names(values)] <- values
  inherited
}

.ctmaRunConvertedModel <- function(model, call, frame, workflow) {
  if (!inherits(model, "CoTiMAModel")) model <- ctmaFitToModel(model)
  if (!identical(model$target, workflow)) stop(
    "Convert explicitly with model = '", workflow,
    "' before changing between pooled and separate-study models.", call. = FALSE)
  fun <- if (workflow == "init") ctmaInit else ctmaFit
  args <- .ctmaModelCallArguments(model, call, frame, fun)
  scale <- args$scaleTime
  if (is.null(scale)) scale <- 1
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) || scale <= 0)
    stop("scaleTime must be one finite positive number.",call.=FALSE)
  if (workflow == "init") {
    prep <- model$primaryStudies
    if (model$preserve) {
      structural <- intersect(names(as.list(call)), c("randomIntercepts", "stationary",
        "determinedDiffusions", "stationaryVariance", "n.latent", "n.manifest", "binaries"))
      changed <- structural[!vapply(structural,function(name)
        isTRUE(all.equal(args[[name]],model$argumentList[[name]])),logical(1))]
      if (length(changed)) stop("Convert with model='init' to rebuild the model structure for: ",
        paste(changed,collapse=", "),call.=FALSE)
      attr(prep, "ctma_models") <- model$ctModels
      attr(prep, "ctma_model_overrides") <- intersect(names(as.list(call)),
        c("indVarying", "indVaryingT0", "manifestMeans", "manifestVars",
          "cint", "T0means", "T0var", "drift", "diff", "lambda"))
    }
    args$primaryStudies <- prep
    return(do.call(ctmaInit, args))
  }
  if (model$preserve) {
    # Preserve model constraints; construction changes require a standard-model conversion.
    changes <- intersect(names(as.list(call)), c("randomIntercepts", "drift", "studyVarying", "manifestMeans", "manifestVars",
      "cint", "T0means", "T0var", "lambda", "stationary", "determinedDiffusions"))
    if (length(changes)) stop("To change pooled model construction, convert with model='fit'. Explicit construction arguments: ",
      paste(changes,collapse=", "),call.=FALSE)
    data <- model$data
    m <- model$ctModels[[1L]]
    explicit <- names(as.list(call))
    if ("indVarying" %in% explicit) {
      selection <- args$indVarying
      if (isTRUE(selection)) selection <- "CINT"
      m$pars$indvarying <- FALSE
      if (is.character(selection) && length(selection)==1L &&
          toupper(selection) %in% c("CINT", "MANIFEST")) {
        family <- if (toupper(selection)=="MANIFEST") "MANIFESTMEANS" else "CINT"
        rows <- which(m$pars$matrix==family)
        if (!length(rows)) stop("The stored model has no ",family," rows.",call.=FALSE)
        fixed <- rows[is.na(m$pars$param[rows])]
        m$pars$param[fixed] <- paste0("MeanV", m$pars$row[fixed], "_",tolower(family))
        m$pars$value[rows] <- NA_real_
        m$pars$indvarying[rows] <- TRUE
      } else if (!is.null(selection) && !identical(selection,FALSE))
        stop("indVarying must be FALSE, NULL, TRUE, 'CINT', or 'MANIFEST'.",call.=FALSE)
    }
    if ("indVaryingT0" %in% explicit) m <- .ctmaApplyT0MeanVariation(m,args$indVaryingT0)
    time <- m$timeName
    if (is.null(time)) time <- "time"
    data[[time]] <- data[[time]] * scale
    setup <- list(data = data, ctModel = m, argumentList = model$argumentList,
      studyList = model$studyList, emprawList = model$emprawList,
      n.latent=model$argumentList$n.latent)
    setup$argumentList$scaleTime <- scale
    if (is.list(setup$argumentList$CoTiMAStanctArgs))
      setup$argumentList$CoTiMAStanctArgs$scaleTime <- scale
    controls <- intersect(names(args), c("backend", "backendArgs", "optimize", "priors",
      "fit", "coresToUse", "chains", "iter", "inits", "verbose"))
    args <- c(list(model = setup), args[controls])
    result <- do.call(ctmaFit, args)
    if (isFALSE(args$fit)) return(result)
    if (inherits(result,c("ctJuliaFit","ctmaJuliaFit")))
      return(.ctmaWrapPooledJuliaFit(result,setup))
    if (inherits(result,"ctStanFit")) {
      container <- setup
      container$studyFitList <- result
      container$summary <- summary(result,parmatrices=TRUE)
      container$n.studies <- 1L
      container$n.latent <- model$argumentList$n.latent
      class(container) <- "CoTiMAFit"
      return(.ctmaStoreOriginalTime(container))
    }
    result
  }
  # A pooled constructor needs study metadata and wide data, but no initial optimization.
  prep <- model$primaryStudies
  k <- prep$n.studies
  n <- as.integer(model$argumentList$n.latent)
  if (is.null(prep$emprawList)) stop(
    "Constructing a pooled model requires cached empirical data; use reUseEmprawData=TRUE.",call.=FALSE)
  manifest_count <- model$argumentList$n.manifest
  if (is.null(manifest_count) || manifest_count <= n) manifest_count <- 0L
  source <- structure(list(n.studies = k, n.latent = n, n.manifest = manifest_count,
    primaryStudyList = prep, studyList = model$studyList, emprawList = prep$emprawList,
    studyFitList = rep(list(NULL),k),
    activeDirectory = prep$activeDirectory,
    modelResults = list(DRIFT = rep(list(rep(0, n*n)), k)),
    statisticsList = list(allSampleSizes = vapply(prep$emprawList,nrow,integer(1)),
      allTpoints = vapply(prep$deltas, function(x) length(x)+1L, integer(1)),
      allDeltas = unlist(prep$deltas), lambda = diag(n)),
    parameterNames = list(DRIFT = as.vector(outer(seq_len(n),seq_len(n),
      function(i,j) paste0("V",j,"toV",i))))), class="CoTiMAFit")
  args$n.latent <- NULL
  args$n.manifest <- NULL
  args$ctmaInitFit <- source
  do.call(ctmaFit, args)
}
