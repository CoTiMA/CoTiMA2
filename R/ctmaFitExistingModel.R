# Public ctmaFit(model=...) entry, sharing the standard native backend route.
.ctmaFitExistingModel <- function(model, call, backend, backendArgs,
    optimize, priors, fit, cores, chains, iter, inits, verbose) {
  allowed <- c("model", "backend", "backendArgs", "optimize", "priors", "fit",
               "coresToUse", "chains", "iter", "inits", "verbose")
  supplied <- setdiff(names(as.list(call)[-1L]), allowed)
  if (length(supplied)) stop(
    "With model=, construction settings come from the stored model. Remove: ",
    paste(supplied, collapse = ", "), ".", call. = FALSE)
  backend <- match.arg(tolower(backend[1L]), c("stan", "julia", "ctsemjulia"))
  if (!is.list(model) || !is.data.frame(model$data) ||
      !is.list(model$ctModel) || !is.data.frame(model$ctModel$pars)) stop(
    "model must contain pooled long-format data and ctModel with a parameter table.",
    call. = FALSE)
  if (!is.logical(fit) || length(fit) != 1L || is.na(fit))
    stop("fit must be TRUE or FALSE.", call. = FALSE)
  if (!is.logical(optimize) || length(optimize) != 1L || is.na(optimize))
    stop("optimize must be TRUE or FALSE.", call. = FALSE)
  # Discard old compiled representations, but retain the model, data and metadata.
  setup <- model[c("data", "ctModel", "argumentList", "studyList", "emprawList")]
  setup <- setup[!vapply(setup, is.null, logical(1))]
  class(setup) <- "CoTiMAFit"
  if (is.null(setup$argumentList)) setup$argumentList <- list()
  if (!identical(backend, "ctsemjulia")) return(.ctmaFitExistingOtherBackend(
    setup, model$initialFit, backend, backendArgs, optimize, priors,
    fit, cores, chains, iter, inits, verbose, call))
  estimation <- if (optimize) "ML" else "Bayes"
  resolved <- .ctmaResolveFitBackend(backend = "ctsemjulia",
    backendArgs = backendArgs, legacy_stan_args = NULL, estimation = estimation)
  args <- resolved$backend_args
  # Preserve stationarity unless explicitly overridden in backendArgs.
  if (is.null(backendArgs$stationary)) {
    stationary <- setup$argumentList$stationary
    if (is.null(stationary)) stationary <- setup$argumentList$CoTiMAStanctArgs$stationary
    if (!is.null(stationary)) args$stationary <- stationary
  }
  if (!optimize) {
    if (!is.null(chains) && is.null(backendArgs$sampleControl$chains))
      args$sampleControl$chains <- as.integer(chains)
    if (!is.null(iter)) {
      if (is.null(backendArgs$sampleControl$warmup)) args$sampleControl$warmup <- as.integer(iter %/% 2)
      if (is.null(backendArgs$sampleControl$draws)) args$sampleControl$draws <- as.integer(iter - iter %/% 2)
    }
  }
  .ctmaVersionNote()
  .ctmaFitCtsemJuliaBackend(setup, args, cores = cores,
    cotima_call = call, inits = inits, priors = priors, verbose = verbose,
    estimation = estimation, fit = fit)
}

# Convert public family/parameter selections to the existing RE masks.
.ctmaExistingREMasks <- function(model, args) {
  families <- c("DRIFT", "DIFFUSION", "T0VAR")
  requested <- toupper(as.character(args$study_random_effects))
  if (length(setdiff(requested, families))) stop(
    "For this backend, study_random_effects supports DRIFT, DIFFUSION and T0VAR.", call. = FALSE)
  pars <- model$pars
  labels <- as.character(args$study_parameters)
  free <- !is.na(pars$param) & is.na(pars$value)
  available <- unique(as.character(pars$param[free & pars$matrix %in% families]))
  if (length(setdiff(labels, available))) stop(
    "Unknown or unsupported study_parameters: ",
    paste(setdiff(labels, available), collapse = ", "), call. = FALSE)
  n <- as.integer(model$n.latent)
  setNames(lapply(families, function(family) {
    mask <- .ctma_bool_matrix(args[[paste0("re",family)]], n, paste0("re",family))
    rows <- which(free & pars$matrix == family &
      (family %in% requested | pars$param %in% labels))
    if (length(rows)) mask[cbind(pars$row[rows],pars$col[rows])] <- TRUE
    mask
  }), paste0("re",families))
}

.ctmaFitExistingOtherBackend <- function(setup, initial_fit, backend, args,
    optimize, priors, fit, cores, chains, iter, inits, verbose, call) {
  if (is.null(args)) args <- list()
  if (!is.list(args) || (length(args) && (is.null(names(args)) ||
      any(!nzchar(names(args))) || anyDuplicated(names(args))))) stop(
    "backendArgs must be a fully named list without duplicate names.", call. = FALSE)
  common <- c("study_random_effects", "study_parameters", "study_covariance",
    "reDRIFT", "reDIFFUSION", "reT0VAR", "study_id", "study_id_ti",
    "study_dummy_ti", "remove_study_fixed_effects", "rePrep", "reFit")
  selected <- intersect(names(args),common)
  re <- args[selected]
  engine <- args[setdiff(names(args),common)]
  masks <- .ctmaExistingREMasks(setup$ctModel,re)
  has_re <- any(unlist(masks))
  stationary <- setup$argumentList$stationary
  if (is.null(stationary)) stationary <- setup$argumentList$CoTiMAStanctArgs$stationary
  if (is.null(stationary)) stationary <- FALSE
  if (!has_re) {
    if (length(selected)) stop("RE controls were supplied but no study effects were selected.", call.=FALSE)
    if (backend == "stan") {
      managed <- c("datalong","model","backend","fit","optimize","priors","cores","chains","iter","inits")
      if (length(intersect(names(engine),managed))) stop(
        "Supply fitting arguments at the top level, not inside backendArgs: ",
        paste(intersect(names(engine),managed),collapse=", "),call.=FALSE)
      defaults <- list(datalong=setup$data,model=setup$ctModel,fit=fit,
        optimize=optimize,priors=priors,cores=cores,inits=inits,verbose=verbose,
        stationary=stationary,sameInitialTimes=isTRUE(setup$argumentList$sameInitialTimes))
      if ("backend" %in% names(formals(ctsem::ctFit))) defaults$backend <- "stan"
      if (!is.null(chains)) defaults$chains <- chains
      if (!is.null(iter)) defaults$iter <- iter
      return(do.call(ctsem::ctFit,utils::modifyList(defaults,engine)))
    }
    resolved <- .ctmaResolveFitBackend("julia",engine,NULL,
      if (optimize) "ML" else "Bayes")$backend_args
    if (!fit) return(setup)
    if (!is.null(inits)) resolved$inits <- inits
    if (!is.null(chains)) resolved$chains <- chains
    if (!is.null(iter)) {
      if (optimize) resolved$maxiters <- iter else {
        if (is.null(engine$warmup)) resolved$warmup <- iter %/% 2L
        if (is.null(engine$samples)) resolved$samples <- iter - iter %/% 2L
      }
    }
    # The bundled ordinary-fitting path does not implement MAP priors.
    .ctmaCheckBundledPriorRequest("julia",optimize,priors)
    out <- .ctmaJuliaFitFromSetup(setup,if(optimize) "ML" else "Bayes",resolved)
    return(.ctmaFinalizeJuliaCoTiMAFit(out,setup,call,resolved))
  }
  if (!optimize && !isTRUE(priors)) stop(
    "Bayesian study-RE fitting with this backend requires priors=TRUE.",call.=FALSE)
  if (backend == "julia" && optimize) stop(
    "The bundled Julia study-RE backend supports Bayesian estimation only; use optimize=FALSE, priors=TRUE.",call.=FALSE)
  covariance <- if (is.null(re$study_covariance)) "full" else re$study_covariance
  covariance <- match.arg(covariance,c("full","independent","equicorrelated","block"))
  if (backend == "julia" && covariance != "independent") stop(
    "The bundled Julia study-RE backend currently requires study_covariance='independent'.",call.=FALSE)
  prep_options <- re$rePrep
  if (is.null(prep_options)) prep_options <- list()
  if (!is.list(prep_options) || (length(prep_options) && is.null(names(prep_options))))
    stop("backendArgs$rePrep must be a named list.",call.=FALSE)
  allowed <- setdiff(names(formals(ctmaRePrep)),c("model",names(masks),"covariance",
    "study_id","study_id_ti","study_dummy_ti","remove_study_fixed_effects"))
  if (length(setdiff(names(prep_options),allowed))) stop("Unknown or conflicting rePrep controls.",call.=FALSE)
  if (is.null(prep_options$parameterization)) prep_options$parameterization <- "noncentered"
  if (backend == "julia" && prep_options$parameterization != "noncentered") stop(
    "The bundled Julia study-RE backend requires noncentered parameterization.",call.=FALSE)
  if (length(engine) && backend == "stan") stop(
    "For Stan study-RE fits use backendArgs$rePrep and backendArgs$reFit for additional controls.",call.=FALSE)
  if (is.null(setup$argumentList$CoTiMAStanctArgs)) setup$argumentList$CoTiMAStanctArgs <- list()
  setup$argumentList$CoTiMAStanctArgs$stationary <- stationary
  generated <- .ctma_generate_ctsem_setup(setup, verbose=verbose>0)
  if (!is.null(initial_fit)) {
    # Only reuse raw estimates when the generated parameterization agrees.
    fields <- c("matrix","row","col","param","value","transform","indvarying")
    old <- initial_fit$ctstanmodelbase$pars
    new <- generated$ctstanmodelbase$pars
    dims <- c("nparams","nindvarying","nindvaryingoffdiagonals","ntipredeffects","nsubjects")
    if (is.null(old) || !isTRUE(all.equal(old[,fields],new[,fields],check.attributes=FALSE)) ||
        !identical(initial_fit$standata[dims],generated$standata[dims])) stop(
      "model$initialFit has an incompatible raw parameterization. Omit it or refit the pooled model with the current ctsem version.",call.=FALSE)
    generated$stanfit <- initial_fit$stanfit
  }
  prep_call <- c(list(model=generated,covariance=covariance),masks,
    re[intersect(names(re),c("study_id","study_id_ti","study_dummy_ti","remove_study_fixed_effects"))],prep_options)
  prepared <- do.call(ctmaRePrep,prep_call)
  prepared$CoTiMASetup <- setup
  if (backend == "julia" && any(prepared$specification$type != 1L)) stop(
    "The bundled Julia backend does not support these covariance-cell study effects.",call.=FALSE)
  if (!fit) return(prepared)
  if (backend == "julia" && is.null(prepared$init_template)) stop(
    "The bundled Julia study-RE backend needs fitted starting values. Supply a compatible ctsem pooled fit as model$initialFit, or use backend='stan' or 'ctsemjulia'.",call.=FALSE)
  fit_options <- re$reFit
  if (is.null(fit_options)) fit_options <- list()
  allowed_fit <- setdiff(names(formals(ctmaReFit)),c("prepared","backend","backendArgs","estimation","method","cores","chains","iter","init"))
  if (!is.list(fit_options) || (length(fit_options) && is.null(names(fit_options))) ||
      length(setdiff(names(fit_options),allowed_fit))) stop("Unknown or conflicting reFit controls.",call.=FALSE)
  if (is.null(chains)) chains <- 4L
  if (is.null(iter)) iter <- 2000L
  # Validate engine controls early, even for the bundled Julia RE path.
  if (backend == "julia") engine <- .ctmaMergeBackendArgs(
    ctmaBackendDefaults("julia","Bayes","ctmaReFit"),engine,"backendArgs")
  out <- do.call(ctmaReFit,c(list(prepared=prepared,backend=backend,
    estimation=if(!optimize) "bayes" else if(priors) "map" else "ml",
    cores=cores,chains=chains,iter=iter,init=inits,backendArgs=engine),fit_options))
  out$source_object <- setup
  out$cotima_call <- call
  out
}
