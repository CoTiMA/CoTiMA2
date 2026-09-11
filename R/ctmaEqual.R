#' ctmaEqual
#'
#' @description test if the two or more invariant drift parameters in the CoTiMAFit object supplied are equal. The supplied CoTiMA fit-object
#' (ctmaInvariantFit) has to be a model fitted with \code{\link{ctmaFit}} where at least two parameters were set invariant across primary studies (e.g., 2 cross
#' effects). All parameters that are set invariant in the supplied model are then constrained to be equal by ctmaEqual (no user action
#' required), the model is fitted, and a log-liklihood ratio test is performed informing about the probability that equality applies.
#'
#' @param ctmaInvariantFit  object to which a CoTiMA fit has been assigned to (i.e., what has been returned by \code{\link{ctmaFit}}).
#' In most cases probably a model in which (only) two effects were specified with invariantDrift.
#' @param activeDirectory defines another active directory than the one used in ctmaInvariantFit
#' @param activateRPB  set to TRUE to receive push messages with CoTiMA notifications on your phone
#' @param digits Number of digits used for rounding (in outputs)
#' @param coresToUse If neg., the value is subtracted from available cores, else value = cores to use
#' @param backend estimation backend, \code{"stan"} (default) or
#'   \code{"julia"}.
#' @param backendArgs named list of controls for the selected backend.
#'
#' @importFrom  RPushbullet pbPost
#' @importFrom  parallel detectCores
#' @importFrom  ctsem ctFit
#' @importFrom  OpenMx vech2full
#'
#' @export ctmaEqual
#'
#' @examples
#' # Fit a CoTiMA with a set of parameters set equal that were set
#' # invariant in a previous model (of which the fit object is
#' # supplied in argument ctmaInvariantFit)
#' \dontrun{
#' CoTiMAFullInv23Fit_6$activeDirectory <- "/Users/tmp/" # adapt!
#' CoTiMAFullInvEq23Fit_6 <- ctmaEqual(ctmaInvariantFit=CoTiMAFullInv23Fit_6)
#' }
#'
#' @return returns a model where two or more parameters were set equal across primary studies and a log-likelihood difference test
#' informing about the probability that the equality assumption is correct.
#'
ctmaEqual <- function(
    ctmaInvariantFit=NULL,
    activeDirectory=NULL,
    activateRPB=FALSE,
    digits=4,
    coresToUse=2,
    backend=c("stan", "julia"),
    backendArgs=list()
)

{  # begin function definition (until end of file)
  .ctmaVersionNote()
  .ctma_timer <- .ctmaRegisterExecutionTimer("ctmaEqual")
  on.exit(.ctma_timer(), add = TRUE)

  .ctmaCheckCtsemJuliaRequest(backend, "ctmaEqual")

  # check if fit object is specified
  if (is.null(ctmaInvariantFit)){
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    ErrorMsg <- "\nA fitted CoTiMA object with more than a single invariant drift effect (fit of ctmaFit) has to be supplied compare the effects. \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  julia_input <- inherits(ctmaInvariantFit, "ctmaJuliaFit")
  supported_model_type <- length(ctmaInvariantFit$model.type) == 1L &&
    !is.na(ctmaInvariantFit$model.type) &&
    ctmaInvariantFit$model.type %in% c("mx", "stanct", "ct")
  # check if mutipleDriftFit object is supplied
  if (!julia_input && !supported_model_type) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    ErrorMsg <- "\nA fitted CoTiMA object with more than a single invariant drift effect (fit of ctmaFit) has to be supplied compare the effects. \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  inherited_stan_args <- ctmaInvariantFit$argumentList$CoTiMAStanctArgs
  inherited_estimation <- ctmaInvariantFit$argumentList$estimation
  if (is.null(inherited_estimation)) {
    inherited_optimize <- inherited_stan_args$optimize
    inherited_estimation <- if (isFALSE(inherited_optimize)) "Bayes" else "ML"
  }
  selected_backend <- match.arg(tolower(backend), c("stan", "julia"))
  backend_resolution <- .ctmaResolveFitBackend(
    backend = selected_backend,
    backendArgs = backendArgs,
    legacy_stan_args = if (identical(selected_backend, "stan") &&
      !length(backendArgs)) inherited_stan_args else NULL,
    estimation = inherited_estimation,
    workflow = "ctmaEqual"
  )
  backend <- backend_resolution$backend
  requested_estimation <- backend_resolution$estimation
  resolved_backend_args <- backend_resolution$backend_args
  CoTiMAStanctArgs <- if (identical(backend, "stan")) {
    .ctmaResolveAllInvStanArgs(overrides = backend_resolution$stan_args)
  } else {
    .ctmaResolveAllInvStanArgs(
      optimize = identical(requested_estimation, "ml"),
      priors = identical(requested_estimation, "bayes")
    )
  }

  input_model <- if (julia_input) {
    ctmaInvariantFit$source_object$ctModel
  } else {
    candidate <- ctmaInvariantFit$studyFitList[[1]]$ctstanmodelbase
    if (is.null(candidate)) candidate <- ctmaInvariantFit$studyFitList$ctstanmodelbase
    candidate
  }
  if (is.null(input_model)) {
    stop("The supplied fit does not contain its ctsem model specification.", call. = FALSE)
  }
  drift_rows_all <- which(input_model$pars$matrix == "DRIFT")
  ti_effect_columns <- paste0(input_model$TIpredNames, "_effect")
  invariant_rows <- if (length(ti_effect_columns)) {
    drift_rows_all[apply(
      input_model$pars[drift_rows_all, ti_effect_columns, drop = FALSE],
      1L, function(x) all(!as.logical(x))
    )]
  } else drift_rows_all
  invariant_names <- as.character(input_model$pars$param[invariant_rows])
  if (length(invariant_names) < 2L) {
    if (activateRPB==TRUE) {RPushbullet::pbPost("note", paste0("CoTiMA (",Sys.time(),")" ), paste0(Sys.info()[[4]], "\n","Data processing stopped.\nYour attention is required."))}
    ErrorMsg <- "\nA fitted CoTiMA object was supplied, but is has to have more than a single invariant drift effect to compare the effects. \nGood luck for the next try!"
    stop(ErrorMsg)
  }



  #######################################################################################################################
  ################################################# Check Cores To Use ##################################################
  #######################################################################################################################
  resource_resolution <- .ctmaResolveCores(
    coresToUse, backend = backend, workflow = "ctmaEqual"
  )
  coresToUse <- resource_resolution$resolved
  if (identical(backend, "julia")) {
    julia_resources <- .ctmaInitJuliaResources(
      coresToUse, requested_estimation, backendArgs, resolved_backend_args,
      workflow = "ctmaEqual"
    )
    resolved_backend_args <- julia_resources$backend_args
  }



  #######################################################################################################################
  ############## Extracting parameters from fitted primary studies created with ctmaFit Function  #######################
  ############### and estimating modified model in which invariant (fixed) effects are set equal. #######################
  #######################################################################################################################

  start.time <- Sys.time(); start.time
  {
    n.latent <- input_model$n.latent; n.latent
    if (is.null(activeDirectory)) {
      activeDirectory <- ctmaInvariantFit$activeDirectory
      if (is.null(activeDirectory)) {
        activeDirectory <- ctmaInvariantFit$source_object$argumentList$activeDirectory
      }
    }
    n.studies <- if (julia_input) input_model$n.TIpred + 1L else unlist(ctmaInvariantFit$n.studies)
    statisticsList <- ctmaInvariantFit$statisticsList
    studyList <- ctmaInvariantFit$studyList
    if (julia_input) {
      if (is.null(studyList)) studyList <- ctmaInvariantFit$source_object$studyList
      parameterNames <- list(
        DRIFT = as.character(input_model$pars$param[drift_rows_all])
      )
    } else parameterNames <- ctmaInvariantFit$parameterNames
    driftNames <- parameterNames$DRIFT
    targetNames <- invariant_names
    # CHD 30.11.2024
    randomIntercepts <- ctmaInvariantFit$argumentList$randomIntercepts
    if (julia_input && is.null(randomIntercepts)) {
      randomIntercepts <- ctmaInvariantFit$source_object$argumentList$randomIntercepts
    }
    if (is.null(randomIntercepts)) randomIntercepts <- FALSE
  }

  # copy previous model
  prevStanctModel <- input_model

  # identify Drift coefficents that were fixed (across all TI, which is just a check)
  tmpRow <- which(prevStanctModel$pars$matrix == "DRIFT"); tmpRow
  # CHD 28.11.2024
  if ( (randomIntercepts != "MANIFEST") & (randomIntercepts != "CINT") ) {
    tmpRow <- invariant_rows
    equalDriftPos <- match(tmpRow, drift_rows_all)
  }
  if ( (randomIntercepts == "MANIFEST") | (randomIntercepts == "CINT") ) {
    equalDriftPos <- grep("DRIFT", prevStanctModel$pars$matrix); equalDriftPos
    equalDriftPos <- equalDriftPos[which(is.na(prevStanctModel$pars$value)[equalDriftPos])]; equalDriftPos
    tmpRow <- equalDriftPos[which(prevStanctModel$pars$TI1_effect[equalDriftPos] == FALSE)]; tmpRow
    equalDriftPos <- which(equalDriftPos %in% tmpRow); equalDriftPos
  }
  tmp1 <- prevStanctModel$pars[tmpRow, paste0(prevStanctModel$TIpredNames,'_effect')]; tmp1
  tmp2 <- apply(tmp1, 1, unique); tmp2
  tmp3 <- unlist(lapply(tmp2, length)); tmp3
  targetDriftRow <- tmpRow[tmp3==1 & tmp2==FALSE]; targetDriftRow

  # new model
  stanctModel <- prevStanctModel
  newDriftLabel <- paste(stanctModel$pars[targetDriftRow, "param"], collapse = "_eq_"); newDriftLabel
  stanctModel$pars[targetDriftRow, "param"] <- newDriftLabel

  previous_fit <- if (julia_input) ctmaInvariantFit else ctmaInvariantFit$studyFitList[[1]]
  if (is.null(previous_fit)) previous_fit <- ctmaInvariantFit$studyFitList
  prevEst <- if (inherits(previous_fit, "ctmaJuliaFit")) {
    as.numeric(coef(previous_fit, type = "raw"))
  } else {
    previous_fit$stanfit$rawest
  }
  # correct for having appropriate inits
  tmp1 <- which(tmpRow %in% targetDriftRow); tmp1
  newInits <- mean(prevEst[tmp1]); newInits
  prevEst[tmp1] <- newInits; prevEst
  prevEst <- prevEst[-tmp1[-1]]; prevEst
  # shortcut
  if ( (randomIntercepts == "MANIFEST") | (randomIntercepts == "CINT") ) prevEst <- NULL


  prevData <- ctmaInvariantFit$data
  if (julia_input && is.null(prevData)) prevData <- ctmaInvariantFit$source_object$data

  priors <- FALSE
  if (!(is.null(CoTiMAStanctArgs$priors))) priors <- CoTiMAStanctArgs$priors else priors <- FALSE

  CoTiMAStanctArgs$optimcontrol$bootstrapUncertainty <- NULL

  if (identical(backend, "julia")) {
    fitStanctModel <- .ctmaInitJuliaStudyFit(
      datalong = prevData,
      model = stanctModel,
      same_initial_times = FALSE,
      stan_args = CoTiMAStanctArgs,
      backend_args = resolved_backend_args,
      estimation = requested_estimation,
      start_values = if (julia_input) NULL else prevEst,
      digits = digits
    )
  } else fitStanctModel <- ctsem::ctFit(
    inits=prevEst,
    datalong = prevData,
    ctstanmodel = stanctModel,
    savesubjectmatrices=CoTiMAStanctArgs$savesubjectmatrices,
    stanmodeltext=CoTiMAStanctArgs$stanmodeltext,
    iter=CoTiMAStanctArgs$iter,
    intoverstates=CoTiMAStanctArgs$intoverstates,
    binomial=CoTiMAStanctArgs$binomial,
    fit=CoTiMAStanctArgs$fit,
    intoverpop=CoTiMAStanctArgs$intoverpop,
    stationary=CoTiMAStanctArgs$stationary,
    plot=CoTiMAStanctArgs$plot,
    #derrind=CoTiMAStanctArgs$derrind, # CHD deprecated, deleted Aug 2023
    optimize=CoTiMAStanctArgs$optimize,
    optimcontrol=CoTiMAStanctArgs$optimcontrol,
    nlcontrol=CoTiMAStanctArgs$nlcontrol,
    #nopriors=CoTiMAStanctArgs$nopriors,
    #nopriors=TRUE,
    #priors=CoTiMAStanctArgs$priors, # added Aug 2023
    priors=priors,
    chains=CoTiMAStanctArgs$chains,
    forcerecompile=CoTiMAStanctArgs$forcerecompile,
    savescores=CoTiMAStanctArgs$savescores,
    gendata=CoTiMAStanctArgs$gendata,
    control=CoTiMAStanctArgs$control,
    verbose=CoTiMAStanctArgs$verbose,
    cores=coresToUse)

  if (identical(backend, "stan") && !(is.null(CoTiMAStanctArgs$resample))) {
    fitStanctModel <- ctmaStanResample(ctmaFittedModel=fitStanctModel) #, CoTiMAStanctArgs=CoTiMAStanctArgs)
  }
  equalDriftStanctFit <- if (inherits(fitStanctModel, "ctmaJuliaFit")) {
    fitStanctModel$ctsem_summary
  } else {
    summary(fitStanctModel, digits=digits)
  }

  # Extract estimates & statistics
  # account for changes in ctsem 3.4.1
  if ("matrix" %in% colnames(equalDriftStanctFit$parmatrices)) ctsem341 <- TRUE else ctsem341 <- FALSE
  tmpMean <- grep("ean", colnames(equalDriftStanctFit$parmatrices)); tmpMean
  tmpSd <- tmpMean+1; tmpSd
  Tvalues <- equalDriftStanctFit$parmatrices[,tmpMean]/equalDriftStanctFit$parmatrices[,tmpSd]; Tvalues
  equalDrift_Coeff <- cbind(equalDriftStanctFit$parmatrices, Tvalues); equalDrift_Coeff
  equalDrift_Coeff[, tmpMean:(dim(equalDrift_Coeff)[2])] <- round(equalDrift_Coeff[, tmpMean:(dim(equalDrift_Coeff)[2])], digits); equalDrift_Coeff

  # re-label
  if (ctsem341) {
    tmp1 <- which(equalDrift_Coeff[, "matrix"] == "DRIFT")
    driftNamesTmp <- c(matrix(driftNames, n.latent, n.latent, byrow=TRUE)); driftNamesTmp
    rownames(equalDrift_Coeff) <- paste0(equalDrift_Coeff[, c("matrix")], "_",
                                         equalDrift_Coeff[, c("row")], "_",
                                         equalDrift_Coeff[, c("col")])
  } else {
    tmp1 <- which(rownames(equalDrift_Coeff) == "DRIFT"); tmp1
    driftNamesTmp <- c(matrix(driftNames, n.latent, n.latent, byrow=FALSE)); driftNamesTmp
  }

  driftNamesTmpBackup <- driftNamesTmp
  driftNamesTmp[equalDriftPos] <- paste0(driftNamesTmp[equalDriftPos], " (invariant & equal)"); driftNamesTmp
  if ( (randomIntercepts != "MANIFEST") & (randomIntercepts != "CINT") ) {
    # A shared equality parameter can be reported once for each affected
    # matrix cell. Keep data-frame row names unique while retaining the common
    # parameter label in the displayed text.
    rownames(equalDrift_Coeff)[tmp1] <- make.unique(driftNamesTmp, sep = "__cell")
  }

  if ( (randomIntercepts == "MANIFEST") | (randomIntercepts == "CINT") ) {
    tmp1 <- tmp2 <- grep("DRIFT_", rownames(equalDrift_Coeff)); tmp1
    allDriftNames <- rownames(equalDrift_Coeff)[tmp1]; allDriftNames
    tmp1 <- grep("dt", allDriftNames); tmp1
    allDriftNames2 <- allDriftNames[-tmp1]; allDriftNames2
    driftNamesTmp2 <- gsub(" \\(invariant & equal\\)", "", driftNamesTmp)
    driftNamesTmp2 <- gsub("V", "", driftNamesTmp2)
    driftNamesTmp2 <- gsub("to", "_", driftNamesTmp2)
    tmp1 <- c()
    for (z in 1:length(driftNamesTmp2)) tmp1 <- c(tmp1, grep(driftNamesTmp2[z], allDriftNames2))
    allDriftNames[tmp1] <- driftNamesTmp
    rownames(equalDrift_Coeff)[tmp2] <- allDriftNames
    rm(tmp1); rm(tmp2)
  }

  tmp2 <- grep("toV", rownames(equalDrift_Coeff)); tmp2
  rownames(equalDrift_Coeff)[tmp2] <- paste0("DRIFT ", rownames(equalDrift_Coeff)[tmp2]); equalDrift_Coeff


  #tmp1 <- which(rownames(equalDrift_Coeff) == "DRIFT"); tmp1
  #rownames(equalDrift_Coeff)[tmp1] <- "DRIFT " # space at the end for later extraction
  #tmp2 <- tmp1[equalDriftPos]; tmp2
  #tmp3 <- paste0(rownames(equalDrift_Coeff)[tmp2] , "(invariant & equal)"); tmp3
  ##tmp1 <- which(rownames(equalDrift_Coeff) %in% newDriftLabel); tmp1
  ##tmp2 <- paste0(rownames(equalDrift_Coeff)[tmp1] , " (invariant)"); tmp2
  ##rownames(equalDrift_Coeff)[tmp1] <- tmp2; equalDrift_Coeff
  #rownames(equalDrift_Coeff)[tmp2] <- tmp3; equalDrift_Coeff
  #
  equalDrift_Minus2LogLikelihood  <- -2*equalDriftStanctFit$loglik; equalDrift_Minus2LogLikelihood
  equalDrift_estimatedParameters  <- equalDriftStanctFit$npars; equalDrift_estimatedParameters
  equalDrift_df <- "deprecated"

  if ( (randomIntercepts != "MANIFEST") & (randomIntercepts != "CINT") ) {
    model_Drift_Coef <- equalDrift_Coeff[grep("DRIFT ", rownames(equalDrift_Coeff)), tmpMean]; model_Drift_Coef
    names(model_Drift_Coef) <- rownames(equalDrift_Coeff)[grep("DRIFT ", rownames(equalDrift_Coeff))]
    names(model_Drift_Coef)[equalDriftPos] <- newDriftLabel; model_Drift_Coef
    names(model_Drift_Coef) <- gsub("DRIFT ", "", names(model_Drift_Coef)); model_Drift_Coef
  }
  if ( (randomIntercepts == "MANIFEST") | (randomIntercepts == "CINT") ) {
    model_Drift_Coef <- equalDrift_Coeff[grep("toV", rownames(equalDrift_Coeff)), tmpMean]; model_Drift_Coef
    names(model_Drift_Coef) <- rownames(equalDrift_Coeff)[grep("toV", rownames(equalDrift_Coeff))]
    #model_Drift_Coef
  }

  model_Diffusion_Coef <- equalDrift_Coeff[grep("DIFFUSIONcov", substr(rownames(equalDrift_Coeff), 1, 12)), tmpMean] ; model_Diffusion_Coef
  if (!(ctsem341)) model_Diffusion_Coef <- c(OpenMx::vech2full(model_Diffusion_Coef)); model_Diffusion_Coef
  names(model_Diffusion_Coef) <- paste0("diff_", driftNames); model_Diffusion_Coef

  if (ctsem341) {
    model_T0var_Coef <- equalDrift_Coeff[grep("T0cov", substr(rownames(equalDrift_Coeff), 1, 5)), tmpMean]; model_T0var_Coef
  } else {
    model_T0var_Coef <- equalDrift_Coeff[grep("T0VAR", substr(rownames(equalDrift_Coeff), 1, 5)), tmpMean]; model_T0var_Coef
    model_T0var_Coef <- c(OpenMx::vech2full(model_T0var_Coef)); model_T0var_Coef
  }
  names(model_T0var_Coef) <- paste0("T0VAR_", driftNames); model_T0var_Coef

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  #st <- paste0("Computation started at: ", start.time); st
  #et <- paste0("Computation ended at: ", end.time); et
  #tt <- paste0("Computation lasted: ", round(time.taken, digits)); tt

  # results is built again after test for equality
  results <- list(activeDirectory=activeDirectory,
                  backend=backend,
                  backendArgs=resolved_backend_args,
                  resource_resolution=resource_resolution,
                  plot.type=NULL, model.type="stanct",
                  #time=list(start.time=start.time, end.time=end.time, time.taken=time.taken),
                  coresToUse=coresToUse, n.studies=n.studies,
                  n.latent=n.latent,
                  studyList=studyList, studyFitList=fitStanctModel,
                  data=prevData, statisticsList=statisticsList,
                  modelResults=list(DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef, CINT=NULL),
                  parameterNames=parameterNames,
                  CoTiMAStanctArgs=if (identical(backend, "stan")) CoTiMAStanctArgs else NULL,
                  argumentList=list(
                    backend=backend,
                    backendArgs=resolved_backend_args,
                    CoTiMAStanctArgs=if (identical(backend, "stan")) CoTiMAStanctArgs else NULL,
                    estimation=requested_estimation,
                    randomIntercepts=randomIntercepts,
                    coresToUse=coresToUse
                  ),
                  equalDrift=newDriftLabel,
                  summary=list(#model=paste(targetNames, collapse=" equal to "),
                    model=newDriftLabel,
                    estimates=equalDrift_Coeff,
                    minus2ll=equalDrift_Minus2LogLikelihood,
                    n.parameters = round(equalDrift_estimatedParameters, digits),
                    scaleTime=if (is.null(ctmaInvariantFit$summary$scaleTime)) 1 else ctmaInvariantFit$summary$scaleTime))
  class(results) <- "CoTiMAFit"

  # model comparison
  comparison_input <- ctmaInvariantFit
  if (is.null(comparison_input$summary$model)) {
    comparison_input$summary$model <- paste(targetNames, collapse = " invariant")
  }
  compFit <- ctmaCompFit(comparison_input, results)[[1]][-7]
  results$summary$ll_diff_test <- compFit

  invisible(results)

} ### END function definition
