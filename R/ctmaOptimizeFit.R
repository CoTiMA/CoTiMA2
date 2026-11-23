#' ctmaOptimizeFit
#'
#' @description Replaces deprecated \code{\link{ctmaOptimizeInit}}, which was limited to initial fitting
#' (i.e., applies \code{\link{ctmaInit}}) of a primary study reFits times to capitalize on chance for obtaining
#' a hard-to-find optimal fit.
#' Now, optimizing a CoTiMA model generated with \code{\link{ctmaFit}} can also be done.
#' Using \code{\link{ctmaOptimizeFit}} could be helpful if a model yields out-of-range estimates, which could happen if the fitting
#' algorithm unfortunately used random start values that resulted in a locally but not globally optimal fit. Essentially, using
#' \code{\link{ctmaOptimizeFit}} is like gambling, hoping that at least one set of starting values (the number it tries is specified in the reFits argument)
#' enables finding the global optimal fit.
#'
#' @details Study shuffles use original study identities and keep study-indexed preparation fields aligned. When no per-study wide-data cache is retained by the best pooled fit, usedStudyList contains study definitions with reUseEmprawData set FALSE; no empirical cache is invented.
#'
#' @param activeDirectory activeDirectory
#' @param coresToUse if neg., the value is subtracted from available cores, else value = cores to use
#' @param CoTiMAStanctArgs parameters that can be set to improve model fitting of the \code{\link[ctsem]{ctStanFit}} Function
#' @param ctmaFitFit a object fitted with \code{\link{ctmaFit}}
#' @param ctmaInitFit the ctmaInitFit object that was used to create the ctmaFitFit object with \code{\link{ctmaFit}}
#' @param ctStanFit a fit object created with ctStanFit (default=NULL)
#' @param ctStanFitArgs list of arguments passed forward to ctStanFit except datalong, ctstanmodel, cores & verbose  (default=NULL)
#' @param ctStanData data to be used for ctStanFit  (default=NULL)
#' @param customPar logical. If set TRUE leverages the first pass using priors and ensure that the drift diagonal cannot easily go too negative (helps since ctsem > 3.4)
#' @param finishsamples number of samples to draw (either from hessian based covariance or posterior distribution) for final results computation (default = 1000).
#' @param iter number of iterations (default = 5000)
#' @param primaryStudies list of primary study information created with \code{\link{ctmaPrep}} or \code{\link{ctmaFitToPrep}}
#' @param problemStudy number (position in list) where the problem study in primaryStudies is found
#' @param randomIV logical (default = FALSE). randomly varies between "indVarying='MANIFEST'" and "indVarying='CINT'"
#' @param randomPar logical (default = FALSE). Overrides arguments used for customPar and randomly sets customPar either TRUE or FALSE
#' @param randomRI logical (default = FALSE). randomly varies between "randomIntercepts='MANIFEST'" and "randomIntercepts='CINT'"
#' @param randomScaleTime two finite positive numbers giving the lower and
#'   upper limits (default = c(1,1)) of the uniform distribution from which
#'   \code{scaleTime} is drawn for each attempt. Draws are rounded to two
#'   decimal places.
#' @param randomScaleTI logical (default = FALSE). If TRUE, randomly sets
#'   \code{scaleTI} to TRUE or FALSE for every \code{ctmaInit()} or
#'   \code{ctmaFit()} attempt.
#' @param reFits positive integer giving the number of complete refitting
#'   attempts.
#' @param saveModelFits save the fit of each Fit attempt (default = FALSE).
#' @param scaleTI scale TI predictors - not recommended until version 0.5.3.1. Does not change aggregated results anyways, just interpretation of effects for dummies representing primary studies.
#' @param scaleTime scale time (interval) - sometimes desirable to improve fitting
#' @param shuffleStudyList (default = FALSE) randomly re-arranges studies in primaryStudyList. We encountered a few cases where this mattered, even though it should not. Only works if ctmaFit is optimized.
#' @param verbose integer from 0 to 2. Higher values print more information during model fit – for debugging
#' @param backend Estimation backend used for complete-model refits. Stan is
#'   the default. Julia is currently available for the \code{ctmaFitFit}
#'   branch and ML estimation only.
#' @param backendArgs Named list of controls for the selected backend. Do not
#'   supply this together with \code{CoTiMAStanctArgs}. Calling
#'   \code{ctmaBackendDefaults()} without arguments returns Stan defaults;
#'   specify the backend and workflow to inspect another supported list.
#'
#' @param randomDeterminedDiffusions Logical, default FALSE. With gambling=TRUE,
#'   randomly choose determinedDiffusions=FALSE or TRUE independently for each
#'   pooled ML refit. Requires a baseline explicitly fitted with stationary=TRUE,
#'   its ctmaInitFit, and backend stan or ctsemjulia. Both variants are checked
#'   with fit=FALSE before fitting; unsupported models are rejected. Model,
#'   time scale, study order and stationarity remain fixed. Old raw starts are
#'   discarded. seed controls reproducible choices and attempts records all
#'   selected variants, failures and objective values. This opt-in path cannot
#'   be combined with other random/model/scaling or systematic controls.
#' @param gambling Logical; TRUE (default) preserves the legacy randomized search.
#'   FALSE selects the systematic pooled ML workflow described below.
#' @param candidateStudies Original study IDs explicitly nominated for exclusion
#'   checks; these are not positions in the study list. Set to "auto" to
#'   nominate up to ten studies from retained study-local weak-direction or
#'   verified primary-fit diagnostics. No exclusions occur unless requested.
#' @param exclusion Character vector: "none" (default), "single", "blocks", or
#'   "all". The last excludes all nominated candidates together. May combine
#'   single, blocks and all; at least two studies must remain in every case.
#' @param studyBlocks List of vectors of original study IDs for block exclusions.
#' @param referenceStudy "current", "largestN", or one original study ID.
#'   If an exclusion removes the reference, use the largest remaining study.
#' @param referenceCandidates Optional additional original study IDs to test as
#'   reference on the full data. Set to "auto" to nominate up to three largest
#'   unflagged alternative studies (distinct observation times break size ties).
#'   Equivalent duplicate cases run only once.
#' @param quickCheck Logical; defaults to TRUE when gambling=FALSE. Run bounded
#'   pilots, not final estimates or validated confidence intervals.
#' @param quickControl Named list: maxiter=100, maxSeconds=60,
#'   fullMaxSeconds=600, maxFullChecks=0, jitter=0.01. Time caps include worker
#'   startup/setup. After all pilots, optionally confirm up to maxFullChecks
#'   cases in planned order. jitter applies only to starts after the first.
#' @param timeBudget Total wall-time budget in seconds for a systematic call.
#' @param maxCases Maximum number of distinct planned cases (default 20).
#' @param seed Nonnegative integer controlling independent systematic starts
#'   or the opt-in randomDeterminedDiffusions search. The latter restores the
#'   caller's RNG state and records a separate seed for each refit.
#' @param checkpointDirectory Directory for manifests, logs, fits and diagnostic
#'   CSV files. Defaults to a unique subdirectory of activeDirectory.
#' @param resume Reuse completed checkpoints only when data, plan, versions and
#'   settings match. Supply the original checkpointDirectory.
#' @param planOnly Logical; FALSE by default. TRUE returns the nomination
#'   report and casePlan without workers, fitting, likelihood evaluation or
#'   checkpoint-file creation. Only available with gambling=FALSE.
#' @param x,object A systematic optimization result for print or summary.
#' @param ... Additional arguments for the print method. In ctmaOptimizeFit(), unknown arguments are rejected with a migration hint where applicable.
#' @details With gambling=FALSE, the first case retains all data and the current
#'   reference. Additional reference and explicitly requested exclusion cases
#'   preserve the original observations, time scale and model constraints.
#'   This initial implementation supports pooled ML (optimize=TRUE, priors=FALSE)
#'   with dummy study coding and backend stan or ctsemjulia; it retains the
#'   baseline backend. Weighted coding, hierarchical study random effects,
#'   reduced-rank covariance, and reference changes of partially constrained
#'   study-effect families require an explicit model-specific transformation.
#'   Bundled julia searches remain available through gambling=TRUE.
#'
#'   Systematic checks require optional packages callr and ps and an installed
#'   CoTiMA2 package. Each numerical case runs in an isolated R subprocess with
#'   a time cap. Published Stan optimizers without an iteration-limit control
#'   retain the wall-time cap and record a warning. Timeout and finite
#'   short-pilot results are inconclusive. Older Stan APIs may need their
#'   standard uncertainty postprocessing to return complete TI-predictor fits;
#'   these pilots retain only 20 approximate draws and remain inconclusive.
#'   No probability of convergence is inferred from a quick check. Warnings,
#'   optimizer convergence, gradients and available curvature diagnostics are
#'   retained separately. Auto exclusion nominations use reported weak directions
#'   containing at most four mapped study-TI coordinates and no shared/unmapped
#'   coefficients, or explicit primary-fit nonconvergence/weak-direction flags
#'   when original study IDs verify alignment. Supply ctmaInitFit if separate
#'   primary fits are available. Small samples or few occasions alone do not
#'   nominate exclusions. Empty automatic nominations skip single/all exclusions;
#'   block exclusions still require explicit studyBlocks. At most nstudies-2
#'   studies are nominated, and maxCases still limits the complete plan.
#'   These descriptive flags are reference-dependent and do not prove which
#'   study causes nonidentification. The diagnostic curvature_assessed flag records whether a finite saved Hessian is available; weak_directions is NA when unassessed.
#'   A largest-N reference can improve numerical scaling
#'   but cannot establish parameter identification.
#' @return The legacy return is unchanged with gambling=TRUE. Systematic results
#'   contain bestFit (best full-data attempt, including the original baseline),
#'   bestConvergedFit, caseDiagnostics, warnings, errors, sensitivity$cases and
#'   checkpointDirectory. With auto nominations, automaticSelection contains
#'   studyDiagnostics, nominated IDs and explanatory notes; the diagnostic
#'   table is also saved as automaticStudyDiagnostics.csv. Excluded-data likelihoods never compete with full-data
#'   likelihoods. Inspect subset fits as sensitivity analyses, not replacements
#'   selected for a preferred effect. summary() returns the diagnostic table, or casePlan in planOnly mode.
#'
#' @importFrom foreach %dopar%
#' @importFrom stats runif
#' @importFrom methods is
#' @importFrom ctsem ctStanFit
#'
#' @examples
#' \dontrun{
#' optimFit313 <- ctmaOptimizeFit(primaryStudies=CoTiMAstudyList_3,
#'                                 activeDirectory="/Users/tmp/",  # adapt!
#'                                 problemStudy=which(CoTiMAstudyList_3$studyNumbers == 313),
#'                                 reFits=10,
#'                                 n.latent=2)
#' summary(optimFit313)
#' }
#'
#' @export ctmaOptimizeFit
#'
#' @return Returns a list with \code{bestFit}, the successful fit with the
#'   smallest finite \code{minus2ll}; \code{all_minus2ll}, containing one value
#'   per attempt and \code{NA} for failed or non-finite attempts; the settings
#'   used across attempts; and captured \code{warnings} and \code{errors}.
#'   If no attempt succeeds, the function stops with an error summarizing the
#'   captured fitting errors.
#'
ctmaOptimizeFit <- function(
                            activeDirectory=NULL,
                            coresToUse=c(2),
                            CoTiMAStanctArgs=NULL,
                            ctmaFitFit=NULL,
                            ctmaInitFit=NULL,
                            ctStanFit=NULL,
                            ctStanData=NULL,
                            customPar=FALSE,
                            finishsamples=NULL,
                            iter=5000,
                            primaryStudies=NULL,
                            problemStudy=NULL,
                            randomIV=FALSE,
                            randomPar=FALSE,
                            randomRI=FALSE,
                            randomScaleTI=FALSE,
                            randomScaleTime=c(1,1),
                            saveModelFits=FALSE,
                            shuffleStudyList=FALSE,
                            reFits=NULL,
                            scaleTime=NULL,
                            scaleTI=NULL,
                            verbose=1,
                            ctStanFitArgs=list(stanmodeltext = NA,
                                               iter = 1000,
                                               intoverstates = TRUE,
                                               binomial = FALSE,
                                               fit = TRUE,
                                               intoverpop = "auto",
                                               sameInitialTimes = FALSE,
                                               stationary = FALSE,
                                               plot = FALSE,
                                               derrind = NA,
                                               optimize = TRUE,
                                               optimcontrol = list(),
                                               nlcontrol = list(),
                                               nopriors = NA,
                                               priors = FALSE,
                                               chains = 2,
                                               #cores = ifelse(optimize, getOption("mc.cores", 2L), "maxneeded"),
                                               inits = NULL,
                                               compileArgs = list(),
                                               forcerecompile = FALSE,
                                               saveCompile = TRUE,
                                               savescores = FALSE,
                                               savesubjectmatrices = FALSE,
                                               saveComplexPars = FALSE,
                                               gendata = FALSE,
                                               control = list(),
                                               #verbose = 0,
                                               vb = FALSE),
                            backend=c("stan", "julia"),
                            backendArgs=list(),
                            gambling=TRUE,
                            candidateStudies=NULL,
                            exclusion="none",
                            studyBlocks=NULL,
                            referenceStudy="current",
                            referenceCandidates=NULL,
                            quickCheck=!gambling,
                            quickControl=list(),
                            timeBudget=3600,
                            maxCases=20L,
                            seed=1234L,
                            checkpointDirectory=NULL,
                            resume=FALSE,
                            planOnly=FALSE,
                            randomDeterminedDiffusions=FALSE, ...)
{
  .ctmaCheckRetiredArguments(match.call(expand.dots = FALSE), "ctmaOptimizeFit")
  if (!is.logical(gambling) || length(gambling)!=1L || is.na(gambling))
    stop("gambling must be TRUE or FALSE.", call.=FALSE)
  randomDeterminedDiffusions <- .ctmaLogicalFlag(randomDeterminedDiffusions, "randomDeterminedDiffusions")
  if (randomDeterminedDiffusions) {
    if (!gambling) stop("randomDeterminedDiffusions requires gambling=TRUE.",call.=FALSE)
    if (!is.null(ctStanFit) || is.null(ctmaFitFit))
      stop("randomDeterminedDiffusions currently requires a pooled ctmaFitFit and its ctmaInitFit; primary-study and raw ctsem searches are not supported.",call.=FALSE)
    if (isTRUE(randomIV) || isTRUE(randomRI) || isTRUE(randomPar) || isTRUE(randomScaleTI) ||
        isTRUE(shuffleStudyList) || !identical(as.numeric(randomScaleTime),c(1,1)) ||
        !is.null(scaleTime) || !is.null(scaleTI) || isTRUE(customPar))
      stop("randomDeterminedDiffusions changes only covariance coordinates; random model/scaling/shuffle options and explicit scale overrides must be disabled.",call.=FALSE)
    if (!is.null(candidateStudies) || !identical(exclusion,"none") || !is.null(studyBlocks) ||
        !identical(referenceStudy,"current") || !is.null(referenceCandidates) || isTRUE(quickCheck) ||
        resume || !missing(quickControl) || !missing(timeBudget) || !missing(maxCases) ||
        !is.null(checkpointDirectory) || isTRUE(planOnly))
      stop("Study/reference/quick-check controls require gambling=FALSE.",call.=FALSE)
    return(.ctmaOptimizeRandomDetermined(ctmaFitFit,ctmaInitFit,reFits,seed,
      backend=if(missing(backend))NULL else backend[1L],backendArgs=backendArgs,
      legacyArgs=CoTiMAStanctArgs,cores=coresToUse,iter=iter,finishsamples=finishsamples,
      activeDirectory=activeDirectory,saveModelFits=saveModelFits,verbose=verbose))
  }
  if (!gambling) {
    .ctmaVersionNote()
    .ctma_timer <- .ctmaRegisterExecutionTimer("ctmaOptimizeFit")
    on.exit(.ctma_timer(), add=TRUE)
    if (!is.null(ctStanFit) || !is.null(primaryStudies) || !is.null(problemStudy))
      stop("gambling=FALSE requires a pooled ctmaFitFit; primary-study and raw ctStanFit searches retain the legacy path.", call.=FALSE)
    if (isTRUE(randomIV) || isTRUE(randomRI) || isTRUE(randomPar) ||
        isTRUE(randomScaleTI) || isTRUE(shuffleStudyList) ||
        !identical(as.numeric(randomScaleTime),c(1,1)))
      stop("Random model/scaling/shuffle options conflict with gambling=FALSE; use explicit referenceStudy and a fixed model.",call.=FALSE)
    if (!is.null(scaleTime) || !is.null(scaleTI))
      stop("gambling=FALSE preserves the saved time/TI scale; fit a new baseline before changing it.",call.=FALSE)
    if (!is.null(CoTiMAStanctArgs))
      stop("Use backendArgs for gambling=FALSE.",call.=FALSE)
    return(.ctmaOptimizeSystematic(ctmaFitFit,ctmaInitFit,
      backend=if(missing(backend))NULL else backend[1L],backendArgs=backendArgs,
      candidateStudies=candidateStudies,exclusion=exclusion,studyBlocks=studyBlocks,
      referenceStudy=referenceStudy,referenceCandidates=referenceCandidates,
      quickCheck=quickCheck,quickControl=quickControl,timeBudget=timeBudget,
      maxCases=maxCases,seed=seed,checkpointDirectory=checkpointDirectory,
      resume=resume,coresToUse=coresToUse,reFits=if(is.null(reFits))1L else reFits,
      iter=iter,verbose=verbose,activeDirectory=activeDirectory,
      finishsamples=if(is.null(finishsamples))1000L else finishsamples,planOnly=planOnly))
  }
  if (!is.null(candidateStudies) || !identical(exclusion,"none") ||
      !is.null(studyBlocks) || !identical(referenceStudy,"current") ||
      !is.null(referenceCandidates) || isTRUE(quickCheck) || resume ||
      !missing(quickControl) || !missing(timeBudget) || !missing(maxCases) ||
      !missing(seed) || !is.null(checkpointDirectory) || isTRUE(planOnly))
    stop("Study/reference/quick-check controls require gambling=FALSE.",call.=FALSE)
  .ctmaVersionNote()
  .ctma_timer <- .ctmaRegisterExecutionTimer("ctmaOptimizeFit")
  on.exit(.ctma_timer(), add = TRUE)

  .ctmaCheckCtsemJuliaRequest(backend, "ctmaOptimizeFit")

  julia_backend_overrides <- backendArgs
  backend <- match.arg(tolower(backend), c("stan", "julia"))
  if (identical(backend, "julia")) {
    if (is.null(ctmaFitFit) || !is.null(ctStanFit)) {
      stop(
        "backend = 'julia' currently requires the complete-model ",
        "'ctmaFitFit' branch of ctmaOptimizeFit().",
        call. = FALSE
      )
    }
    if (!isTRUE(ctmaFitFit$argumentList$optimize)) {
      stop(
        "ctmaOptimizeFit() currently supports backend = 'julia' only for ML estimation.",
        call. = FALSE
      )
    }
    resolved_optimize_backend <- .ctmaResolveFitBackend(
      backend = "julia",
      backendArgs = backendArgs,
      legacy_stan_args = CoTiMAStanctArgs,
      estimation = "ML"
    )
    backendArgs <- resolved_optimize_backend$backend_args
    CoTiMAStanctArgs <- NULL
  } else if (length(backendArgs)) {
    if (!is.null(CoTiMAStanctArgs)) {
      stop(
        "Supply either 'CoTiMAStanctArgs' or 'backendArgs', not both. ",
        "Use backendArgs for new code.",
        call. = FALSE
      )
    }
    .ctmaValidateBackendArgs(
      ctmaBackendDefaults(backend = "stan"),
      backendArgs,
      path = "backendArgs"
    )
    CoTiMAStanctArgs <- backendArgs
  }

  #######################################################################################################################
  ################################################# Check Cores To Use ##################################################
  #######################################################################################################################
  {
    {
      resource_resolution <- .ctmaResolveCores(
        coresToUse, backend = backend, workflow = "ctmaOptimizeFit"
      )
      coresToUse <- resource_resolution$resolved
      if (identical(backend, "julia")) {
        julia_resources <- .ctmaInitJuliaResources(
          cores = coresToUse,
          estimation = "ML",
          backendArgs = julia_backend_overrides,
          resolved_backend_args = backendArgs,
          workflow = "ctmaOptimizeFit"
        )
        backendArgs <- julia_resources$backend_args
      }
    }

    if (!(is.null(ctStanFit))) {
      if (!(is(ctStanFit, "ctStanFit"))) {
        ErrorMsg <- "\nThe ctStanFit object provided was not created with ctStanFit! \nGood luck for the next try!"
        stop(ErrorMsg)
      }
      if (is.null(ctStanData)) {
        ErrorMsg <- "\nA ctStanFit object was provided. A data frame for the argument ctStanData is required, too! \nGood luck for the next try!"
        stop(ErrorMsg)
      }
    }

    # Resolve stored and user-supplied Stan controls by name. User overrides
    # take precedence over settings inherited from the supplied CoTiMA fit.
    if (is.null(ctStanFit) && identical(backend, "stan")) {
      inherited_stan_args <- NULL
      if (!is.null(ctmaFitFit)) {
        inherited_stan_args <- ctmaFitFit$argumentList$CoTiMAStanctArgs
      } else if (!is.null(ctmaInitFit)) {
        inherited_stan_args <- ctmaInitFit$argumentList$CoTiMAStanctArgs
      }
      CoTiMAStanctArgs <- .ctmaResolveOptimizeStanArgs(
        overrides = CoTiMAStanctArgs,
        inherited = inherited_stan_args,
        finishsamples = finishsamples
      )
    }

  }

  ########################################################################################################################

  #'%dopar%' <- foreach::'%dopar%' deprecated

  if (!(is.null(scaleTime))) {
    randomScaleTime[1] <- randomScaleTime[2] <- scaleTime
    Msg <- paste0("You provided scaleTime. This overrides randomScaleTime; ",
      "both randomScaleTime limits are set to ", scaleTime, ".\n")
    message(Msg)
  }

  if (!is.numeric(randomScaleTime) || length(randomScaleTime) != 2L ||
      anyNA(randomScaleTime) || any(!is.finite(randomScaleTime)) ||
      any(randomScaleTime <= 0)) {
    stop("randomScaleTime must contain two finite positive numbers.",
      call. = FALSE)
  }

  if (randomScaleTime[2] < randomScaleTime[1]) {

    ErrorMsg <- "\nrandomScaleTime[1] has to be <= randomScaleTime[2]! \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  if( (!(is.null(ctmaFitFit))) & (!(is.null(primaryStudies))) ) {
    ErrorMsg <- "Arguments for both ctmaFitFit and primaryStudies were provided. Only one out of the two can be chosen!"
    stop(ErrorMsg)
  }

  if( (!(is.null(ctmaFitFit))) & ((is.null(ctmaInitFit))) ) {
    ErrorMsg <- "Argument for ctmaFitFit was provided but not for ctmaInitFit. Need a ctmaInitFit, too!"
    stop(ErrorMsg)
  }

  if( (!(is.null(ctmaFitFit))) & (!(is.null(ctmaInitFit))) ) {
    if (ctmaFitFit$argumentList$ctmaInitFit != deparse(substitute(ctmaInitFit)))  {
      ErrorMsg <- paste0("The wrong ctmaInitFit object was provided. I need ",  ctmaFitFit$argumentList$ctmaInitFit, "!")
      stop(ErrorMsg)
    }
  }

  if((!is.null(ctStanFit)) & ( (!(is.null(ctmaFitFit))) | ((!is.null(ctmaInitFit))) ) ) {
    ErrorMsg <- "A ctStanFit was provided together with a ctmaFitFit or ctmaIniFit object. Make a decision!"
    stop(ErrorMsg)
  }

  if ( is.null(ctStanFit) & is.null(ctmaFitFit) & is.null(ctmaInitFit) )  {
    ErrorMsg <- "The ctmaFitFit or ctmaIniFit (or ctStanFit) argument has to be used."
    stop(ErrorMsg)
  }

  if (length(reFits) != 1L || is.na(reFits) || reFits < 1L ||
      reFits != as.integer(reFits)) {
    stop("reFits must be a positive integer.", call. = FALSE)
  }
  reFits <- as.integer(reFits)

  if (!is.logical(randomScaleTI) || length(randomScaleTI) != 1L ||
      is.na(randomScaleTI)) {
    stop("randomScaleTI must be TRUE or FALSE.", call. = FALSE)
  }


  # Moderator Checks Moved to Sectioon where ctmaInit is optimized (not relevant if ctmaFit is optimized) # CHD Auf 2023


  # INIT Fit
  if ((is.null(ctmaFitFit)) & (is.null(ctStanFit))) {
    # CHD changed 21 SEP 2022
    ErrorMsg <- "argument primaryStudies is missing"
    if (is.null(primaryStudies))  stop(ErrorMsg)
    ErrorMsg <- "argument problemStudy is missing"
    if (is.null(problemStudy)) stop(ErrorMsg)
    ErrorMsg <- "argument activeDirectory is missing"
    if (is.null(activeDirectory)) stop(ErrorMsg)

    # CHD 14.7.2025
    if (is.null(ctmaInitFit$argumentList$indVarying)) ctmaInitFit$argumentList$indVarying <- FALSE
    if (is.null(ctmaInitFit$argumentList$randomIntercepts)) ctmaInitFit$argumentList$randomIntercepts <- FALSE
    indVarying <- ctmaInitFit$argumentList$indVarying
    randomIntercepts <- ctmaInitFit$argumentList$randomIntercepts

    # CHD 14.7.2025
    if (!is.null(randomIV) & !is.null(randomRI)) {
      if ((randomIV == FALSE) & (randomRI == FALSE)) {
        #if (is.null(ctmaInitFit$argumentList$randomIntercepts)) ctmaInitFit$argumentList$randomIntercepts <- FALSE
        #if (is.null(ctmaInitFit$argumentList$indVarying)) ctmaInitFit$argumentList$indVarying <- FALSE
        Msg <- paste0("Both \"randomIV\" and \"randomRI\" are set FALSE (default). I use the same settings as used for ctmaInit(), which was
      \"randomIntercepts = ", ctmaInitFit$argumentList$randomIntercepts, "\", and
      \"indVarying = ", ctmaInitFit$argumentList$indVarying, ".\n")
        indVarying <- ctmaInitFit$argumentList$indVarying
        randomIntercepts <- ctmaInitFit$argumentList$randomIntercepts
        message(Msg)
      }
    }


    # CHD 14.7.2025
    if (!is.null(randomIV) & !is.null(randomRI)) {
      if ((randomIV == TRUE) & (randomRI == TRUE)) {
        Msg <- "Both \"randomIV\" and \"randomRI\" are set TRUE. Both are statistically equivalent, but since the latter is more difficult to fit I set \"randomRI = FALSE\". DO \"randomIV = FALSE\" to prevent this.\n"
        message(Msg)
        randomRI <- FALSE
      }
    }

    # CHD 14.7.2025
    if (is.null(randomIV)) {
      Msg <- "\"randomIV\" was set to NULL, which overrides all argument set previously in ctmaInit() and sets  \"indVarying = FALSE\". Take care.\n"
      randomIV <- FALSE
      indVarying <- FALSE
      message(Msg)
    }

    # CHD 14.7.2025
    if (is.null(randomRI)) {
      Msg <- "\"randomRI\" was set to NULL, which overrides all argument set previously in ctmaInit() and sets  \"randomIntercepts = FALSE\". Take care.\n"
      randomRI <- FALSE
      randomIntercepts <- FALSE
      message(Msg)
    }


    # create new study list with a single problem study only
    listElements <- names(primaryStudies); listElements
    newStudyList <- as.list(listElements)
    validElements <- c("deltas", "sampleSizes", "pairwiseNs", "empcovs", "moderators", "startValues", "studyNumbers", "rawData", "empMeans", "empVars",
                       "source", "ageM", "malePercent", "occupation", "country", "alphas", "targetVariables", "recodeVariables", "combineVariables",
                       "combineVariablesNames", "missingVariables", "inits", "emprawList") #, "n.studies", "summary", "excelSheets", "plot.type")
    counter <- 0
    for (i in listElements) {
      counter <- counter + 1
      if (i %in% validElements) {
        if (i %in% c("pairwiseNs", "empcovs", "rawData", "deltas", "emprawList")) {
          newStudyList[[counter]] <- primaryStudies[[counter]][problemStudy]
        } else {
          newStudyList[[counter]] <- list(unlist(primaryStudies[[counter]][problemStudy], recursive=TRUE))
        }
      } else {
        newStudyList[[counter]] <- unlist(primaryStudies[[counter]])
      }
      if (is.logical(newStudyList[[counter]])) newStudyList[[counter]] <- NA
    }
    names(newStudyList) <- names(primaryStudies)
    newStudyList$n.studies <- 1

    currentLL <- 10^20; currentLL
    all_minus2ll <- c()
    all_scaleTime <- all_customPar <- c()
    all_scaleTI <- logical()
    all_usedStudyList <- primaryStudies
    all_randomIV <- all_randomRI <- c()
    warns <- errs <- list()


    for (i in 1:reFits) {
      #i <- 1
      scaleTime <- round(stats::runif(1, min=randomScaleTime[1], max=randomScaleTime[2]), 2)
      if (randomPar == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        customPar = c(TRUE, FALSE)[tmp1]
      }
      if (!(is.null(randomScaleTime))) {
        Msg <- paste0("Argument scaleTime is set to: ", scaleTime, ".")
        message(Msg)
      }
      if (randomPar == TRUE) {
        Msg <- paste0("Argument customPar is set to: ", customPar, ".")
        message(Msg)
      }
      all_scaleTime <- c(all_scaleTime, scaleTime)
      all_customPar <- c(all_customPar, customPar)

      if (randomScaleTI == TRUE) {
        scaleTI <- sample(c(TRUE, FALSE), size = 1L)
        message("Argument scaleTI is set to: ", scaleTI, ".")
      } else if (is.null(scaleTI)) {
        scaleTI <- ctmaInitFit$argumentList$scaleTI
        if (is.null(scaleTI)) scaleTI <- CoTiMAStanctArgs$scaleTI
      }
      all_scaleTI <- c(all_scaleTI, scaleTI)

      if (randomRI == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        randomIntercepts <- c("MANIFEST", "CINT")[tmp1]
        all_randomRI <- c(all_randomRI, randomIntercepts)
      } else {
        if (randomIntercepts != FALSE) randomIntercepts <- ctmaInitFit$argumentList$randomIntercepts
        all_randomRI <- c(all_randomRI, randomIntercepts)
      }

      if (randomIV == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        indVarying <- c("MANIFEST", "CINT")[tmp1]
        all_randomIV <- c(all_randomIV, indVarying)
      } else {
        if (indVarying != FALSE) indVarying <- ctmaInitFit$argumentList$indVarying
        all_randomIV <- c(all_randomIV, indVarying)
      }
      #

      if (is.null(finishsamples)) finishsamples <- ctmaInitFit$argumentList$finishsamples

      # CHD 12.4.24
      #if (is.null(indVarying)) indVarying <- ctmaFitFit$argumentList$indVarying

      problem <- FALSE
      #print(randomIntercepts)
      #print(indVarying)
      fit <- tryCatch(ctmaInit(primaryStudies=newStudyList,
                               coresToUse = coresToUse, # changed Aug 2023
                               scaleTime = scaleTime,
                               scaleTI=scaleTI,
                               customPar=customPar,
                               finishsamples=finishsamples,
                               iter=iter,
                               activeDirectory = activeDirectory,
                               CoTiMAStanctArgs=CoTiMAStanctArgs,
                               n.latent=ctmaInitFit$argumentList$n.latent,
                               n.manifest=ctmaInitFit$argumentList$n.manifest,
                               #indVarying = ctmaInitFit$argumentList$indVarying,
                               indVarying = indVarying,
                               checkSingleStudyResults=FALSE,
                               T0means=ctmaInitFit$argumentList$T0means,
                               manifestMeans=ctmaInitFit$argumentList$manifestMeans,
                               manifestVars=ctmaInitFit$argumentList$manifestVars,
                               chains=ctmaInitFit$argumentList$chains,
                               cint=ctmaInitFit$argumentList$cint,
                               diff=ctmaInitFit$argumentList$diff,
                               digits=ctmaInitFit$argumentList$digits,
                               drift=ctmaInitFit$argumentList$drift,

                               indVaryingT0=ctmaInitFit$argumentList$indVaryingT0,
                               lambda=ctmaInitFit$argumentList$lambda,
                               #loadSingleStudyModelFit=loadSingleStudyModelFit,
                               #nopriors=nopriors,
                               optimize=ctmaInitFit$argumentList$optimize,
                               #primaryStudies=primaryStudies,
                               priors=ctmaInitFit$argumentList$priors,
                               sameInitialTimes=ctmaInitFit$argumentList$sameInitialTimes,
                               #saveRawData=saveRawData,
                               #saveSingleStudyModelFit=saveSingleStudyModelFit,
                               #silentOverwrite=silentOverwrite,
                               T0var=ctmaInitFit$argumentList$T0var,
                               useSV=ctmaInitFit$argumentList$useSV,
                               verbose=verbose,
                               #randomIntercepts=ctmaInitFit$argumentList$randomInterceptsSettings,
                               randomIntercepts=randomIntercepts),
                      error = function(e) {
                        problem <<- TRUE
                        errs[[i]] <<- e
                        NULL
                      }
      )


      fit_value <- .ctmaOptimizeFitValue(fit)
      if (!problem && is.finite(fit_value)) {
        all_minus2ll <- c(all_minus2ll, fit_value)

        if (saveModelFits != FALSE) {
          saveRDS(fit, paste0(activeDirectory, saveModelFits, " ", i, " .rds"))
        }

        if (fit_value < currentLL) {
          currentLL <- fit_value
          bestFit <- fit
          usedStudyList <- ctmaInitFit$primaryStudyList
          usedTimeScale <- scaleTime
          usedScaleTI <- scaleTI
        }
      } else {
        all_minus2ll <- c(all_minus2ll, NA_real_)
      }


    }
  }

  # ctmaFitFit
  if (!(is.null(ctmaFitFit))) {
    if (!(is(ctmaFitFit, "CoTiMAFit")) &&
        !inherits(ctmaFitFit, "ctmaJuliaFit")) {
      ErrorMsg <- "The ctmaFitFit object must have been created with ctmaFit()."
      stop(ErrorMsg)
    }

    # CHD 14.7.2025
    if (!is.null(randomIV) & !is.null(randomRI)) {
      if ((randomIV == FALSE) & (randomRI == FALSE)) {
        Msg <- paste0("Both \"randomIV\" and \"randomRI\" are set FALSE (default). I use the same settings as used for ctmaFit(), which was
      \"randomIntercepts = ", ctmaFitFit$argumentList$randomIntercepts, "\", and
      \"indVarying = ", ctmaFitFit$argumentList$indVarying, ".\n")
        indVarying <- ctmaFitFit$argumentList$indVarying
        randomIntercepts <- ctmaFitFit$argumentList$randomIntercepts
        message(Msg)
      }
    }

    # CHD 14.7.2025
    if (!is.null(randomIV) & !is.null(randomRI)) {
      if ((randomIV == TRUE) & (randomRI == TRUE)) {
        Msg <- "Both \"randomIV\" and \"randomRI\" are set TRUE. Both are statistically equivalent, but since the latter is more difficult to fit I set \"randomRI = FALSE\". DO \"randomIV = FALSE\" to prevent this.\n"
        message(Msg)
        randomRI <- FALSE
      }
    }

    # CHD 14.7.2025
    if (is.null(randomIV)) {
      Msg <- "\"randomIV\" was set to NULL, which overrides all argument set previously in ctmaInit() and sets  \"indVarying = FALSE\". Take care.\n"
      randomIV <- FALSE
      indVarying <- FALSE
      message(Msg)
    }

    # CHD 14.7.2025
    if (is.null(randomRI)) {
      Msg <- "\"randomRI\" was set to NULL, which overrides all argument set previously in ctmaInit() and sets  \"randomIntercepts = FALSE\". Take care.\n"
      randomRI <- FALSE
      randomIntercepts <- FALSE
      message(Msg)
    }


    currentLL <- 10^20; currentLL
    all_minus2ll <- c()
    all_scaleTime <- all_customPar <- c()
    all_scaleTI <- c()
    all_usedStudyList <- c()
    all_randomIV <- all_randomRI <- c()
    warns <- errs <- list()

    for (i in 1:reFits) {
      cat(paste0("This is fit attempt #", i, " out of ", reFits, "re-fits."))
      scaleTime <- round(stats::runif(1, min=randomScaleTime[1], max=randomScaleTime[2]), 2)
      all_scaleTime <- c(all_scaleTime, scaleTime)
      if (randomPar == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        customPar <- c(TRUE, FALSE)[tmp1]
      } else {
        if (is.null(customPar)) customPar <- ctmaFitFit$argumentList$customPar
      }
      #
      all_customPar <- c(all_customPar, customPar)
      if (randomScaleTI == TRUE) {
        scaleTI <- sample(c(TRUE, FALSE), size = 1L)
      } else {
        if (is.null(scaleTI)) scaleTI <- ctmaFitFit$argumentList$scaleTI
      }
      all_scaleTI <- c(all_scaleTI, scaleTI)
      #
      if (randomIV == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        indVarying <- c("MANIFEST", "CINT")[tmp1]
        all_randomIV <- c(all_randomIV, indVarying)
      } else {
        indVarying <- ctmaFitFit$argumentList$indVarying
        all_randomIV <- c(all_randomIV, indVarying)
      }
      #
      if (randomRI == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        randomIntercepts <- c("MANIFEST", "CINT")[tmp1]
        all_randomRI <- c(all_randomRI, randomIntercepts)
      } else {
        randomIntercepts <- ctmaFitFit$argumentList$randomIntercepts
        all_randomRI <- c(all_randomRI, randomIntercepts)
      }
      #
      if (shuffleStudyList == TRUE) {
        #
        tmpStudyList <- ctmaInitFit$studyList; length(tmpStudyList)
        studyNumbers <- unlist(lapply(tmpStudyList, function(x) x$originalStudyNo)); studyNumbers
        newStudyOrder <- sample(studyNumbers, length(studyNumbers), replace=FALSE); newStudyOrder
        ctmaInitFit <- .ctmaOptimizeShuffleInit(ctmaInitFit, newStudyOrder)
        all_usedStudyList <- c(all_usedStudyList, newStudyOrder)
      }

      if (!(is.null(randomScaleTime))) {
        Msg <- paste0("Argument scaleTime is set to: ", scaleTime, ".")
        message(Msg)
      }
      if (randomPar == TRUE) {
        Msg <- paste0("Argument customPar is set to: ", customPar, ".")
        message(Msg)
      }
      if (randomScaleTI == TRUE) {
        Msg <- paste0("Argument scaleTI is set to: ", scaleTI, ".")
        message(Msg)
      }
      if (shuffleStudyList == TRUE ) {
        tmp <- unlist(ctmaInitFit$primaryStudyList$studyNumbers)
        #tmp <- tmp[-length(tmp)]
        tmp <- paste(tmp, collapse=" ")
        Msg <- paste0("Order of studies in the shuffled study list is: ", tmp, ".")
        message(Msg)
      }

      if (is.null(finishsamples)) finishsamples <- ctmaFitFit$argumentList$finishsamples
      if (is.null(iter)) iter <- 5000

      fit <- tryCatch(withCallingHandlers(
        expr = ctmaFit(ctmaInitFit=ctmaInitFit,
                       primaryStudyList=ctmaInitFit$primaryStudyList,
                       cluster=ctmaFitFit$argumentList$cluster,
                       activeDirectory=activeDirectory,

                       digits=ctmaFitFit$argumentList$digits,
                       drift=ctmaFitFit$argumentList$drift,
                       invariantDrift=ctmaFitFit$argumentList$invariantDrift,
                       moderatedDrift=ctmaFitFit$argumentList$moderatedDrift,
                       equalDrift=ctmaFitFit$argumentList$equalDrift,
                       mod.number=ctmaFitFit$argumentList$mod.number,
                       mod.type=ctmaFitFit$argumentList$mod.type,
                       mod.names=ctmaFitFit$argumentList$mod.names,
                       #indVarying=ctmaFitFit$argumentList$indVarying,
                       indVarying=indVarying,
                       coresToUse=coresToUse, # changed Aug 2023
                       sameInitialTimes=ctmaFitFit$argumentList$sameInitialTimes,
                       scaleTI=scaleTI,
                       scaleMod=ctmaFitFit$argumentList$scaleMod,
                       transfMod=ctmaFitFit$argumentList$transfMod,
                       scaleClus=ctmaFitFit$argumentList$scaleClus,
                       #scaleTime=ctmaFitFit$argumentList$scaleTime,
                       scaleTime=scaleTime,
                       optimize=ctmaFitFit$argumentList$optimize,
                       #nopriors=ctmaFitFit$argumentList$nopriors,
                       finishsamples=finishsamples,
                       iter=iter,
                       chains=ctmaFitFit$argumentList$chains,
                       verbose=verbose,
                       allInvModel=ctmaFitFit$argumentList$allInvModel,
                       customPar=customPar,
                       inits=ctmaFitFit$argumentList$inits,
                       modsToCompare=ctmaFitFit$argumentList$modsToCompare,
                       catsToCompare=ctmaFitFit$argumentList$catsToCompare,
                       driftsToCompare=ctmaFitFit$argumentList$driftsToCompare,
                       useSampleFraction=ctmaFitFit$argumentList$useSampleFraction,
                       T0means=ctmaFitFit$argumentList$T0means,
                       manifestMeans=ctmaFitFit$argumentList$manifestMeans,
                       CoTiMAStanctArgs=CoTiMAStanctArgs,
                       #randomIntercepts=ctmaFitFit$argumentList$randomIntercepts,
                       randomIntercepts=randomIntercepts,
                       manifestVars=ctmaFitFit$argumentList$manifestVars,
                       WEC=FALSE,
                       studyCoding=if (is.null(ctmaFitFit$argumentList$studyCoding)) {
                         if (isTRUE(ctmaFitFit$argumentList$WEC)) "weighted" else "dummy"
                       } else ctmaFitFit$argumentList$studyCoding,
                       studyVarying=if (is.null(ctmaFitFit$argumentList$studyVarying)) {
                         if (isTRUE(ctmaFitFit$argumentList$WEC)) {
                           c("CINT", "DIFFUSION", "T0VAR", "DRIFT")
                         } else c("CINT", "DIFFUSION", "T0VAR")
                       } else ctmaFitFit$argumentList$studyVarying,
                       priors=ctmaFitFit$argumentList$priors,
                       binaries=ctmaFitFit$argumentList$binaries,
                       T0var=ctmaFitFit$argumentList$T0var,
                       ind.mod.names=ctmaFitFit$argumentList$ind.mod.names,
                       ind.mod.number=ctmaFitFit$argumentList$ind.mod.number,
                       ind.mod.type=ctmaFitFit$argumentList$ind.mod.type,
                       cint=ctmaFitFit$argumentList$cint,
                       indVaryingT0=ctmaFitFit$argumentList$indVaryingT0,
                       # ctmaOptimizeFit() always performs complete refits,
                       # including when ctmaFitFit is a fit=FALSE template.
                       fit=TRUE,
                       backend=backend,
                       backendArgs=if (identical(backend, "julia")) backendArgs else list()
        ),
        warning = function(w) {warns[[i]] <<- w}),
        error = function(e) {
          errs[[i]] <<- e
          NULL
        }
      )

      fit_value <- .ctmaOptimizeFitValue(fit)
      all_minus2ll <- c(all_minus2ll, fit_value)
      if (is.finite(fit_value)) {

        if (saveModelFits != FALSE) {
          saveRDS(fit, paste0(activeDirectory, saveModelFits, " ", i, " .rds"))
        }

        if (fit_value < currentLL) {
          currentLL <- fit_value
          bestFit <- fit
          usedStudyList <- ctmaInitFit$primaryStudyList
          usedTimeScale <- scaleTime
          usedScaleTI <- scaleTI
        }
      } # end if(is.list(fit))
    }
  }

  if (!(is.null(ctStanFit))) {
    currentLL <- 10^20; currentLL
    all_minus2ll <- c()
    all_scaleTime <- all_customPar <- c()
    all_scaleTI <- "Available only if CoTiMA models are optimized."
    all_usedStudyList <- "Available only if CoTiMA models are optimized."
    all_randomIV <- all_randomRI <-
      "Available only if CoTiMA models are optimized."
    warns <- errs <- list()
    for (i in 1:reFits) {
      #i <- 1
      scaleTime <- round(stats::runif(1, min=randomScaleTime[1], max=randomScaleTime[2]), 2)
      all_scaleTime <- c(all_scaleTime, scaleTime)
      if (randomPar == TRUE) {
        tmp1 <- round(stats::runif(1, min=1, max=2), 0); tmp1
        customPar <- c(TRUE, FALSE)[tmp1]
      } else {
        if (is.null(customPar)) customPar <- FALSE # ctmaFitFit$argumentList$customPar
      }
      all_customPar <- c(all_customPar, customPar)
      scaleTI <- FALSE
      #
      if (!(is.null(randomScaleTime))) {
        Msg <- paste0("Argument scaleTime is set to: ", scaleTime, ".")
        message(Msg)
      }
      if (randomPar == TRUE) {
        Msg <- paste0("Argument customPar is set to: ", customPar, ".")
        message(Msg)
      }
      #if (randomScaleTI == TRUE) {
      #  Msg <- paste0("Argument scaleTI is set to: ", scaleTI, ".")
      #  message(Msg)
      #}
      #if (shuffleStudyList == TRUE ) {
      #  tmp <- unlist(ctmaInitFit$primaryStudyList$studyNumbers)
      #  #tmp <- tmp[-length(tmp)]
      #  tmp <- paste(tmp, collapse=" ")
      #  Msg <- paste0("Order of studies in the shuffled study list is: ", tmp, ".")
      #  message(Msg)
      #}

      if (is.null(finishsamples)) finishsamples <- 1000
      if (is.null(iter)) iter <- 5000

      datalong <- cbind(ctStanFit$data$subject, matrix(ctStanFit$data$time, ncol=1),
                        ctStanFit$data$Y, ctStanFit$data$tdpreds)

      if (ctStanFit$standata$ntipred == 0) {
        colnames(datalong) <- colnames(ctStanFit$ctdatastruct)[1:ncol(datalong)]
      }

      if (ctStanFit$standata$ntipred > 0) {
        datalong <- cbind(datalong, matrix(ctStanFit$data$tipreds,
                                           nrow=length(ctStanFit$data$subject)))
        colnames(datalong) <- colnames(ctStanFit$ctdatastruct)
      }

      datalong[,2] <- datalong[,2] * scaleTime
      ctStanModel <- ctStanFit$ctstanmodelbase
      ctStanFitArgs$optimcontrol$finishsamples <- finishsamples

      fit <- tryCatch(withCallingHandlers(
        expr = ctsem::ctStanFit(datalong=datalong,
                                ctstanmodel=ctStanModel,
                                coresToUse=coresToUse,
                                finishsamples=finishsamples,
                                stanmodeltext = ctStanFitArgs$stanmodeltext,
                                iter=iter,
                                intoverstates = ctStanFitArgs$intoverstates,
                                binomial = ctStanFitArgs$binomial,
                                fit = ctStanFitArgs$fit,
                                intoverpop = ctStanFitArgs$intoverpop,
                                sameInitialTimes = ctStanFitArgs$sameInitialTimes,
                                stationary = ctStanFitArgs$stationary,
                                plot = ctStanFitArgs$plot,
                                derrind = ctStanFitArgs$derrind,
                                optimize = ctStanFitArgs$optimize,
                                optimcontrol = ctStanFitArgs$optimcontrol,
                                nlcontrol = ctStanFitArgs$nlcontrol,
                                nopriors = ctStanFitArgs$nopriors,
                                priors = ctStanFitArgs$priors,
                                chains = ctStanFitArgs$chains,
                                cores = coresToUse,
                                inits = ctStanFitArgs$inits,
                                compileArgs = ctStanFitArgs$compileArgs,
                                forcerecompile = ctStanFitArgs$forcerecompile,
                                saveCompile = ctStanFitArgs$saveCompile,
                                savescores = ctStanFitArgs$savescores,
                                savesubjectmatrices = ctStanFitArgs$savesubjectmatrices,
                                saveComplexPars = ctStanFitArgs$saveComplexPars,
                                gendata = ctStanFitArgs$gendata,
                                control = ctStanFitArgs$control,
                                verbose = verbose
        ),
        warning = function(w) {warns[[i]] <<- w}),
        error = function(e) {
          errs[[i]] <<- e
          NULL
        }
      )
      if (is.list(fit) && !inherits(fit, "condition")) {
        fit$summary <- summary(fit)
        fit$summary$minus2ll <- fit$summary$logposterior * -2
      }
      fit_value <- .ctmaOptimizeFitValue(fit)
      all_minus2ll <- c(all_minus2ll, fit_value)

      if (is.finite(fit_value)) {
        if (saveModelFits != FALSE) {
          saveRDS(fit, paste0(activeDirectory, saveModelFits, " ", i, " .rds"))
        }

        if (fit_value < currentLL) {
          currentLL <- fit_value
          bestFit <- fit
          usedStudyList <- "Available only if CoTiMA models are optimized."
          usedTimeScale <- scaleTime
          usedScaleTI <- "Available only if CoTiMA models are optimized."
        }
      }
    }
    resultsSummary <- "Available only if CoTiMA models are optimized."
  } else if (exists("bestFit", inherits = FALSE)) {
    resultsSummary <- if (inherits(bestFit, "ctmaJuliaFit")) {
      "Julia backend fit; use summary(bestFit) for transformed results."
    } else {
      bestFit$studyFitList[[1]]$resultsSummary
    }
  }

  if (!exists("bestFit", inherits = FALSE)) {
    error_messages <- vapply(errs, function(e) {
      if (inherits(e, "condition")) conditionMessage(e) else as.character(e)
    }, character(1L))
    detail <- if (length(error_messages)) paste0(
      " Errors: ", paste(unique(error_messages), collapse = " | ")
    ) else ""
    stop(
      "All ", reFits, " optimization attempts failed or returned a ",
      "non-finite fit criterion.", detail,
      call. = FALSE
    )
  }

  if (is.null(ctStanFit)) {
    fit_for_prep <- if (inherits(bestFit, "ctmaJuliaFit")) {
      bestFit$source_object
    } else {
      bestFit
    }
    result_used_studies <- .ctmaOptimizePrepResult(fit_for_prep)
    result_time_scale <- bestFit$argumentList$scaleTime
    result_scale_ti <- bestFit$argumentList$scaleTI
  } else {
    result_used_studies <- usedStudyList
    result_time_scale <- usedTimeScale
    result_scale_ti <- usedScaleTI
  }

  results <- list(bestFit=bestFit, all_minus2ll=all_minus2ll, summary=bestFit$summary,
                  usedStudyList=result_used_studies,
                  usedTimeScale=result_time_scale,
                  usedScaleTI=result_scale_ti,
                  randomIV = randomIV,
                  randomRI = randomRI,
                  resultsSummary=resultsSummary,
                  all_scaleTime = all_scaleTime,
                  all_customPar = all_customPar,
                  all_scaleTI = all_scaleTI,
                  all_randomIV = all_randomIV,
                  all_randomRI = all_randomRI,
                  all_usedStudyList = all_usedStudyList,
                  warnings = warns,
                  errors = errs,
                  backend = backend,
                  backendArgs = if (identical(backend, "julia")) backendArgs else CoTiMAStanctArgs

  )
  class(results) <- "CoTiMAFit"

  invisible(results)
}


# Reorder by original study identity, not by unrelated prep-cache lengths.
# Historical prep objects may retain a full-data startValues cache after a
# study subset; the per-study records are authoritative in that case.
.ctmaOptimizeShuffleInit <- function(init, order) {
  studies <- init$studyList
  ids <- vapply(studies,function(z) as.character(z$originalStudyNo),character(1))
  order <- as.character(order)
  if (!length(ids) || anyNA(ids) || anyDuplicated(ids) ||
      length(order)!=length(ids) || anyDuplicated(order) || !setequal(ids,order))
    stop("Study shuffle requires a permutation of unique original study IDs.",call.=FALSE)
  ix <- match(order,ids); n <- length(ids)
  init$studyList <- studies[ix]
  for(j in seq_len(n))init$studyList[[j]]$studyNumber <- studies[[ix[j]]]$originalStudyNo
  map <- c(studyNumbers="originalStudyNo",empcovs="empcov",deltas="delta_t",
    sampleSizes="sampleSize",moderators="moderators",pairwiseNs="pairwiseN",
    rawData="rawData",source="source",startValues="startValues")
  for(name in intersect(names(map),names(init$primaryStudyList))) {
    old <- init$primaryStudyList[[name]]
    if(length(old)==n) init$primaryStudyList[[name]] <- old[ix]
    else if(length(old)) {
      key <- unname(map[name])
      values <- lapply(init$studyList,function(z) z[[key]])
      if(any(vapply(values,is.null,logical(1))))stop(
        "Cannot align prep field '",name,"' to original studies; regenerate the initialization object.",call.=FALSE)
      init$primaryStudyList[[name]] <- values
    }
  }
  for(name in c("emprawList","studyFitList")) {
    old <- init[[name]]
    if(length(old)==n)init[[name]]<-old[ix]
    else if(length(old))stop("Cannot align ",name," with original study IDs.",call.=FALSE)
  }
  if(length(init$primaryStudyList$emprawList))init$primaryStudyList$emprawList<-init$emprawList
  init
}


.ctmaOptimizePrepResult <- function(fit) {
  cached <- !is.null(fit$emprawList)
  result <- ctmaFitToPrep(fit,reUseEmprawData=cached)
  result$reUseEmprawData <- cached
  if(!cached).ctma_print_message("No study-wise wide-data cache was retained in the best pooled fit. usedStudyList contains study definitions; cached-data reuse is disabled. The best fit retains its long data separately.")
  result
}
