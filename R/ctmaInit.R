#' ctmaInit
#'
#' @description Fits ctsem models to each primary study in the supplied list of primary studies prepared by \code{\link{ctmaPrep}}.
#'
#' @param model Optional reusable model from \code{ctmaFitToModel()}. Explicit arguments override inherited settings.
#' @param activeDirectory defines another active directory than the one used in \code{\link{ctmaPrep}}
#' @param binaries which manifest is a binary. Still experimental
#' @param chains number of chains to sample, during HMC or post-optimization importance sampling.
#' @param checkSingleStudyResults Displays estimates from single study ctsem models and waits for user input to continue. Useful to check estimates before they are saved.
#' @param cint default 'auto' (= 0). Are set free if random intercepts model with varying cints is requested (by indVarying='cint')
#' @param coresToUse if negative, the value is added to the number of available
#'   cores; thus \code{-1} uses all but one core. With the Julia backend this is
#'   automatically translated into an ML or Bayesian thread allocation unless
#'   an explicit Julia thread control is supplied in \code{backendArgs}.
#'   Explicit controls must remain within the resolved CPU budget.
#' @param CoTiMAStanctArgs parameters that can be set to improve model fitting of the \code{\link[ctsem]{ctFit}} Function
#' @param customPar logical. If set TRUE leverages the first pass using priors and ensure that the drift diagonal cannot easily go too negative (helps since ctsem > 3.4)
#' @param diff Diffusion-matrix specification with exactly
#'   \code{n.latent^2} entries. Character entries are freely estimated labels
#'   (identical labels impose equality constraints); numeric entries are fixed
#'   values, for example zero to omit a covariance. Entries follow R's
#'   column-major matrix order; the diffusion matrix is symmetric and only its
#'   lower triangle is parameterized.
#' @param digits number of digits used for rounding (in outputs)
#' @param doPar parallel and multiple fitting if single studies. A value > 1 will fit each study doPar times in parallel mode during which no output is generated (screen remains silent). Useful to obtain best fit.
#' @param drift Drift-matrix specification with exactly
#'   \code{n.latent^2} entries. Character entries are freely estimated labels
#'   (identical labels impose equality constraints); numeric entries are fixed
#'   values, such as zero to exclude an effect.
#' @param finishsamples number of samples to draw (either from Hessian-based
#'   covariance or the posterior distribution) for final results computation
#'   (default = 1000). With the Julia backend this sets
#'   \code{backendArgs$transformed_draws}, unless that control is supplied
#'   explicitly in \code{backendArgs}.
#' @param fit TRUE (default) fits the requested model. FALSE returns the generated
#'   long-format data, the editable ctsem model in \code{ctModel}, and a
#'   copyable, executable R specification of that model in \code{ctModelText}.
#'   The model can then be modified and fitted using
#'   \code{\link[ctsem]{ctFit}}.
#' @param indVarying control for unobserved heterogeneity through randomly
#'   (inter-individually) varying manifest means. With
#'   \code{backend="ctsemjulia"}, \code{"MANIFEST"} or \code{"CINT"}
#'   activates ctsem's native Laplace integration. This is distinct from
#'   \code{randomIntercepts}, which adds explicit latent intercept processes.
#' @param indVaryingT0 NULL (default) or TRUE retains the automatic initial-mean
#'   random-effect specification selected by indVarying. FALSE disables individual
#'   T0MEANS random effects while retaining manifest/CINT random effects and
#'   the within-person T0VAR specification. In particular, the supported
#'   stationary determined-diffusion manifest-RI model requires FALSE.
#' @param iter number of interation (defaul = 1000). Sometimes larger values could be required fom Bayesian estimation
#' @param lambda R-type matrix with pattern of fixed (=1) or free (any string) loadings.
#' @param loadSingleStudyModelFit load the fit of single study ctsem models
#' @param manifestMeans Scalar 0 fixes all manifest means; "auto" assigns free
#'   labels. A vector mixes labels and fixed values, e.g. c("mmV1", 0).
#'   Named lists for manifestMeans, cint and T0means must name every manifest
#'   or latent process, respectively. Numbers fix values and labels free them.
#'   Explicit fixed entries remain fixed with indVarying and randomIntercepts.
#'   A binary manifest offset differs from a dynamic CINT and can represent
#'   nonzero baseline log odds with centered initial latent means.
#' @param manifestVars define the error variances of the manifests within a single time point using R-type lower triangular matrix with nrow=n.manifest & ncol=n.manifest.
#' @param n.latent number of latent variables of the model (hast to be specified)!
#' @param n.manifest number of manifest variables of the model (if left empty it will assumed to be identical with n.latent).
#' @param optimize if set to FALSE, Stan's Hamiltonian Monte Carlo sampler is used (default = TRUE = maximum a posteriori / importance sampling) .
#' @param primaryStudies list of primary study information created with \code{\link{ctmaPrep}}
#' @param priors With optimization, TRUE requests MAP for the stan and ctsemjulia backends.
#'   The bundled julia optimizer supports ML only and rejects TRUE when optimize=TRUE.
#'   With sampling, Bayesian priors are handled by the selected sampler. If FALSE, any priors are disabled – sometimes desirable for optimization
#' @param randomIntercepts (default = FALSE) Experimental. Adds explicit
#'   latent intercept processes and overrides ctsem's default representation
#'   of individually varying intercepts. It is not synonymous with
#'   \code{indVarying}; use \code{indVarying="MANIFEST"} or \code{"CINT"}
#'   to request person-level random effects integrated by Laplace under
#'   \code{backend="ctsemjulia"}.
#' @param sameInitialTimes Only important for raw data. If TRUE (default=FALSE), T0MEANS occurs for every subject at the same time, rather than just at the earliest observation.
#' @param saveRawData save (synthesised pseudo) raw date. List: saveRawData$studyNumbers, $fileName, $row.names, col.names, $sep, $dec
#' @param saveSingleStudyModelFit save the fit of single study ctsem models (could save a lot of time afterwards if the fit is loaded)
#' @param scaleTI scale TI predictors
#' @param scaleTime scale time (interval) - sometimes desirable to improve fitting
#' @param silentOverwrite overwrite old files without asking
#' @param T0means Default 0 (assuming standardized variables). Can be assigned labels to estimate them freely.
#' @param T0var (default = 'auto')
#' @param useSV if TRUE (default=FALSE) start values will be used if provided in the list of primary studies
#' @param verbose integer from 0 to 2. Higher values print more information during model fit - for debugging
#' @param backend Estimation backend: \code{"stan"} (the default), CoTiMA2's
#'   experimental \code{"julia"} backend, or ctsem's experimental native
#'   \code{"ctsemjulia"} backend for ML or Bayesian estimation.
#' @param compute Optional execution plan from \code{ctmaCompute()}.
#' @param backendArgs Named list of controls for the selected backend. CoTiMA2
#'   Julia controls follow \code{\link{ctmaJuliaFit}}. For
#'   \code{backend="ctsemjulia"}, use \code{inits}, \code{optimcontrol},
#'   \code{priors}, and \code{verbose}. For Bayesian estimation, use
#'   \code{backendArgs$sampleControl} and choose \code{sampler="nuts"},
#'   \code{"saem"}, \code{"map"}, or \code{"smc"}, or omit sampler for
#'   ctsem's automatic choice. Set \code{backendArgs$intoverpop} explicitly
#'   to choose joint (\code{"none"}) or marginal (\code{"laplace"} or
#'   \code{"augmented"}) sampling, or keep \code{"auto"}. Common controls include \code{chains},
#'   \code{warmup}, \code{draws}, \code{maxDraws}, \code{processes}, and
#'   \code{seed}. SMC instead uses \code{particles} and \code{coreset}; omit
#'   chains and warmup. Placement is delegated to ctsem unless explicitly
#'   supplied. Primary studies are processed
#'   sequentially with either Julia backend. Calling
#'   \code{ctmaBackendDefaults()} without arguments displays Stan defaults;
#'   use \code{ctmaBackendDefaults("ctsemjulia", "ML", "ctmaInit")} for this
#'   workflow or replace \code{"ML"} by \code{"Bayes"} for sampler controls.
#'   workflow. It explicitly returns
#'   \code{optimcontrol$uncertainty = "hessian"}, matching ctsem's native Julia
#'   ML default and providing standard errors and intervals. Fast development
#'   diagnostics can skip this post-fit work with
#'   \code{backendArgs = list(optimcontrol = list(estonly = TRUE))}. Stan has
#'   no separate \code{uncertainty} entry in its backend-default list.

#'
#' @importFrom crayon red blue
#' @importFrom parallel detectCores
#' @importFrom ctsem ctDeintervalise ctLongToWide ctIntervalise ctWideToLong ctModel ctFit ctExtract ctCollapse
#' @importFrom utils read.table write.table packageDescription
#' @importFrom openxlsx addWorksheet writeData createWorkbook openXL saveWorkbook
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster
#' @importFrom foreach %dopar%
#' @importFrom stats cov2cor
#' @importFrom OpenMx expm
#'
#' @export ctmaInit
#'
#' @examples
#' # Fit a ctsem model to all three primary studies summarized in
#' # CoTiMAstudyList_3 and save the three fitted models
#' \dontrun{
#' CoTiMAInitFit_3 <- ctmaInit(primaryStudies=CoTiMAstudyList_3,
#'                             n.latent=2,
#'                             checkSingleStudyResults=FALSE,
#'                             activeDirectory="/Users/tmp/") # adapt!
#' summary(CoTiMAInitFit_3)
#' }
#'
#' @details
#' Random-intercept models with at most two observed waves per person receive
#' a design warning: unrestricted within-person dynamics and random-intercept
#' covariance generally require additional identifying constraints. Unequal
#' intervals and model restrictions can change identification; the warning
#' does not alter the model. Native ctsem no-width uncertainty flags are
#' retained as NA in primary-study parameter tables and combined summaries.
#'
#' The returned \code{studyDesign} records original study numbers, observed
#' wave-count ranges and subject counts. \code{studyLists$twoWaves} and
#' \code{studyLists$threeOrMoreWaves} are ctmaPrep-compatible subsets of the
#' input, preserving order and metadata. Empty groups are NULL. Grouping uses
#' the maximum observed waves per person; mixed observation counts are shown
#' in \code{minWaves} and \code{maxWaves}. Fewer-than-two-wave or unknown
#' designs remain in the complete result and are identified in the overview.
#' No study is excluded from estimation. Cached subset Excel workbooks contain
#' only the available subset summary.
#'
#' @return ctmaFit returns a list containing some arguments supplied, the fitted models, different elements summarizing the main results,
#' model type, and the type of plot that could be performed with the returned object. The arguments in the returned object are activeDirectory,
#' coresToUse, n.latent, n.manifest, and primaryStudyList. The study count is returned as n.studies, the created matrix of loadings of
#' manifest on latent factors is returned as lambda, and a re-organized list of primary studies with some information ommited is returned as
#' studyList. The fitted models for each primary study are found in studyFitList, which is a large list with many elements (e.g., the ctsem
#' model specified by CoTiMA, the rstan model created by ctsem, the fitted rstan model etc.). Further results returned are emprawList
#' (containing the synthesised raw data created), statisticsList (comprising baisc stats such as average sample size, no. of measurement points,
#' etc.), a list with modelResults (i.e., DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef,
#' CINT=model_Cint_Coef), and the paramter names internally used. The summary list,  which is printed if the summary function is applied to the
#' returned object, comprises "estimates" (the aggregated effects), possible randomIntercepts,confidenceIntervals, the
#' minus2ll value and its n.parameters, and possible warning messages (message). Plot type is plot.type=c("drift") and model.type="stanct"
#' ("omx" was deprecated).
#' With \code{fit=FALSE}, \code{ctModelText} contains executable R code that
#' reconstructs the returned \code{ctModel}, including its exact parameter
#' table, binary/ordinal manifest types, individual-variation settings, and TI
#' predictor effects. Evaluate it with \code{eval(parse(text = x$ctModelText))}
#' or copy it into an R script.
#'
#' @param stationary Logical, default FALSE. Start the latent within-person
#'   process in its stationary distribution. This forwards ctsem's stationary
#'   initialization option; it does not fix observed or within-person variances
#'   to one. In CoTiMA2 use the native \code{backend="ctsemjulia"} or a
#'   compatible Stan ctsem version; the bundled Julia backend is not supported.
#'   If omitted, an explicitly supplied legacy/backend stationary setting is
#'   retained. Conflicting explicit settings are errors.
#' @param stationaryVariance Optional fixed positive stationary process variance.
#'   Currently supports one Gaussian latent process with zero intercepts and no
#'   person-level random effects, using backend='ctsemjulia', stationary=TRUE,
#'   determinedDiffusions=TRUE. Set to 1 for a unit-variance
#'   process. Study-specific DRIFT is supported; diffusion is derived as
#'   Q=-2*DRIFT*stationaryVariance. This is distinct from the default estimated
#'   stationary covariance. Do not select separate DIFFUSION/T0VAR study effects.
#' @param determinedDiffusions Logical, default FALSE. Experimental ML/Bayesian
#'   stable parameterization deriving diffusion covariance as
#'   \eqn{Q=-(A\Sigma+\Sigma A^T)} from drift and an estimated stationary
#'   within-person covariance. Currently requires \code{stationary=TRUE},
#'   one or two unrestricted
#'   Gaussian latent processes. Use ctmaFit for native study random DRIFT effects.
#'   Random manifest intercepts are allowed with \code{indVaryingT0=FALSE}.
#'   Their native z covariance construction and correlation safeguards are
#'   preserved; unsupported legacy covariance constructions are rejected.
#'   Fixed/equal/custom process constraints, phantom factors, old raw inits,
#'   and prior equivalence to unrestricted diffusion coordinates are not supported. The free
#'   parameter count is not reduced merely by changing parameterization.
#'   The original drift parameters and their exact transforms are retained,
#'   including the usual negative-diagonal constraints. Covariance coordinates
#'   describe total stationary within-person variance and admissible shape.
#'   Under identical stationarity assumptions this preserves the model domain;
#'   numerical optimization can nevertheless reach different solutions.
#'   For pooled common processes, use \code{allInvModel=TRUE} or, in
#'   CoTiMA2, explicitly choose compatible \code{studyVarying}. Unsupported
#'   restrictions are rejected instead of being silently dropped.
#'   Published ctsem 3.11.1 disables its stationary switch. For the supported
#'   determined model with CINT fixed to zero, an exact compatibility path
#'   explicitly sets the estimated initial covariance and zero mean instead.
#'   Binary/ordinal measurements, random CINT and person-level process random effects
#'   are not yet supported by the determined option.
#'
#'   Bayesian sampling is supported by \code{backend="ctsemjulia"} with \code{optimize=FALSE, priors=TRUE}. Native ctsem priors apply directly to raw drift, log stationary variance/shape and study-contrast coordinates. They are not the transformed priors of an unrestricted DIFFUSION model; no additional Jacobian is added. Native study-level random DRIFT effects are mapped to the shared auxiliary coordinates; derived DIFFUSION/T0VAR cannot have independent random effects. Joint NUTS is requested with \code{backendArgs=list(intoverpop="none", sampleControl=list(sampler="nuts"))}.
#' @param ... Reserved. Unknown arguments are rejected with a migration hint where applicable.
ctmaInit <- function(

    activeDirectory=NULL,
    binaries=NULL,
    chains=NULL,
    checkSingleStudyResults=FALSE,
    cint=0,
    coresToUse=c(2),
    CoTiMAStanctArgs=NULL,
    customPar=FALSE,
    diff=NULL,
    digits=4,
    doPar=1,
    drift=NULL,

    finishsamples=NULL,
    fit=TRUE,
    indVarying=FALSE,
    indVaryingT0=NULL,
    iter=NULL,
    lambda=NULL,
    loadSingleStudyModelFit=c(),
    manifestMeans=0,
    manifestVars=NULL,
    n.latent=NULL,
    n.manifest=0,
    optimize=TRUE,
    primaryStudies=NULL,
    priors=FALSE,
    randomIntercepts=FALSE,
    sameInitialTimes=FALSE,
    saveRawData=list(),
    saveSingleStudyModelFit=c(),
    scaleTI=NULL,
    scaleTime=NULL,
    silentOverwrite=FALSE,
    T0means=0,
    T0var='auto',
    useSV=FALSE,
    verbose=0,
    backend=c("stan", "julia", "ctsemjulia"),
    backendArgs=list(),
    stationary=FALSE,
    determinedDiffusions=FALSE,
    stationaryVariance=NULL,
    compute=NULL, model=NULL, ...)

{  # begin function definition (until end of file)
  .ctmaCheckRetiredArguments(match.call(expand.dots = FALSE), "ctmaInit")
  if (!is.null(model)) return(.ctmaRunConvertedModel(model,
    match.call(expand.dots = TRUE), parent.frame(), "init"))
  .ctma_job <- .ctmaComputeEntry(compute, "ctmaInit", match.call(), parent.frame(), backend)
  if (!is.null(.ctma_job)) return(.ctma_job)
  backendArgs <- .ctmaComputeBind(compute, backend, backendArgs, "ctmaInit")
  if (!is.null(backendArgs) && !is.list(backendArgs)) stop(
    "backendArgs must be a fully named list.",call.=FALSE)
  .ctma_stationary <- .ctmaStationaryOptions(
    stationary, determinedDiffusions, missing(stationary),
    backend_stationary = if (!is.null(backendArgs$stationary)) backendArgs$stationary else
      CoTiMAStanctArgs$stationary,
    backend = tolower(backend[1L]), optimize = optimize,
    priors = if (!is.null(backendArgs$priors)) backendArgs$priors else priors,
    inits = backendArgs$inits, stationaryVariance = stationaryVariance )
  stationary <- .ctma_stationary$stationary
  determinedDiffusions <- .ctma_stationary$determinedDiffusions

  .ctmaVersionNote()
  .ctma_timer <- .ctmaRegisterExecutionTimer("ctmaInit")
  on.exit(.ctma_timer(), add = TRUE)
  .ctma_wave_warning_seen <- FALSE
  .ctma_restore_threads <- .ctmaLimitThreadEnvironment()
  on.exit(.ctma_restore_threads(), add = TRUE)

  .ctmaCheckCtsemJuliaRequest(
    backend, "ctmaInit",
    supported = identical(tolower(backend[1L]), "ctsemjulia")
  )

  .ctmaCheckBundledPriorRequest(backend, optimize, priors)

  requested_estimation <- if (isTRUE(as.logical(optimize)[1L])) {
    "ML"
  } else {
    "Bayes"
  }
  backend_resolution <- .ctmaResolveFitBackend(
    backend = backend,
    backendArgs = backendArgs,
    legacy_stan_args = CoTiMAStanctArgs,
    estimation = requested_estimation,
    workflow = "ctmaInit"
  )
  backend <- backend_resolution$backend
  resolved_backend_args <- backend_resolution$backend_args
  if (identical(backend, "ctsemjulia")) resolved_backend_args$stationary <- stationary
  CoTiMAStanctArgs <- backend_resolution$stan_args
  if (!is.null(CoTiMAStanctArgs)) CoTiMAStanctArgs$stationary <- stationary
  if (identical(backend, "ctsemjulia") &&
      identical(requested_estimation, "Bayes")) {
    if (!is.null(chains) &&
        is.null(backendArgs$sampleControl$chains)) {
      resolved_backend_args$sampleControl$chains <- as.integer(chains)[1L]
    }
    if (!is.null(iter)) {
      if (is.null(backendArgs$sampleControl$warmup)) {
        resolved_backend_args$sampleControl$warmup <-
          as.integer(iter %/% 2L)
      }
      if (is.null(backendArgs$sampleControl$draws)) {
        resolved_backend_args$sampleControl$draws <-
          as.integer(iter - iter %/% 2L)
      }
    }
  }
  if (identical(backend, "julia")) {
    resolved_backend_args <- .ctmaInitJuliaFinishsamples(
      finishsamples = finishsamples,
      backendArgs = backendArgs,
      resolved_backend_args = resolved_backend_args
    )
  }
  if (backend %in% c("julia", "ctsemjulia") && doPar > 1L) {
    stop(
      "ctmaInit(backend='", backend, "') currently requires doPar = 1. ",
      "This prevents nested R/Julia parallelism while the native route is ",
      "being validated.",
      call. = FALSE
    )
  }


  #start.time <- Sys.time()

  #original.options <- options("scipen"); original.options
  original.options <- getOption("scipen"); original.options
  options(scipen = 999); options("scipen") # turn scientific notation off.
  on.exit(options(scipen = original.options), add = TRUE)  # restore and retain timing handler

  catchWarnings <- FALSE
  if (catchWarnings) { # function definition to handle possible errors and warnings during fitting
    run_ctStanFit_logged <- function(expr) {
      flag_warn_hessinv <- FALSE
      warn_msgs <- character(0)

      out <- tryCatch(
        withCallingHandlers(
          expr,
          warning = function(w) {
            msg <- conditionMessage(w)
            warn_msgs <<- c(warn_msgs, msg)

            if (grepl("Hessian inversion failed", msg, fixed = TRUE)) {
              flag_warn_hessinv <<- TRUE
            }

            #invokeRestart("muffleWarning")  # optional: Konsole sauber halten
          }
        ),
        error = function(e) {
          # Error sauber zurückgeben statt Abbruch
          err_msg <- conditionMessage(e)

          # optional: "inits = c(...)" aus dem Text ziehen (wenn ctsem das ausgibt)
          inits_txt <- NA_character_
          m <- regexpr("inits\\s*=\\s*c\\([^\\)]*\\)", err_msg)
          if (m[1] != -1) inits_txt <- regmatches(err_msg, m)

          return(list(
            ok = FALSE,
            fit = NULL,
            warn_hessinv = flag_warn_hessinv,
            warnings = warn_msgs,
            error = err_msg,
            inits_suggestion = inits_txt
          ))
        }
      )

      # Wenn kein Error: out ist der Fit
      if (!is.list(out) || isTRUE(inherits(out, "ctFit"))) {
        return(list(
          ok = TRUE,
          fit = out,
          warn_hessinv = flag_warn_hessinv,
          warnings = warn_msgs,
          error = NA_character_,
          inits_suggestion = NA_character_
        ))
      }

      # Wenn error-handler schon list() geliefert hat:
      out
    }
  }


  #######################################################################################################################
  ########################################### Check Model Specification #################################################
  #######################################################################################################################
  {
    Msg <- "################################################################################# \n########################## Check Model Specification ############################ \n#################################################################################"
    message(Msg)

    if (utils::packageDescription("ctsem")$Version > "3.10.2") type <- "ct" else type <- "stanct"

    if (  (indVarying == "Cint") | (indVarying == "cint")) indVarying <- "CINT"
    if ( (indVarying == "Manifest") | (indVarying == "manifest")) indVarying <- "MANIFEST"
    if (indVarying == TRUE) indVarying <- "MANIFEST"
    #
    if ( (randomIntercepts == 'cint')  | (randomIntercepts == 'Cint')) randomIntercepts <- "CINT"
    if ( (randomIntercepts == "Manifest") | (randomIntercepts == "manifest") ) randomIntercepts <- "MANIFEST"
    if (randomIntercepts == TRUE)  randomIntercepts <- "MANIFEST"
    # For separately fitted primary studies, CoTiMA's pseudo-process random
    # intercept representation supplies the Julia-compatible formulation of
    # individually varying manifest means or continuous-time intercepts. This
    # translation is intentionally confined to ctmaInit(); in pooled models
    # the two options imply different covariance structures across studies.
    iv_resolution <- .ctmaInitJuliaIndividualVariation(
      indVarying = indVarying,
      randomIntercepts = randomIntercepts,
      backend = backend,
      backendArgs = backendArgs,
      resolved_backend_args = resolved_backend_args
    )
    indVarying <- iv_resolution$indVarying
    randomIntercepts <- iv_resolution$randomIntercepts
    resolved_backend_args <- iv_resolution$backend_args
    randomInterceptsSettings <- randomIntercepts
    #
    if ((randomIntercepts == "MANIFEST") | (randomIntercepts == "CINT")) {
      indVarying <- FALSE
    }
    err <- FALSE
    if ( (!(randomIntercepts %in% c("MANIFEST", "CINT", FALSE))) |
         (!(indVarying %in% c("MANIFEST", "CINT", FALSE))) ) {
      err <- TRUE
    }
    if (err){

      ErrorMsg <- "\n The arguments \"indVarying\" and \"randomIntercepts\" have to be TRUE or FALSE or \"MANIFEST\" or \"CINT\". "
      stop(ErrorMsg)
    }
    #}


    if (is.null(verbose)) tmp1 <- 1 else tmp1 <- 0
    if (is.null(verbose) & (optimize == FALSE) )  {verbose <- 0}
    if (tmp1 == 1) {verbose <- .ctmaStanBackendDefaults()$verbose}

    if (is.null(primaryStudies)) {

      cat(crayon::red$bold(" List with lists of primary study information not specified!", sep="\n"))
      cat(crayon::red$bold(" ", " ", sep="\n"))
      cat(crayon::red$bold("Should I try to get primary study information (e.g., empcov1, delta_t22) from the global environment?", sep="\n"))
      cat(crayon::red$bold(" ", " ", sep="\n"))
      cat(crayon::red$bold("(This is the way how earlier version of CoTiMA were run but no longer recommended.)", sep="\n"))
      cat(crayon::red$bold(" ", " ", sep="\n"))
      cat(crayon::blue("Press 'q' to quit and specify or 'c' to give it a try. Press ENTER afterwards ", "\n"))
      char <- readline(" ")
      while (!(char == 'c') & !(char == 'C') & !(char == 'q') & !(char == 'Q')) {
        cat((crayon::blue("Please press 'q' to quit and specify primaryStudies or 'c' to try reading from global environment Press ENTER afterwards.", "\n")))
        char <- readline(" ")
      }
      if (char == 'q' | char == 'Q') {
        ErrorMsg <- "\nGood luck for the next try!"
        stop(ErrorMsg)
      } else {
        if (!(is.numeric(tryCatch(get(paste("sampleSize", 1, sep = "")), error = function(e) e)))) {
          ErrorMsg <- "\nGetting primary study information from the global environment failed \nTo test, I searched for sampleSize1 and could not find it! \nGood luck for the next try!"
          stop(ErrorMsg)
        } else {
          cat(crayon::blue("Please type the number of primary studies to read from global environment. Press ENTER afterwards ", "\n"))
          char <- as.numeric(readline(""))
          while ((is.na(char))) {
            cat((crayon::blue("Please type the number of primary studies to read from global environment. Press ENTER afterwards ", "\n")))
            char <- as.numeric(readline(""))
          }
          primaryStudies <- ctmaPrep(selectedStudies=1:char)
        } # END else (catch failed)
      } # END else
    } #


    if (is.null(n.latent)) {

      ErrorMsg <- "\nNumber of variables (n.latent) not specified! \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    if (is.null(binaries)) binaries <- rep(0, max(n.manifest, n.latent))
    if (length(binaries) != max(n.manifest, n.latent)) {

      ErrorMsg <- "\nThe number of binaries provided is incorrect! \nGood luck for the next try!"
      stop(ErrorMsg)
    }


    if (is.null(activeDirectory)) {

      ErrorMsg <- "\nNo active directory has been specified! \nGood luck for the next try!"
      stop(ErrorMsg)
    }

    resource_resolution <- .ctmaResolveCores(
      coresToUse, backend = backend, workflow = "ctmaInit"
    )
    coresToUse <- resource_resolution$resolved

    if (identical(backend, "julia")) {
      julia_resources <- .ctmaInitJuliaResources(
        cores = coresToUse,
        estimation = requested_estimation,
        backendArgs = backendArgs,
        resolved_backend_args = resolved_backend_args,
        workflow = "ctmaInit"
      )
      resolved_backend_args <- julia_resources$backend_args
    }

    # CHD changed on 19 Sep 2022
    if (doPar > 1) {
      myCluster <- parallel::makeCluster(coresToUse)
      on.exit(parallel::stopCluster(myCluster), add = TRUE)
      doParallel::registerDoParallel(myCluster)
      '%dopar%' <- foreach::'%dopar%'
    }

    if (n.manifest > n.latent ) {
      if (is.null(lambda)) {

        ErrorMsg <- "\nManifest variables specified, but not matrix of loadings (lambda) specified, which is required. \nGood luck for the next try!"
        stop(ErrorMsg)
      } else {
        if ( (dim(lambda)[1] != n.manifest) | (dim(lambda)[2] != n.latent) ) {

          ErrorMsg <- "\nDimensions of loadings matrix (lambda) do not match n.latent and n.manifest. \nGood luck for the next try!"
          stop(ErrorMsg)
        }
      }
    }

    if ((n.manifest <= n.latent ) &  (is.null(lambda)))  lambda=diag(n.latent)


    { # fitting params

      invariantDrift <- FALSE
      moderatedDrift <- NULL

      CoTiMAStanctArgs <- .ctmaResolveInitStanArgs(
        overrides = CoTiMAStanctArgs,
        scaleTI = scaleTI,
        scaleTime = scaleTime,
        optimize = optimize,
        priors = priors,
        finishsamples = finishsamples,
        chains = chains,
        iter = iter,
        verbose = verbose
      )
      CoTiMAStanctArgs$stationary <- stationary
    }

    # correction of misspecified user input
    if (is.null(T0means)) T0means <- 0
    if (is.null(manifestMeans)) manifestMeans <- 0


  } ### END Check Model Specification ###


  #######################################################################################################################
  ##################### Read user provided data and create list with all study information ##############################
  #######################################################################################################################
  {
    Msg <- "################################################################################# \n###### Read user provided data and create list with all study information ####### \n#################################################################################"
    message(Msg)

    studyList <- list()
    minInterval <- .0001 # replace missing values for time lags (cannot be NA because it is a definition variable)
    maxLengthModeratorVector <- 0
    loadRawDataStudyNumbers <- unlist(lapply(primaryStudies$rawData,
                                             function(extract) extract$studyNumbers)); loadRawDataStudyNumbers

    # resorting primary study information
    if (!(exists("moderatorNumber"))) moderatorNumber <- 1; moderatorNumber
    # determine number of studies (if list is created by ctmaInit.R it is 1 element too long)
    # Option 1: Applies for older versions of ctmaPrep
    tmp <- unlist(primaryStudies$studyNumbers); tmp
    tmp2 <- which(is.na(tmp)); tmp2
    for (i in tmp2) primaryStudies$studyNumbers[[i]] <- NULL
    tmp <- length(tmp[!(is.na(tmp))]); tmp
    if  ( is.na(primaryStudies$deltas[tmp]) &
          is.na(primaryStudies$sampleSizes[tmp]) &
          is.na(primaryStudies$deltas[tmp]) &
          is.null(dim(primaryStudies$pairwiseNs[tmp])) &
          is.null(dim(primaryStudies$emcovs[tmp])) &
          all(is.na(unlist(primaryStudies$moderators[tmp]))) &
          is.null(primaryStudies$rawData$fileName[tmp])
          # 9. Aug. 2022:
          & (length(tmp) > 1) ) {
      if (tmp != 1) n.studies <- tmp-1 else n.studies <- 1
      primaryStudies$studyNumbers[[tmp]] <- NULL
    } else {
      n.studies <- tmp
    }
    # Option 2: versions of ctmaPrep after 14. August 2020
    if (!(is.null(primaryStudies$n.studies))) n.studies <- primaryStudies$n.studies; n.studies

    # delete empty list entries
    if (n.studies > 1) { # may not apply if a single study is fitted (e.g., ctmaOptimizeInit) # RECENT CHANGE
      tmp1 <- which(names(primaryStudies) == "n.studies"); tmp1
      for (i in 1:(tmp1-1)) primaryStudies[[i]][[n.studies+1]] <- NULL
    }

    for (i in 1:n.studies) {
      studyList[[i]] <- list(studyNumber=i, empcov=primaryStudies$empcovs[[i]], delta_t=primaryStudies$deltas[[i]],
                             sampleSize=primaryStudies$sampleSizes[[i]], originalStudyNo=primaryStudies$studyNumber[[i]],
                             timePoints=sum(length(primaryStudies$deltas[[i]]), 1), moderators=primaryStudies$moderators[[i]],
                             maxModerators=length(primaryStudies$moderators[[i]]), startValues=primaryStudies$startValues[[i]],
                             rawData=primaryStudies$rawData[[i]], pairwiseN=primaryStudies$pairwiseNs[[i]],
                             source=paste(primaryStudies$source[[i]], collapse=", "))
      if (useSV == FALSE) studyList[[i]]$startValues <- NULL
      if (is.null(studyList[[i]]$startValues)) studyList[[i]]$startValues <- NA


      if (length(primaryStudies$moderators[[i]]) > maxLengthModeratorVector) maxLengthModeratorVector <- length(primaryStudies$moderators[[i]])
      # check matrix symmetry if matrix is provided
      if (!(primaryStudies$studyNumbers[i] %in% loadRawDataStudyNumbers)) {
        if (isSymmetric(primaryStudies$empcovs[[i]]) == FALSE) {
          ErrorMsg <- paste0("\nThe correlation matrix of study no.", i, " is not symmetric. Check and re-start! \nGood luck for the next try!")
          stop(ErrorMsg)
        }
      }
    }

    ### synthesise pseudo raw data for all studies or load raw data if available & specified
    empraw <- lags <- moderators <- emprawMod <- allSampleSizes <- lostN <- overallNDiff <- relativeNDiff <- list()
    emprawLong <- list()
    empraw.ind.mod <- list() # CHD 19.6.2023

    if (n.studies > 1) {
      if (is.na(primaryStudies$deltas[length(primaryStudies$deltas)])) primaryStudies$deltas <- primaryStudies$deltas[-length(primaryStudies$deltas)]
    }
    if (is.na(primaryStudies$deltas[length(primaryStudies$deltas)])) primaryStudies$deltas[length(primaryStudies$deltas)] <- NULL
    allTpoints <-unlist(lapply(primaryStudies$deltas, function(extract) length(extract)+1)); allTpoints # more complicated than expected

    # manifests
    if (n.manifest > n.latent) {
      manifestNames <- paste0("y", 1:n.manifest); manifestNames
      latentNames <- paste0("V", 1:n.latent); latentNames
    } else {
      manifestNames <- paste0("V", 1:n.latent); manifestNames
      latentNames <- paste0("V", 1:n.latent); latentNames
    }

    for (i in 1:n.studies) {
      #i <- 1
      if (!(studyList[[i]]$originalStudyNo %in% loadRawDataStudyNumbers)) {
        currentSampleSize <- (lapply(studyList, function(extract) extract$sampleSize))[[i]]; currentSampleSize
        currentTpoints <- (lapply(studyList, function(extract) extract$timePoints))[[i]]; currentTpoints
        currentEmpcov <- (lapply(studyList, function(extract) extract$empcov))[[i]]; currentEmpcov
        currentLags <- (lapply(studyList, function(extract) extract$delta_t))[[i]]; currentLags
        currentPairwiseN <- (lapply(studyList, function(extract) extract$pairwiseN))[[i]]; currentPairwiseN
        currentModerators <- (lapply(studyList, function(extract) extract$moderators))[[i]]; currentModerators

        currentVarnames <- c()
        for (j in 1:(currentTpoints)) {
          if (n.manifest == 0) {
            for (h in 1:n.latent) {
              currentVarnames <- c(currentVarnames, paste0("V",h,"_T", (j-1)))
            }
          } else {
            for (h in 1:n.manifest) {
              currentVarnames <- c(currentVarnames, paste0("y",h,"_T", (j-1)))
            }
          }
        }

        # Synthesise Pseudo Raw Data
        if (!(studyList[[i]]$originalStudyNo %in% loadSingleStudyModelFit)) {
          tmp1 <- paste0(" Create Pseudo Raw Data for Study No. ", i, ".    Could take long !!! ")
          tmp2 <- nchar(tmp1); tmp2
          tmp3 <- (81 - tmp2)/2; tmp3
          tmp4 <- strrep("#", round(tmp3 + 0.45, 0)); tmp4
          tmp5 <- strrep("#", round(tmp3 - 0.45, 0)); tmp5
          tmp6 <- paste0(tmp4, tmp1, tmp5); tmp6
          Msg <- paste0("################################################################################# \n", tmp6, "\n#################################################################################")
          message(Msg)
        }


        # reduce computation time by synthesising pseudo Raw Data with small sample size because the data are loaded anyway
        currentSampleSizeTmp <- currentSampleSize
        currentPairwiseNTmp <- currentPairwiseN
        currentEmpcovTmp <- currentEmpcov
        # CHD 3.12.2024
        x1 <- paste0(activeDirectory, loadSingleStudyModelFit[1], " singleStudyFits/",loadSingleStudyModelFit[1], " studyFit", studyList[[i]]$originalStudyNo, ".rds"); x1
        #if (length(loadSingleStudyModelFit) > 0) {
        if ( (length(loadSingleStudyModelFit) > 0) & exists(x1) ) {
          tmp <- as.numeric(loadSingleStudyModelFit[2:length(loadSingleStudyModelFit)]); tmp
          if (studyList[[i]]$originalStudyNo %in% tmp) {
            currentSampleSizeTmp <- (n.latent * currentTpoints)^2; currentSampleSizeTmp
            currentPairwiseNTmp <- matrix(0,0,0)
            currentEmpcovTmp <- matrix(0, n.latent * currentTpoints, n.latent * currentTpoints); currentEmpcovTmp
          }
        }

        tmp <- suppressWarnings(ctmaPRaw(empCovMat=currentEmpcovTmp, empN=currentSampleSizeTmp, empNMat=currentPairwiseNTmp))

        empraw[[i]] <- tmp$data
        lostN[[i]] <- tmp$lostN
        overallNDiff[[i]] <- tmp$overallLostN
        relativeNDiff[[i]] <- tmp$relativeLostN
        lags[[i]] <- matrix(currentLags, nrow=dim(empraw[[i]])[1], ncol=currentTpoints-1, byrow=TRUE)
        empraw[[i]] <- cbind(empraw[[i]], lags[[i]]);
        colnames(empraw[[i]]) <- c(c(currentVarnames, paste0("dT", seq(1:(currentTpoints-1)))))
        empraw[[i]] <- as.data.frame(empraw[[i]])

        # ADDED TO DEAL WITH MANIFEST VARIABLES
        n.var <- max(c(n.manifest, n.latent)); n.var
        ## START correction of current lags if entire time point is missing for a case
        emprawLongTmp <- ctsem::ctWideToLong(empraw[[i]], Tpoints=currentTpoints, n.manifest=n.var, manifestNames=manifestNames)
        emprawLongTmp <- suppressMessages(ctsem::ctDeintervalise(datalong = emprawLongTmp, id='id', dT='dT'))
        ## eliminate rows where ALL latents are NA
        if (n.manifest > n.latent) targetCols <- paste0("y", 1:n.manifest) else targetCols <- paste0("V", 1:n.latent)
        emprawLongTmp <- emprawLongTmp[, ][ apply(emprawLongTmp[, targetCols, drop=FALSE], 1, function(x) sum(is.na(x)) != n.var ), ]
        # eliminate rows where time is NA
        emprawLongTmp <- emprawLongTmp[which(!(is.na(emprawLongTmp[, "time"]))), ]
        # make wide format
        emprawWide <- suppressMessages(ctsem::ctLongToWide(emprawLongTmp, id='id', time='time', manifestNames=manifestNames))
        # intervalise
        emprawWide <- suppressMessages(ctsem::ctIntervalise(emprawWide, Tpoints=currentTpoints, n.manifest=n.var, manifestNames=manifestNames))
        # restore
        empraw[[i]] <- as.data.frame(emprawWide)
        # END correction
        # add potential moderators
      }

      # load raw data on request
      preProcessed <- FALSE
      if (studyList[[i]]$originalStudyNo %in% loadRawDataStudyNumbers) {
        if ( (!(is.null(primaryStudies$emprawList[[i]]))) &
             (
               (isTRUE(primaryStudies$reUseEmprawData) |
                (is.null(primaryStudies$empcovs[[i]]))
                |
                (length(primaryStudies$empcovs[[i]]) <= 1) # could be NA or 0x0 matrix
               )
             )
        ) {
          # if the function list of primary studies is already post-processed (ctmaSV) and called from ctmaOptimizeINit)
          preProcessed <- TRUE
          .ctmaValidateCachedStudyWide(primaryStudies$emprawList[[i]],
            studyList[[i]], max(c(n.latent,n.manifest)))
          empraw[[i]] <- primaryStudies$emprawList[[i]]
          if (!(exists("n.var"))) n.var <- max(c(n.latent, n.manifest))
          tmp1 <- dim(empraw[[i]])[2]; tmp1
          currentTpoints <- (tmp1 + 1)/(n.var+1); currentTpoints
          currentTpointsBackup <- currentTpoints # in case function is called from ctmaOptimizeInit
          currentVarnames <- c()
          for (j in 1:(currentTpoints)) {
            if (n.manifest == 0) {
              for (h in 1:n.latent) {
                currentVarnames <- c(currentVarnames, paste0("V",h,"_T", (j-1)))
              }
            } else {
              for (h in 1:n.manifest) {
                currentVarnames <- c(currentVarnames, paste0("y",h,"_T", (j-1)))
              }
            }
          }
          tmp1 <- empraw[[i]]
          tmp2 <- tmp1[,grep("dT", colnames(tmp1))] == minInterval
          tmp1[ ,grep("dT", colnames(tmp1))][tmp2] <- NA
          for (h in 1:(currentTpoints-1)) studyList[[i]]$delta_t[h] <- mean(tmp1[, paste0("dT", h)], na.rm=TRUE)
          #studyList[[i]]$delta_t
        } else {
          if (!(exists("n.var"))) n.var <- max(c(n.latent, n.manifest))
          currentTpoints <- length((lapply(studyList, function(extract) extract$delta_t))[[i]])+1; currentTpoints

          tmp1 <- studyList[[i]]$rawData$fileName; tmp1
          tmp2 <- studyList[[i]]$rawData$header; tmp2
          tmp3 <- studyList[[i]]$rawData$dec; tmp3
          tmp4 <- studyList[[i]]$rawData$sep; tmp4
          if (!(is.null(primaryStudies$activeDirectory))) {
            if (activeDirectory != primaryStudies$activeDirectory) {
              tmp1 <- gsub(primaryStudies$activeDirectory, activeDirectory, tmp1, fixed=TRUE)
              #print(tmp1)
            }
          } else {
            #tmp1 <- gsub(primaryStudies$activeDirectory, activeDirectory, tmp1, fixed=TRUE)
          }
          tmpData <- utils::read.table(file=tmp1,
                                       header=tmp2,
                                       dec=tmp3,
                                       sep=tmp4)
          currentVarnames <- c()
          for (j in 1:(currentTpoints)) {
            if (n.manifest == 0) {
              for (h in 1:n.latent) {
                currentVarnames <- c(currentVarnames, paste0("V",h,"_T", (j-1)))
              }
            } else {
              for (h in 1:n.manifest) {
                currentVarnames <- c(currentVarnames, paste0("y",h,"_T", (j-1)))
              }
            }
          }

          # CHD 19.6.2023 extract possible ind level moderators 28.6.2023
          if ("n.ind.mod" %in% names(studyList[[i]]$rawData)) {
            if (studyList[[i]]$rawData$n.ind.mod != 0) {
              if (!(exists("n.var"))) n.var <- max(c(n.latent, n.manifest))
              targetCols <- (n.var * currentTpoints + currentTpoints -1 + 1): ncol(tmpData); targetCols
              empraw.ind.mod[[i]] <- list()
              empraw.ind.mod[[i]] <- tmpData[, targetCols]
              empraw.ind.mod[[i]][empraw.ind.mod[[i]] == studyList[[i]]$rawData$missingValues] <- NA
              tmpData <- tmpData[, -targetCols]
            } else {
              empraw.ind.mod[[i]] <- NA
            }
          }
          #}


          # replace missing values
          tmpData <- as.matrix(tmpData) # important: line below will not work without having data as a matrix
          tmpData[tmpData %in% studyList[[i]]$rawData$missingValues] <- NA
          empraw[[i]] <- as.data.frame(tmpData)

          ## START correction of current lags if entire time point is missing for a case
          # if called from ctmaOptimize
          if (preProcessed == TRUE) {
            if (!(exists("n.var"))) n.var <- max(n.latent, n.manifest)
            # change variable names
            tmp1 <- dim(empraw[[i]])[2]; tmp1
            currentTpoints <- (tmp1 + 1)/(n.var+1); currentTpoints
          }

          # CHD 19.6.2023 this is an error check (n.ind.mod set to 0 but data set includes TI at the end)
          #if ((!is.na(empraw.ind.mod[[i]]))  & (currentTpoints != round(currentTpoints, 0))) {
          if (currentTpoints != round(currentTpoints, 0)) {

            ErrorMsg <- "\nI have problems with the raw data set. Possibly you specified n.ind.mod incorrectly before doing ctmaPrep. The n.ind.mod should by > 0 if the raw data set does not have the last time interval as last column.\nGood luck for the next try!"
            stop(ErrorMsg)
          }

          currentTpointsBackup <- currentTpoints # in case function is called from ctmaOptimizeInit
          if (n.manifest > n.latent) {
            colnames(empraw[[i]])[1:(currentTpoints * n.manifest)] <- paste0(paste0("x", 1:n.manifest), "_T", rep(0:(currentTpoints-1), each=n.manifest))
          } else {
            colnames(empraw[[i]])[1:(currentTpoints * n.latent)] <- paste0(paste0("V", 1:n.latent), "_T", rep(0:(currentTpoints-1), each=n.latent))
          }

          # wide to long
          emprawLongTmp <- ctsem::ctWideToLong(empraw[[i]], Tpoints=currentTpoints, n.manifest=n.var, manifestNames=manifestNames)
          emprawLongTmp <- suppressMessages(ctsem::ctDeintervalise(datalong = emprawLongTmp, id='id', dT='dT'))

          # eliminate rows where time is NA
          emprawLongTmp <- emprawLongTmp[which(!(is.na(emprawLongTmp[, "time"]))), ]
          # 9. Aug. 2022: handle possibly reduced Tpoints (if some rows were eliminated for all cases)
          currentTpoints <- max(table(emprawLongTmp[,1])); currentTpoints

          # make wide format
          emprawWide <- suppressMessages(ctsem::ctLongToWide(emprawLongTmp, id='id', time='time', manifestNames=manifestNames))
          # intervalise
          emprawWide <- suppressMessages(ctsem::ctIntervalise(emprawWide,
                                                              Tpoints=currentTpoints,
                                                              n.manifest=n.var,
                                                              manifestNames=manifestNames,
                                                              mininterval=minInterval))
          # restore
          empraw[[i]] <- as.data.frame(emprawWide)

          # Change the NAs provided for deltas if raw data are loaded
          for (h in 1:(currentTpoints-1)) {
            # temporarily replace dT = .001 by NA
            tmp1 <- grep("dT", colnames(empraw[[i]])); tmp1
            colnamesTmp <- colnames(empraw[[i]])[tmp1]; colnamesTmp
            emprawTmp <- as.matrix(empraw[[i]][, tmp1])
            colnames(emprawTmp) <- colnamesTmp
            emprawTmp[emprawTmp == minInterval] <- NA
            studyList[[i]]$delta_t[h] <- mean(emprawTmp[, paste0("dT", h)], na.rm=TRUE)
          }


          # change sample size if entire cases were deleted
          #studyList[[i]]$sampleSize <- (dim(empraw[[i]]))[1]
          #CHD 2.1.2025
          studyList[[i]]$sampleSize <- (dim(empraw[[i]]))[1] / mean(table(emprawLongTmp[,1]))
          allSampleSizes[[i]] <- dim(empraw[[i]])[1]; allSampleSizes[[i]]
          currentSampleSize <- (lapply(studyList, function(extract) extract$sampleSize))[[i]]; currentSampleSize
          currentTpoints <- allTpoints[[i]]; currentTpoints
          if (is.null(currentTpoints)) currentTpoints <- currentTpointsBackup

          #currentVarnames
          colnames(empraw[[i]]) <- c(c(currentVarnames, paste0("dT", seq(1:(currentTpoints-1)))))

          # standardize (variables - not time lags) if option is chosen
          if (studyList[[i]]$rawData$standardize == TRUE) empraw[[i]][, currentVarnames] <- scale(empraw[[i]][, currentVarnames])

          # replace missing values for time lags dTi by minInterval (has to be so because dTi are definition variables)
          tmpData <- empraw[[i]][, paste0("dT", seq(1:(currentTpoints-1)))]
          tmpData[is.na(tmpData)] <- minInterval
          empraw[[i]][, paste0("dT", seq(1:(currentTpoints-1)))] <- tmpData

          # add moderators to loaded raw data
          # Save raw data  on request
          if ( i %in% saveRawData$studyNumbers ) {
            x1 <- paste0(saveRawData$fileName, i, ".dat"); x1
            utils::write.table(empraw[[i]], file=x1, row.names=saveRawData$row.names, col.names=saveRawData$col.names,
                               sep=saveRawData$sep, dec=saveRawData$dec)
          }
        }
      }

      # augment synthesised pseudo raw data for stanct model
      {
        dataTmp <- empraw[[i]]
        dataTmp2 <- ctsem::ctWideToLong(dataTmp, Tpoints=currentTpoints, n.manifest=n.var, #n.TIpred = (n.studies-1),
                                        manifestNames=manifestNames)
        dataTmp3 <- suppressMessages(ctsem::ctDeintervalise(dataTmp2))
        dataTmp3[, "time"] <- dataTmp3[, "time"] * CoTiMAStanctArgs$scaleTime

        # eliminate rows where ALL latents are NA
        if (n.manifest > n.latent) {
          dataTmp3 <- dataTmp3[, ][ apply(dataTmp3[, paste0("y", 1:n.manifest), drop=FALSE], 1, function(x) sum(is.na(x)) != n.manifest ), ]
        } else {
          dataTmp3 <- dataTmp3[, ][ apply(dataTmp3[, paste0("V", 1:n.latent), drop=FALSE], 1, function(x) sum(is.na(x)) != n.latent ), ]
        }
        emprawLong[[i]] <- dataTmp3

      }
    } ### END for i ...

    # Check if sample sizes specified in prep file deviate from cases provided in possible raw data files
    N1 <- (unlist((lapply(empraw, function(extract) dim(extract)[1])))); N1
    N2 <- (unlist(primaryStudies$sampleSizes)); N2
    N1 <- sum(unlist((lapply(empraw, function(extract) dim(extract)[1]))) , na.rm=TRUE); N1
    N2 <- sum(unlist(primaryStudies$sampleSizes), na.rm=TRUE); N2
    if (!(N1 == N2)) {
      tmp1 <- unlist(lapply(empraw, function(extract) dim(extract)[1])); tmp1
      tmp2 <- unlist(primaryStudies$sampleSizes); tmp2
    }
  } ### END Read user provided data and create list with all study information ###


  #######################################################################################################################
  ################################################### Some Statistics ###################################################
  #######################################################################################################################
  {
    Msg <- "################################################################################# \n################# Compute Summary Statistics of Primary Studies ################# \n#################################################################################"
    message(Msg)

    ### some stats
    # Sample size
    allSampleSizes <-unlist(lapply(studyList, function(extract) extract$sampleSize)); allSampleSizes
    overallSampleSize <- N1; overallSampleSize
    meanSampleSize <- mean(allSampleSizes, na.rm=TRUE); meanSampleSize
    maxSampleSize <- suppressWarnings(max(allSampleSizes, na.rm=TRUE)); maxSampleSize
    minSampleSize <- suppressWarnings(min(allSampleSizes, na.rm=TRUE)); minSampleSize
    # lags (will be computed again below if raw data are loaded #
    allDeltas <-unlist(lapply(studyList, function(extract) extract$delta_t)); allDeltas
    if (is.na(allDeltas[length(allDeltas)])) allDeltas <- allDeltas[-length(allDeltas)]; allDeltas
    meanDelta <- mean(allDeltas, na.rm=TRUE); meanDelta
    maxDelta <- max(allDeltas, na.rm=TRUE); maxDelta
    minDelta <- min(allDeltas, na.rm=TRUE); minDelta
    # Time points
    allTpoints; length(allTpoints)
    overallTpoints <- sum(allTpoints, na.rm=TRUE); overallTpoints
    meanTpoints  <- mean(allTpoints, na.rm=TRUE); meanTpoints
    maxTpoints  <- max(allTpoints, na.rm=TRUE); maxTpoints
    minTpoints  <- min(allTpoints, na.rm=TRUE); minTpoints
    # Moderators Information
    statisticsList <- list(originalStudyNumbers=unlist(primaryStudies$studyNumber),
                           allSampleSizes=allSampleSizes, overallSampleSize=overallSampleSize, meanSampleSize=meanSampleSize,
                           maxSampleSize=maxSampleSize, minSampleSize=minSampleSize,
                           allDeltas=allDeltas, meanDelta=meanDelta, maxDelta=maxDelta, minDelta=minDelta,
                           allTpoints=allTpoints, overallTpoints=overallTpoints, meanTpoints=meanTpoints,
                           maxTpoints=maxTpoints, minTpoints=minTpoints,
                           n.studies=n.studies)

  } ### END Some Statistics ###

  #######################################################################################################################
  ############################# Create ctsem model template to fit all primary studies ##################################
  #######################################################################################################################
  {
    Msg <- "################################################################################# \n############# Set ctsem Model Parameters to fit all Primary Studies ############# \n#################################################################################"
    message(Msg)

    namedCint <- .ctma_named_location(cint, latentNames, "cint")
    namedManifest <- .ctma_named_location(manifestMeans, manifestNames, "manifestMeans")
    namedT0 <- .ctma_named_location(T0means, latentNames, "T0means")
    if (!is.null(namedCint)) cint <- namedCint
    if (!is.null(namedManifest)) manifestMeans <- namedManifest
    if (!is.null(namedT0)) T0means <- namedT0
    location_specs <- .ctma_normalize_location_specs(
      cint, manifestMeans, n.latent, n.var)
    cint <- location_specs$cint
    manifestMeans <- location_specs$manifestMeans

    namesAndParams <- ctmaLabels(
      n.latent=n.latent,
      n.manifest=n.manifest,
      lambda=lambda,
      drift=drift,
      diff=diff, # CHD 13.3.24
      invariantDrift=invariantDrift,
      moderatedDrift=NULL,
      equalDrift=NULL,
      T0means=T0means,
      manifestMeans=manifestMeans,
      manifestVars=manifestVars
    )
    #
    driftNames <- namesAndParams$driftNames; driftNames
    driftFullNames <- namesAndParams$driftFullNames; driftFullNames
    driftParams <- namesAndParams$driftParams; driftParams
    diffNames <- namesAndParams$diffNames; diffNames
    diffParams <- namesAndParams$diffParams; diffParams
    diffFullNames <- namesAndParams$diffFullNames; diffFullNames
    invariantDriftNames <- namesAndParams$invariantDriftNames; invariantDriftNames
    invariantDriftParams <- namesAndParams$invariantDriftParams; invariantDriftParams
    lambdaParams <- namesAndParams$lambdaParams; lambdaParams
    T0VARParams <- namesAndParams$T0VARParams; T0VARParams
    manifestmeansParams <- namesAndParams$manifestMeansParams; manifestmeansParams
    manifestMeansParams <- namesAndParams$manifestMeansParams; manifestMeansParams
    T0meansParams=namesAndParams$T0meansParams; T0meansParams
    manifestVarsParams <- namesAndParams$manifestVarsParams; manifestVarsParams
    #CHD 13.6.2023
    T0VARParams <- T0var; T0VARParams
    CINTParams <- cint; CINTParams

  } ### END Create ctsem model template to fit all primary studies ###


  #######################################################################################################################
  ##################################### Check Specification of Primary Studies ##########################################
  #######################################################################################################################
  {
    Msg <- "################################################################################# \n#################### Check Specification of Primary Studies ##################### \n#################################################################################"
    message(Msg)

    if (length(saveSingleStudyModelFit) == 1){

      if (silentOverwrite==FALSE) {
        cat(crayon::red$bold("You have indicated that you want to save SingleStudyModelFits, but have not selected any study/studies to save the fit for!","\n"))
        cat(crayon::red("Would you like to save the SingleStudyModelFits for ALL studies??","\n"))
        cat(crayon::blue("Press 'y' to save ALL singleStudyModelFits or 's' to continue and","\n"))
        cat(crayon::blue("select the study/studies you whish to save the fits for. If you wish to quite, press 'q'. Press ENTER afterwards","\n"))
        char <- readline(" ")
        while (!(char == 's') & !(char == 'S') & !(char == 'y') & !(char == 'Y') & !(char == 'q') & !(char == 'Q')) {
          cat((crayon::blue("Please press 'y' to save ALL, 's' to specify the study/studies to save, or 'q' to quit. Press ENTER afterwards.", "\n")))
          char <- readline(" ")
        }
        if (char == 'y' | char == 'Y') {
          for (i in unlist(primaryStudies$studyNumber)) saveSingleStudyModelFit <- c(saveSingleStudyModelFit, i)
        } else if (char == 's' | char == 'S') {
          cat(crayon::blue("Which SingleStudyModelFits would you like to save?", "\n"))
          cat(crayon::blue("Please enter the no. of study/studies separated by commas!", "\n"))
          chars <- readline(" ")
          chars <- gsub(" ", "", chars, fixed = TRUE)
          chars <- unlist(strsplit(chars, ","))
          chars
          saveSingleStudyModelFit <- c(saveSingleStudyModelFit, chars)
        } else if (char == 'q' | char == 'Q'){
          ErrorMsg <- "Good luck for the next try!"
          stop(ErrorMsg)
        }
      } else {
        saveSingleStudyModelFit <- c(saveSingleStudyModelFit, unlist(primaryStudies$studyNumber))
      }
    }

    if (length(loadSingleStudyModelFit) == 1){

      if (silentOverwrite==FALSE) {
        cat(crayon::red$bold("You have indicated that you want to load SingleStudyModelFits, but have not selected any study/studies to load the fit for!","\n"))
        cat(crayon::red("Would you like to load the SingleStudyModelFits for ALL studies??","\n"))
        cat(crayon::blue("Press 'y' to load ALL singleStudyModelFits or 's' to continue and","\n"))
        cat(crayon::blue("select the study/studies you whish to load the fits for. If you wish to quite, press 'q'. Press ENTER afterwards","\n"))
        char <- readline(" ")
        while (!(char == 's') & !(char == 'S') & !(char == 'y') & !(char == 'Y') & !(char == 'q') & !(char == 'Q')) {
          cat((crayon::blue("Please press 'y' to load ALL, 's' to specify the study/studies to load, or 'q' to quit. Press ENTER afterwards.", "\n")))
          char <- readline(" ")
        }
        if (char == 'y' | char == 'Y') {
          for (i in unlist(primaryStudies$studyNumber)) loadSingleStudyModelFit <- c(loadSingleStudyModelFit, i)
        } else if (char == 's' | char == 'S') {
          cat(crayon::blue("Which SingleStudyModelFits would you like to load?", "\n"))
          cat(crayon::blue("Please enter the no. of study/studies separated by commas!", "\n"))
          chars <- readline(" ")
          chars <- gsub(" ", "", chars, fixed = TRUE)
          chars <- unlist(strsplit(chars, ","))
          loadSingleStudyModelFit <- c(loadSingleStudyModelFit, chars)
        } else if (char == 'q' | char == 'Q') {
          ErrorMsg <- "Good luck for the next try!"
          stop(ErrorMsg)
        }
      } else {
        loadSingleStudyModelFit <- c(loadSingleStudyModelFit, unlist(primaryStudies$studyNumber))
      }
    }
  } ### END check specification of primary studies ###


  #######################################################################################################################
  ##################################### Fit ctsem Model to each Primary Study ###########################################
  #######################################################################################################################
  {
    # loop through all primary studies
    studyFit <- studyFitCI <- studyFit_Minus2LogLikelihood <- studyFit_estimatedParameters <- list()
    model_Drift_Coef <- model_Drift_SE <- model_Drift_CI <- list()
    model_Diffusion_Coef <- model_Diffusion_SE <- model_Diffusion_CI <- list()
    model_T0var_Coef <- model_T0var_SE <- model_T0var_CI <- list()
    model_Cint_Coef <- model_Cint_SE <- model_Cint_CI <- list()
    model_popsd <- list() # model_popcov <- model_popcor <- list()
    resultsSummary <- list()
    model_popcov_m <- model_popcov_sd <- model_popcov_T <- model_popcov_025 <- model_popcov_50 <- model_popcov_975 <- list()
    model_popcor_m <- model_popcor_sd <- model_popcor_T <- model_popcor_025 <- model_popcor_50 <- model_popcor_975 <- list()
    estProb <- list()
    hessianWarning <- list()

    studyCompute <- identical(compute$tasks, "studies") && isTRUE(fit)
    studyRequests <- vector("list", n.studies)
    loadedStudy <- logical(n.studies)
    collectStudy <- quote({
      if (isTRUE(fit) && !loadedStudy[i]) {
          gc()
          #
          studyFit[[i]]$empraw <- empraw[[i]]
          studyFit[[i]]$cotima_study <- studyList[[i]]
          if (is.null(studyFit[[i]]$standata$priors)) studyFit[[i]]$standata$priors <- 0 # CHD added Sep 2023
          # The native ctsem Julia dispatcher has already requested the
          # parameter-matrix summary needed by the legacy CoTiMA extraction.
          # Replacing it with summary(fit) here would drop `parmatrices` for
          # this backend.
          if (!identical(backend, "ctsemjulia")) {
            studyFit[[i]]$resultsSummary <- .ctmaWithInitStudyWarnings(
              summary(studyFit[[i]]), i, n.studies, studyList[[i]]
            )
          }

          df <- "deprecated"
          studyFit[[i]]$resultsSummary$'df (CoTiMA)' <- df
          #} # END if (!(studyList[[i]]$originalStudyNo %in% ...

          # SAVE
          if ( (length(saveSingleStudyModelFit) > 1) & (studyList[[i]]$originalStudyNo %in% saveSingleStudyModelFit[-1]) ) {
            x1 <- paste0(saveSingleStudyModelFit[1], " studyFit", studyList[[i]]$originalStudyNo, ".rds"); x1
            x2 <- paste0(saveSingleStudyModelFit[1], " singleStudyFits/"); x2
            tmp <- studyFit[[i]]
            tmp$empraw <- empraw[[i]]
            tmp$primaryStudyList=primaryStudies
            tmp$studyList=studyList
            tmp$argumentListOfCtmaInit=list(primaryStudies=deparse(substitute(primaryStudies)),

                                            activeDirectory=activeDirectory,
                                            chains=chains,
                                            checkSingleStudyResults=checkSingleStudyResults,
                                            cint=cint,
                                            coresToUse=coresToUse,
                                            customPar=customPar,
                                            diff=diff,
                                            digits=digits,
                                            doPar=doPar,
                                            drift=drift,

                                            finishsamples=finishsamples,
                                            indVarying=indVarying,
                                            indVaryingT0=indVaryingT0,
                                            iter=iter,
                                            lambda=lambda,
                                            loadSingleStudyModelFit=loadSingleStudyModelFit,
                                            manifestMeans=manifestMeans,
                                            manifestVars=manifestVars,
                                            n.latent=n.latent,
                                            n.manifest=n.manifest,
                                            optimize=optimize,
                                            primaryStudies=primaryStudies,
                                            priors=priors,
                                            sameInitialTimes=sameInitialTimes,
                                            saveRawData=saveRawData,
                                            saveSingleStudyModelFit=saveSingleStudyModelFit,
                                            scaleTI=scaleTI,
                                            scaleTime=scaleTime,
                                            silentOverwrite=silentOverwrite,
                                            T0means=T0means,
                                            T0var=T0var,
                                            useSV=useSV,
                                            verbose=verbose,
                                            randomIntercepts=randomInterceptsSettings,
                                            CoTiMAStanctArgs=CoTiMAStanctArgs,
                                            backend=backend,
                                            backendArgs=resolved_backend_args)

            ctmaSaveFile( activeDirectory, tmp, x1, x2, silentOverwrite=silentOverwrite)
          }
      }

      if (fit == TRUE) {
        studyFit[[i]]$resultsSummary <- .ctmaInitMaskNativeSummary(studyFit[[i]])
        resultsSummary <- studyFit[[i]]$resultsSummary; resultsSummary
        if (!(exists("currenModel"))) currentModel <- studyFit[[i]]$ctstanmodelbase

        studyFit_Minus2LogLikelihood[[i]] <- -2 * resultsSummary$loglik
        studyFit_estimatedParameters[[i]] <- resultsSummary$npars

        # Subsequent if ... else ... became necessary with ctsem 3.4.1 where the namens of the rows (rownames) were moved into a column named "matrix"
        tmp <- grep("toV", rownames(resultsSummary$popmeans)); tmp
        if (!(length(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "Mean"]) == 0)) {
          model_Drift_Coef[[i]] <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "Mean"], n.latent)); model_Drift_Coef[[i]]
        } else { # new version of ctsem
          model_Drift_Coef[[i]] <- c(matrix(resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "DRIFT", "Mean"], n.latent)); model_Drift_Coef[[i]]
        }
        names(model_Drift_Coef[[i]]) <- driftFullNames; model_Drift_Coef[[i]]

        if (!(length(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "Sd"]) == 0)) {
          model_Drift_SE[[i]] <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "Sd"], n.latent, byrow=FALSE)); model_Drift_SE[[i]]
        } else { # new version of ctsem
          model_Drift_SE[[i]] <- c(matrix(resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "DRIFT", "sd"], n.latent, byrow=FALSE)); model_Drift_SE[[i]]
        }
        names(model_Drift_SE[[i]]) <- c(driftFullNames); model_Drift_SE[[i]]

        if (!(length(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "2.5%"]) == 0 )) {
          tmp1 <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "2.5%"], n.latent, byrow=FALSE)); tmp1
          tmp2 <- c(matrix(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DRIFT", "97.5%"], n.latent, byrow=FALSE)); tmp2
        } else {
          tmp1 <- c(matrix(resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "DRIFT", "2.5%"], n.latent, byrow=FALSE)); tmp1
          tmp2 <- c(matrix(resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "DRIFT", "97.5%"], n.latent, byrow=FALSE)); tmp2
        }

        #tmp1; tmp2
        # CHD 19. Nov. 2023 Check if LL == UL or LL == 0  or UL == 0, indicating estimation problems (estProb)
        tmp3a <- which(tmp1 - tmp2 == 0); tmp3a
        tmp3b <- which(tmp1 == 0); tmp3b
        tmp3c <- which(tmp2 == 0); tmp3c
        tmp4 <- unique(c(tmp3a, tmp3b, tmp3c)); tmp4
        tmp5 <- grep("toV", rownames(resultsSummary$popmeans)); tmp5
        tmp6 <- rownames(resultsSummary$popmeans)[tmp]; tmp6
        estProb <- .ctmaInitRecordIntervalProblem(
          estProb, studyList[[i]]$originalStudyNo, "DRIFT", tmp4, driftFullNames
        )

        #estProb
        model_Drift_CI[[i]] <- c(rbind(tmp1, tmp2)); model_Drift_CI[[i]]
        tmp3 <- c(rbind(paste0(driftFullNames, "LL"),
                        paste0(driftFullNames, "UL"))); tmp3
        names(model_Drift_CI[[i]]) <- tmp3; model_Drift_CI[[i]]

        tmp <- grep("diff", rownames(resultsSummary$popmeans)); tmp
        #(!(length(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "Mean"]) == 0))
        if (!(length(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "Mean"]) == 0)) {
          model_Diffusion_Coef[[i]] <- (resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "Mean"])
          names(model_Diffusion_Coef[[i]]) <- rownames(resultsSummary$popmeans)[tmp]
        } else {
          model_Diffusion_Coef[[i]] <- (resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "DIFFUSIONcov", "Mean"]); model_Diffusion_Coef[[i]]
          # CHD 3.3.2026
          #(randomInterceptsSettings != FALSE)
          if (randomInterceptsSettings != FALSE) {
            tmp9 <- matrix(model_Diffusion_Coef[[i]], length(model_Diffusion_Coef[[i]])^.5, length(model_Diffusion_Coef[[i]])^.5)
            model_Diffusion_Coef[[i]] <- c(tmp9[1:n.latent, 1:n.latent]); model_Diffusion_Coef[[i]]
            names(model_Diffusion_Coef[[i]]) <- c(OpenMx::vech2full(rownames(resultsSummary$popmeans)[grep("diff_", rownames(resultsSummary$popmeans))]))
          } else {
            model_Diffusion_Coef[[i]] <- (resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "DIFFUSIONcov", "Mean"])
            names(model_Diffusion_Coef[[i]]) <- c(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp]))
          }
        }
        #model_Diffusion_Coef[[i]]

        #(!(is.null(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "Sd"])))
        if (!(is.null(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "Sd"]))) {
          model_Diffusion_SE[[i]] <- (resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "sd"]) #; model_Diffusion_SE[[i]]
          names(model_Diffusion_SE[[i]]) <- rownames(resultsSummary$popmeans)[tmp]
        } else {
          model_Diffusion_SE[[i]] <- (resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "DIFFUSIONcov", "sd"]); model_Diffusion_SE[[i]]
          if (randomInterceptsSettings != FALSE) {
            tmp9 <- matrix(model_Diffusion_SE[[i]], length(model_Diffusion_SE[[i]])^.5, length(model_Diffusion_SE[[i]])^.5)
            model_Diffusion_SE[[i]] <- c(tmp9[1:n.latent, 1:n.latent]); model_Diffusion_SE[[i]]
            names(model_Diffusion_SE[[i]]) <- c(OpenMx::vech2full(rownames(resultsSummary$popmeans)[grep("diff_", rownames(resultsSummary$popmeans))]))
          } else {
            model_Diffusion_SE[[i]] <- resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "DIFFUSIONcov", "sd"] #; model_Diffusion_SE[[i]]
            names(model_Diffusion_SE[[i]]) <- c(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp]))
          }
        }
        #model_Diffusion_SE[[i]]

        if (!(length(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "2.5%"])) == 0) {
          tmp1 <- resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "2.5%"]; tmp1
          tmp2 <- resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "DIFFUSIONcov", "97.5%"]; tmp2
          model_Diffusion_CI[[i]] <- c(rbind(tmp1, tmp2)); model_Diffusion_CI[[i]]
          tmp3 <- c(rbind(paste0(rownames(resultsSummary$popmeans)[tmp], "LL"),
                          paste0(rownames(resultsSummary$popmeans)[tmp], "UL"))); tmp3
          names(model_Diffusion_CI[[i]]) <- tmp3; model_Diffusion_CI[[i]]
        } else {
          tmp1 <- resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "DIFFUSIONcov", "2.5%"]; tmp1
          tmp2 <- resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "DIFFUSIONcov", "97.5%"]; tmp2
          model_Diffusion_CI[[i]] <- c(rbind(tmp1, tmp2)); model_Diffusion_CI[[i]]
          tmp3 <- c(rbind(paste0(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp]), "LL"),
                          paste0(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp]), "UL"))); tmp3
          names(model_Diffusion_CI[[i]]) <- tmp3; model_Diffusion_CI[[i]]
        }
        #model_Diffusion_CI[[i]]

        # CHD 19. Nov. 2023 Check if LL == UL or LL == 0  or UL == 0, indicating estimation problems (estProb)
        tmp3a <- which(tmp1 - tmp2 == 0); tmp3a
        tmp3b <- which(tmp1 == 0); tmp3b
        tmp3c <- which(tmp2 == 0); tmp3c
        tmp4 <- unique(c(tmp3a, tmp3b, tmp3c)); tmp4
        tmp5 <- grep("diff", rownames(resultsSummary$popmeans)); tmp5
        tmp6 <- c(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp5])); tmp6
        estProb <- .ctmaInitRecordIntervalProblem(
          estProb, studyList[[i]]$originalStudyNo, "DIFFUSION", tmp4, names(model_Diffusion_Coef[[i]])
        )

        #estProb
        #if ( (randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST") )  target <- "T0cov_" else target <- "0var"
        if ( (randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST") | (randomIntercepts == "CINT") )  target <- "T0cov_" else target <- "0var"
        tmp <- grep(target, rownames(resultsSummary$popmeans)); tmp
        if (length(tmp) == 0) {
          tmp <- grep("0cov", resultsSummary$parmatrices$matrix)
          tmp2 <- (resultsSummary$parmatrices[tmp, c("matrix", "row", "col")]); tmp2
          T0covNames <- c()
          for (m in 1:nrow(tmp2)) T0covNames[m] <- paste0(tmp2[m, 1], tmp2[m, 2], tmp2[m, 3])
        } else {
          T0covNames <- c(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp]))
        }

        if (!(length(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "Mean"])) == 0 ) {
          model_T0var_Coef[[i]] <- (resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "Mean"])
          names(model_T0var_Coef[[i]]) <- rownames(resultsSummary$popmeans)[tmp]; model_T0var_Coef[[i]]
        }  else {
          model_T0var_Coef[[i]] <- (resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "T0cov", "Mean"]); model_T0var_Coef[[i]]
          if (randomInterceptsSettings != FALSE) {
            model_T0var_Coef[[i]] <- (resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "T0cov", "Mean"]);  model_T0var_Coef[[i]]
            model_T0var_Coef[[i]] <- c(matrix(model_T0var_Coef[[i]], n.latent * 2L, n.latent * 2L)[1:n.latent, 1:n.latent]); model_T0var_Coef[[i]]
            names(model_T0var_Coef[[i]]) <- T0covNames
          } else {
            if (length(model_T0var_Coef[[i]]) != n.latent^2) {
              names(model_T0var_Coef[[i]]) <- c(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp]))
            } else {
              names(model_T0var_Coef[[i]]) <- T0covNames
            }
          }
          #if ( (randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST") )  {
          #  model_T0var_Coef[[i]] <- (resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "T0cov", "Mean"]);  model_T0var_Coef[[i]]
          #  model_T0var_Coef[[i]] <- c(matrix(model_T0var_Coef[[i]], n.latent^2, n.latent^2)[1:n.latent, 1:n.latent]); model_T0var_Coef[[i]]
          #  names(model_T0var_Coef[[i]]) <- T0covNames
          #}
        }
        #model_T0var_Coef[[i]]

        if (!(is.null(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "Sd"]))) {
          model_T0var_SE[[i]] <- (resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "Sd"]); model_T0var_SE[[i]]
          names(model_T0var_SE[[i]]) <- rownames(resultsSummary$popmeans)[tmp]; model_T0var_SE[[i]]
        } else {
          model_T0var_SE[[i]] <- (resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "T0cov", "sd"])
          if (length(model_T0var_SE[[i]]) != n.latent^2) {
            names(model_T0var_SE[[i]]) <- c(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp]))
          } else {
            names(model_T0var_SE[[i]]) <- T0covNames
          }
          #if ( (randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST") ) {
          if ( (randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST") | (randomIntercepts == "CINT") ) {
            model_T0var_SE[[i]] <- (resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "T0cov", "sd"])
            model_T0var_SE[[i]] <- c(matrix(model_T0var_SE[[i]], n.latent * 2L, n.latent * 2L)[1:n.latent, 1:n.latent])
            names(model_T0var_SE[[i]]) <- T0covNames
          }
        }
        #model_T0var_SE[[i]]

        #(!(length(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "2.5%"]) == 0))
        if (!(length(resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "2.5%"]) == 0)) {
          tmp1 <- resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "2.5%"]; tmp1
          tmp2 <- resultsSummary$parmatrices[rownames(resultsSummary$parmatrices) == "T0VAR", "97.5%"]; tmp2
          model_T0var_CI[[i]] <- c(rbind(tmp1, tmp2)); model_T0var_CI[[i]]
          tmp3 <- c(rbind(paste0(rownames(resultsSummary$popmeans)[tmp], "LL"),
                          paste0(rownames(resultsSummary$popmeans)[tmp], "UL"))); tmp3
          names(model_T0var_CI[[i]]) <- tmp3; model_T0var_CI[[i]]
        } else {
          tmp1 <- resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "T0cov", "2.5%"]; tmp1
          tmp2 <- resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "T0cov", "97.5%"]; tmp2
          model_T0var_CI[[i]] <- c(rbind(tmp1, tmp2)); model_T0var_CI[[i]]
          if (length(tmp1) != n.latent^2) {
            tmp3 <- c(rbind(paste0(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp]), "LL"),
                            paste0(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp]), "UL"))); tmp3
            names(model_T0var_CI[[i]]) <- tmp3; model_T0var_CI[[i]]
          } else {
            tmp3 <- c(rbind(paste0(T0covNames, "LL"),
                            paste0(T0covNames, "UL"))); tmp3
            names(model_T0var_CI[[i]]) <- tmp3; model_T0var_CI[[i]]
          }

          #if ( (randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST") )  {
          if ( (randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST") | (randomIntercepts == "CINT") )  {
            tmp1 <- resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "T0cov", "2.5%"]; tmp1
            tmp1 <- c(matrix(tmp1, n.latent * 2L, n.latent * 2L)[1:n.latent, 1:n.latent])
            tmp2 <- resultsSummary$parmatrices[resultsSummary$parmatrices[, "matrix"] == "T0cov", "97.5%"]; tmp2
            tmp2 <- c(matrix(tmp2, n.latent * 2L, n.latent * 2L)[1:n.latent, 1:n.latent])
            model_T0var_CI[[i]] <- c(rbind(tmp1, tmp2)); model_T0var_CI[[i]]
            tmp3 <- c(rbind(paste0(T0covNames, "LL"),
                            paste0(T0covNames, "UL"))); tmp3
            names(model_T0var_CI[[i]]) <- tmp3; model_T0var_CI[[i]]
          }
        }
        #model_T0var_CI[[i]]


        # CHD 19. Nov. 2023 Check if LL == UL or LL == 0  or UL == 0, indicating estimation problems (estProb)
        tmp3a <- which(tmp1 - tmp2 == 0); tmp3a
        tmp3b <- which(tmp1 == 0); tmp3b
        tmp3c <- which(tmp2 == 0); tmp3c
        tmp4 <- unique(c(tmp3a, tmp3b, tmp3c)); tmp4
        tmp5 <- grep("0var", rownames(resultsSummary$popmeans)); tmp5
        if (length(tmp5) == 0) tmp5 <- grep("0cov", rownames(resultsSummary$popmeans))
        tmp5
        tmp6 <- c(OpenMx::vech2full(rownames(resultsSummary$popmeans)[tmp5])); tmp6
        estProb <- .ctmaInitRecordIntervalProblem(
          estProb, studyList[[i]]$originalStudyNo, "T0VAR", tmp4, names(model_T0var_Coef[[i]])
        )

        #estProb
        if (indVarying %in% c("MANIFEST", "CINT") &&
            !randomIntercepts %in% c("CINT", "MANIFEST")) {
          if (identical(backend, "ctsemjulia")) {
            native_re <- .ctmaInitCtsemJuliaRandomEffectSummary(
              resultsSummary, digits = digits
            )
            model_popsd[[i]] <- native_re$popsd
            model_popcov_m[[i]] <- native_re$popcov_mean
            model_popcov_sd[[i]] <- native_re$popcov_sd
            model_popcov_T[[i]] <- native_re$popcov_T
            model_popcov_025[[i]] <- native_re$popcov_025
            model_popcov_50[[i]] <- native_re$popcov_50
            model_popcov_975[[i]] <- native_re$popcov_975
            model_popcor_m[[i]] <- native_re$popcor_mean
            model_popcor_sd[[i]] <- native_re$popcor_sd
            model_popcor_T[[i]] <- native_re$popcor_T
            model_popcor_025[[i]] <- native_re$popcor_025
            model_popcor_50[[i]] <- native_re$popcor_50
            model_popcor_975[[i]] <- native_re$popcor_975
          } else {
            e <- ctsem::ctExtract(studyFit[[i]])
            model_popsd_tmp <- resultsSummary$popsd
            if (dim(model_popsd_tmp)[1] != n.latent) {
            model_popsd[[i]] <- resultsSummary$popsd
            model_popcov_m[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, mean), digits = digits)
            model_popcov_sd[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, stats::sd), digits = digits)
            model_popcov_T[[i]] <- round(ctsem::ctCollapse(e$popcov, 1, mean)/ctsem::ctCollapse(e$popcov, 1, stats::sd), digits)
            model_popcov_025[[i]] <- ctsem::ctCollapse(e$popcov, 1, function(x) .ctmaInitDrawQuantile(x, .025))
            model_popcov_50[[i]] <- ctsem::ctCollapse(e$popcov, 1, function(x) .ctmaInitDrawQuantile(x, .50))
            model_popcov_975[[i]] <- ctsem::ctCollapse(e$popcov, 1, function(x) .ctmaInitDrawQuantile(x, .975))
            # convert to correlations and do the same (array to list then list to array)
            e$popcor <- lapply(seq(dim(e$popcov)[1]), function(x) e$popcov[x , ,])
            e$popcor <- lapply(e$popcor, stats::cov2cor)
            e$popcor <- array(unlist(e$popcor), dim=c(n.latent*2, n.latent*2, length(e$popcor)))
            if (any(!is.finite(e$popcor))) warning(paste0("Study ",
              studyList[[i]]$originalStudyNo,
              ": random-effect correlation draws contain undefined values; affected quantiles are reported as NA."), call. = FALSE)
            model_popcor_m[[i]] <- round(ctsem::ctCollapse(e$popcor, 3, mean), digits = digits)
            model_popcor_sd[[i]] <- round(ctsem::ctCollapse(e$popcor, 3, stats::sd), digits = digits)
            model_popcor_T[[i]] <- round(ctsem::ctCollapse(e$popcor, 3, mean)/ctsem::ctCollapse(e$popcor, 3, stats::sd), digits)
            model_popcor_025[[i]] <- ctsem::ctCollapse(e$popcor, 3, function(x) .ctmaInitDrawQuantile(x, .025))
            model_popcor_50[[i]] <- ctsem::ctCollapse(e$popcor, 3, function(x) .ctmaInitDrawQuantile(x, .50))
            model_popcor_975[[i]] <- ctsem::ctCollapse(e$popcor, 3, function(x) .ctmaInitDrawQuantile(x, .975))
            #model_popcor <- stats::cov2cor(model_popcov_m)
            }
          }
        }

        # CHD 22.1.2024
        if ( (indVarying != "MANIFEST") & (indVarying != 'CINT') & (randomIntercepts != "CINT") & (randomIntercepts != "MANIFEST") ) {
          model_popsd[[i]] <- "no random intercepts estimated"
          model_popcov_m[[i]] <- model_popcov_sd[[i]] <- model_popcov_T[[i]] <- model_popcov_025[[i]] <- model_popcov_50[[i]] <- model_popcov_975[[i]] <- "no random intercepts estimated"
          model_popcor_m[[i]] <- model_popcor_sd[[i]] <- model_popcor_T[[i]] <- model_popcor_025[[i]] <- model_popcor_50[[i]] <- model_popcor_975[[i]] <- "no random intercepts estimated"
        }

        if ( (randomIntercepts == "CINT") | (randomIntercepts == "MANIFEST")  ) {
          extraction_fit <- if (inherits(studyFit[[i]], "ctmaJuliaFit")) {
            studyFit[[i]]$ctsem_fit
          } else {
            studyFit[[i]]
          }
          if (is.null(extraction_fit)) {
            stop(
              "The Julia fit has no retained ctsem transformation object for random-intercept extraction.",
              call. = FALSE
            )
          }
          e <- ctsem::ctExtract(extraction_fit)
          ctsem::ctCollapse(e$pop_T0cov, 1, mean)
          model_popcov_m[[i]] <- round(ctsem::ctCollapse(e$pop_T0cov, 1, mean), digits = digits)
          model_popcov_sd[[i]] <- round(ctsem::ctCollapse(e$pop_T0cov, 1, stats::sd), digits = digits)
          model_popcov_T[[i]] <- round(ctsem::ctCollapse(e$pop_T0cov, 1, mean)/ctsem::ctCollapse(e$pop_T0cov, 1, stats::sd), digits)
          model_popcov_025[[i]] <- ctsem::ctCollapse(e$pop_T0cov, 1, function(x) .ctmaInitDrawQuantile(x, .025))
          model_popcov_50[[i]] <- ctsem::ctCollapse(e$pop_T0cov, 1, function(x) .ctmaInitDrawQuantile(x, .50))
          model_popcov_975[[i]] <- ctsem::ctCollapse(e$pop_T0cov, 1, function(x) .ctmaInitDrawQuantile(x, .975))
          # convert to correlations and do the same (array to list then list to array)
          e$popcor <- lapply(seq(dim(e$pop_T0cov)[1]), function(x) e$pop_T0cov[x , ,])
          e$popcor <- lapply(e$popcor, stats::cov2cor)
          e$popcor <- array(unlist(e$popcor), dim=c(n.latent*2, n.latent*2, length(e$popcor)))
          if (any(!is.finite(e$popcor))) warning(paste0("Study ",
            studyList[[i]]$originalStudyNo,
            ": random-intercept correlation draws contain undefined values; affected quantiles are reported as NA."), call. = FALSE)
          model_popcor_m[[i]] <- round(ctsem::ctCollapse(e$popcor, 3, mean), digits = digits)
          model_popcor_sd[[i]] <- round(ctsem::ctCollapse(e$popcor, 3, stats::sd), digits = digits)
          model_popcor_T[[i]] <- round(ctsem::ctCollapse(e$popcor, 3, mean)/ctsem::ctCollapse(e$popcor, 3, stats::sd), digits)
          model_popcor_025[[i]] <- ctsem::ctCollapse(e$popcor, 3, function(x) .ctmaInitDrawQuantile(x, .025))
          model_popcor_50[[i]] <- ctsem::ctCollapse(e$popcor, 3, function(x) .ctmaInitDrawQuantile(x, .50))
          model_popcor_975[[i]] <- ctsem::ctCollapse(e$popcor, 3, function(x) .ctmaInitDrawQuantile(x, .975))
        } # end     if ( (randomIntercepts == TRUE) | (randomIntercepts == "MANIFEST")  )

      } # end if (fit == TRUE)

    })

    for (i in seq_len(n.studies)) {
      #i <- 12
      notLoadable <- TRUE
      if ( (length(loadSingleStudyModelFit) > 1) & (studyList[[i]]$originalStudyNo %in% loadSingleStudyModelFit[-1]) ) {
        tmp1 <- paste0(" LOADING SingleStudyFit ", i, " of ", n.studies, " (Study: ", studyList[[i]]$originalStudyNo, ") ")
        tmp2 <- nchar(tmp1); tmp2
        tmp3 <- (81 - tmp2)/2; tmp3
        tmp4 <- strrep("#", round(tmp3 + 0.45, 0)); tmp4
        tmp5 <- strrep("#", round(tmp3 - 0.45, 0)); tmp5
        tmp6 <- paste0(tmp4, tmp1, tmp5); tmp6
        Msg <- paste0("################################################################################# \n" , tmp6 ,"\n#################################################################################")
        message(Msg)
        x1 <- paste0(activeDirectory, loadSingleStudyModelFit[1], " singleStudyFits/",loadSingleStudyModelFit[1], " studyFit", studyList[[i]]$originalStudyNo, ".rds"); x1
        #file.exists(x1)
        if (file.exists(x1)) {
          notLoadable <- FALSE
          loadedStudy[i] <- TRUE
          studyFit[[i]] <- readRDS(file=x1)
          # CHD added 21 Sep 2022
          empraw[[i]] <- studyFit[[i]]$empraw
        } else {
          notLoadable <- TRUE
        }
      } # END if ( (length(loadSingleStudyModelFit) > 1) &  ...

      tmpLogic <- 0
      if ((studyList[[i]]$originalStudyNo %in% loadSingleStudyModelFit[-1]) & (notLoadable == TRUE) ) tmpLogic <- 1
      if (!(studyList[[i]]$originalStudyNo %in% loadSingleStudyModelFit[-1]) ) tmpLogic <- 1
      if (tmpLogic == 1) {
        tmp1 <- paste0(" Fitting SingleStudyModel ", i, " of ", n.studies, " (Study: ", studyList[[i]]$originalStudyNo, ") ")
        tmp2 <- nchar(tmp1); tmp2
        tmp3 <- (81 - tmp2)/2; tmp3
        tmp4 <- strrep("#", round(tmp3 + 0.45, 0)); tmp4
        tmp5 <- strrep("#", round(tmp3 - 0.45, 0)); tmp5
        tmp6 <- paste0(tmp4, tmp1, tmp5); tmp6
        Msg <- paste0("################################################################################# \n", tmp6, "\n#################################################################################")
        message(Msg)

        if (!(optimize)) {
          customPar <- FALSE


          tmp1a <- paste0(" Bayesian sampling was selected, which does require appropriate scaling of time. ")
          tmp2 <- nchar(tmp1a); tmp2
          tmp3 <- (81 - tmp2)/2; tmp3
          tmp4 <- strrep("#", round(tmp3 + 0.45, 0)); tmp4
          tmp5 <- strrep("#", round(tmp3 - 0.45, 0)); tmp5
          tmp6a <- paste0(tmp4, tmp1a, tmp5); tmp6

          tmp1b <- paste0(" See the end of the summary output ")
          tmp2 <- nchar(tmp1b); tmp2
          tmp3 <- (81 - tmp2)/2; tmp3
          tmp4 <- strrep("#", round(tmp3 + 0.45, 0)); tmp4
          tmp5 <- strrep("#", round(tmp3 - 0.45, 0)); tmp5
          tmp6b <- paste0(tmp4, tmp1b, tmp5); tmp6

          Msg <- paste0("################################################################################# \n", tmp6a, "\n", tmp6b, "\n#################################################################################")
          message(Msg)
        }

        # select correct template
        currentTpoints <- (lapply(studyList, function(extract) extract$timePoints))[[i]]; currentTpoints
        driftParamsTmp <- driftParams
        diffParamsTmp <- diffParams
        longestLag <- max((lapply(studyList, function(extract) extract$delta_t))[[i]]); longestLag
        if ((longestLag > 6)  & (customPar)) {
          counter <- 0
          for (h in 1:(n.latent)) {
            for (j in 1:(n.latent)) {
              counter <- counter + 1
              if (h == j) {
                driftParamsTmp[counter] <- paste0(driftParams[counter], paste0("|-log1p_exp(-param *.1 -2)"))
                diffParamsTmp[counter] <- paste0(diffParams[counter], paste0("|log1p_exp(-param *.1 -2)"))
              }
            }
          }
        }

        if ((indVarying == 'cint') | (indVarying == 'Cint')) indVarying <- 'CINT'
        if ( (randomIntercepts != "CINT") & (randomIntercepts != "MANIFEST") ) {

          # CHD 27.11.2024
          if (indVarying == 'CINT')  {
            print(paste0("#################################################################################"))
            print(paste0("######## Note: Individually varying intercepts model requested.  #########"))
            print(paste0("#################################################################################"))

            # CHD 2.12.2024
            #T0meansParams <- 0
            T0meansParams <- 'auto'


            print(paste0("#################################################################################"))
            print(paste0("####################### CT intercepts are set free.  ########################"))
            print(paste0("#################################################################################"))

            CINTParams <- c()
            for (c in 1:n.latent) {
              CINTParams <- c(CINTParams, paste0("cintV", c))
            }
          }

          # CHD 27.11.2024
          if (indVarying == "MANIFEST") {
            print(paste0("#################################################################################"))
            print(paste0("###### Note: Individually varying manifest means model requested.  #######"))
            print(paste0("#################################################################################"))

            # CHD 2.12.2024
            T0meansParams <- 'auto'

            print(paste0("#################################################################################"))
            print(paste0("######### Manifest intercepts (MANIFESTMEANS, not CINT) are set free.  #########"))
            print(paste0("#################################################################################"))

            manifestMeansParams <- 'auto'
          }
        } # end  if ( (randomIntercepts != TRUE) & (randomIntercepts != "MANIFEST") )

        if (utils::packageDescription("ctsem")$Version > "3.10.2") type <- "ct" else type <- "stanct"

    if (!is.null(namedCint)) CINTParams <- namedCint
    if (!is.null(namedManifest)) manifestMeansParams <- namedManifest
    if (!is.null(namedT0)) T0meansParams <- namedT0

        currentModel <- suppressMessages(
          ctsem::ctModel(n.latent=n.latent, n.manifest=n.var, Tpoints=currentTpoints, manifestNames=manifestNames,    # 2 waves in the template only
                         DIFFUSION=matrix(diffParamsTmp, nrow=n.latent, ncol=n.latent), #, byrow=TRUE),
                         DRIFT=matrix(driftParamsTmp, nrow=n.latent, ncol=n.latent),
                         LAMBDA=lambdaParams,
                         T0VAR=T0VARParams,
                         type=type,
                         CINT=matrix(CINTParams, nrow=n.latent, ncol=1),
                         T0MEANS = matrix(c(T0meansParams), nrow = n.latent, ncol = 1),
                         MANIFESTMEANS = matrix(manifestMeansParams, nrow = n.var, ncol = 1),
                         MANIFESTVAR=matrix(manifestVarsParams, nrow=n.var, ncol=n.var)
          )
        )

        if (indVarying == FALSE) currentModel$pars[, "indvarying"] <- FALSE
        # CHD 27.11.2024
        if ( (indVarying != "MANIFEST") & (indVarying != "CINT") ) {
          currentModel$pars[currentModel$pars$matrix %in% 'T0MEANS','indvarying'] <- FALSE
        }

        if (indVarying == 'CINT') {
          currentModel$pars[currentModel$pars$matrix %in% 'CINT','indvarying'] <- TRUE
        } else {
          currentModel$pars[currentModel$pars$matrix %in% 'CINT','indvarying'] <- FALSE
        }

        if (indVarying == "MANIFEST") {
          currentModel$pars[currentModel$pars$matrix %in% 'MANIFESTMEANS','indvarying'] <- TRUE
        } else {
          currentModel$pars[currentModel$pars$matrix %in% 'MANIFESTMEANS','indvarying'] <- FALSE
        }

        currentModel$manifesttype <- binaries

        if ((randomIntercepts == "CINT") | (randomIntercepts == "MANIFEST") ) {
          print(paste0("#################################################################################"))
          print(paste0("# Note: Random intercepts modelled as processes, not as cint or manifest means. #"))
          print(paste0("#################################################################################"))

          nullMat <- matrix(0, n.latent, n.latent); nullMat
          DIFFUSIONtmp <- matrix(diffParamsTmp, nrow=n.latent, ncol=n.latent); DIFFUSIONtmp
          DIFFUSIONtmp <- rbind(cbind(DIFFUSIONtmp, nullMat),
                                cbind(nullMat, nullMat)); DIFFUSIONtmp
          DRIFTtmp <- matrix(driftParamsTmp, nrow=n.latent, ncol=n.latent); DRIFTtmp
          DRIFTtmp <- rbind(cbind(DRIFTtmp, diag(n.latent)),
                            cbind(nullMat, nullMat)); DRIFTtmp
          if (randomIntercepts == "MANIFEST") {
            DRIFTtmp <- matrix(driftParamsTmp, nrow=n.latent, ncol=n.latent); DRIFTtmp
            DRIFTtmp <- rbind(cbind(DRIFTtmp, nullMat),
                              cbind(nullMat, nullMat)); DRIFTtmp
          }
          T0VARtmp = gsub("diff", "T0cov", matrix(diffParamsTmp, n.latent, n.latent)); T0VARtmp
          cint_cint_cov <- gsub("eta", "cint", T0VARtmp); cint_cint_cov
          cint_cint_cov <- gsub("T0c", "C", cint_cint_cov); cint_cint_cov
          T0eta_cint_cov <- matrix(NA, n.latent, n.latent); T0eta_cint_cov
          for (ii in 1:n.latent) {
            for (ll in 1:n.latent) {
              T0eta_cint_cov[ii, ll] <- paste0("Cov_T0eta", ll, "_cint", ii)
            }
          }

          T0VARtmp <- rbind(cbind(T0VARtmp, nullMat),
                            cbind(T0eta_cint_cov, cint_cint_cov)); T0VARtmp
          LAMBDAtmp <- cbind(lambdaParams, nullMat); LAMBDAtmp
          if (randomIntercepts == "MANIFEST") {
            LAMBDAtmp <- cbind(lambdaParams, lambdaParams); LAMBDAtmp
          }
          manifestNamesTmp <- manifestNames; manifestNamesTmp
          latentNamesTmp <- c(latentNames, paste0(latentNames, "_cint")); latentNamesTmp
          T0MEANStmp <- matrix(paste0("Mean", latentNamesTmp), n.latent*2, 1); T0MEANStmp
    if (!is.null(namedT0)) T0MEANStmp[seq_len(n.latent)] <- namedT0
    namedRandom <- if (randomIntercepts == "MANIFEST") namedManifest else namedCint
    if (!is.null(namedRandom)) {
      selected <- .ctma_select_phantom_intercepts(
        DRIFTtmp, LAMBDAtmp, T0VARtmp, T0MEANStmp,
        n.latent, randomIntercepts, namedRandom)
      DRIFTtmp <- selected$drift
      LAMBDAtmp <- selected$lambda
      T0VARtmp <- selected$covariance
      T0MEANStmp <- selected$means
    }

          if (utils::packageDescription("ctsem")$Version > "3.10.2") type <- "ct" else type <- "stanct"

          currentModel <- (
            ctsem::ctModel(n.latent=n.latent*2,
                           n.manifest=n.var,
                           manifestNames=manifestNamesTmp,
                           latentNames = latentNamesTmp,
                           DIFFUSION=DIFFUSIONtmp,
                           DRIFT=DRIFTtmp,
                           LAMBDA=LAMBDAtmp,
                           CINT=matrix(0, nrow=n.latent*2, ncol=1),
                           T0MEANS = T0MEANStmp,
                           MANIFESTMEANS = matrix(manifestMeansParams, nrow=n.latent, ncol=1),
                           MANIFESTVAR=matrix(manifestVarsParams, nrow=n.var, ncol=n.var),
                           T0VAR = T0VARtmp,
                           type = type)
          )
          currentModel$pars$indvarying <- FALSE

          if (!(optimize)) {
            customPar <- FALSE

            tmp1a <- paste0(" Bayesian sampling was selected, which does require appropriate scaling of time. ")
            tmp2 <- nchar(tmp1a); tmp2
            tmp3 <- (81 - tmp2)/2; tmp3
            tmp4 <- strrep("#", round(tmp3 + 0.45, 0)); tmp4
            tmp5 <- strrep("#", round(tmp3 - 0.45, 0)); tmp5
            tmp6a <- paste0(tmp4, tmp1a, tmp5); tmp6a

            tmp1b <- paste0(" See the end of the summary output ")
            tmp2 <- nchar(tmp1b); tmp2
            tmp3 <- (81 - tmp2)/2; tmp3
            tmp4 <- strrep("#", round(tmp3 + 0.45, 0)); tmp4
            tmp5 <- strrep("#", round(tmp3 - 0.45, 0)); tmp5
            tmp6b <- paste0(tmp4, tmp1b, tmp5); tmp6b

            Msg <- paste0("################################################################################# \n", tmp6a, "\n", tmp6b, "\n#################################################################################")
            message(Msg)
            #Msg <- "Bayesian sampling was selected, which does require appropriate scaling of time. See the end of the summary output\n"
            #message(Msg)
          }

        }

        currentModel <- .ctma_restore_manifesttype(currentModel, binaries)
        currentModel <- .ctmaApplyT0MeanVariation(currentModel, indVaryingT0)
        currentModel <- .ctmaDeterminedDiffusionModel(currentModel, determinedDiffusions, stationaryVariance = stationaryVariance)

        saved_models <- attr(primaryStudies, "ctma_models")
        if (!is.null(saved_models)) {
          saved <- saved_models[[i]]
          overrides <- attr(primaryStudies, "ctma_model_overrides")
          columns <- character()
          if ("indVarying" %in% overrides) columns <- c(columns, "indvarying")
          if ("indVaryingT0" %in% overrides) columns <- c(columns, "indvarying")
          families <- c(manifestMeans="MANIFESTMEANS", manifestVars="MANIFESTVAR",
            cint="CINT", T0means="T0MEANS", T0var="T0VAR", drift="DRIFT",
            diff="DIFFUSION", lambda="LAMBDA")
          changed <- unname(families[intersect(names(families), overrides)])
          key <- function(p) paste(p$matrix,p$row,p$col,sep="/")
          match_rows <- match(key(saved$pars), key(currentModel$pars))
          if (length(columns) || length(changed)) {
            if (anyNA(match_rows)) stop("The requested override changes the model dimensions.",call.=FALSE)
            for (column in unique(columns)) saved$pars[[column]] <- currentModel$pars[[column]][match_rows]
            rows <- which(saved$pars$matrix %in% changed)
            for (column in intersect(c("param","value","transform"),names(saved$pars)))
              saved$pars[rows,column] <- currentModel$pars[match_rows[rows],column]
          }
          currentModel <- saved
        }

        wave_note <- .ctmaInitWaveDiagnostic(emprawLong[[i]], currentModel, randomIntercepts)
        if (!is.null(wave_note)) {
          warning(paste0(.ctmaInitStudyWarningLabel(i, n.studies, studyList[[i]]),
            .ctmaInitWaveMessage(wave_note, first = !.ctma_wave_warning_seen)),
            call. = FALSE)
          .ctma_wave_warning_seen <- TRUE
        }

        locationFixed <- currentModel$pars$matrix %in% c("CINT", "MANIFESTMEANS", "T0MEANS") &
    is.na(currentModel$pars$param)
  currentModel$pars$indvarying[locationFixed] <- FALSE
    .ctma_check_binary_intercepts(
          currentModel, binaries,
          study_label = .ctmaInitStudyWarningLabel(i, n.studies, studyList[[i]]),
          process_intercepts = randomIntercepts %in% c("CINT", "MANIFEST") ||
            indVarying %in% c("CINT", "MANIFEST")
        )

        if (fit == FALSE) {
          print(paste0("#################################################################################"))
          print(paste0("#############  No model is fitted, only data and code are generated. ############"))
          print(paste0("#################################################################################"))
        }

        # FIT STANCT MODEL
        if (fit == TRUE) {
          CoTiMAStanctArgs$optimcontrol$bootstrapUncertainty <- NULL
          if (doPar < 2) {
            # CHD changed 7 Oct 2022
            if (any(is.na(studyList[[i]]$startValues))) inits <- NULL else inits <- studyList[[i]]$startValues
            if (determinedDiffusions) {
              if (!is.null(inits)) .ctma_print_message("Stored primary-study starts use old raw coordinates; new starts will be generated.")
              inits <- NULL
            }

            hessianWarning[[i]] <- FALSE

            if (identical(backend, "julia")) {
              results <- .ctmaWithInitStudyWarnings(.ctmaInitJuliaStudyFit(
                datalong = emprawLong[[i]],
                model = currentModel,
                same_initial_times = sameInitialTimes,
                stan_args = CoTiMAStanctArgs,
                backend_args = resolved_backend_args,
                estimation = requested_estimation,
                start_values = inits,
                digits = digits
              ), i, n.studies, studyList[[i]])
              results_summary <- results$resultsSummary
            } else if (identical(backend, "ctsemjulia")) {
              if (studyCompute) {
                studyRequests[[i]] <- list(
                datalong = emprawLong[[i]],
                model = currentModel,
                same_initial_times = sameInitialTimes,
                backend_args = resolved_backend_args,
                start_values = inits,
                priors = priors,
                finishsamples = finishsamples,
                cores = coresToUse,
                verbose = verbose,
                digits = digits,
                estimation = requested_estimation
                )
                next
              }
              results <- .ctmaWithInitStudyWarnings(.ctmaInitCtsemJuliaStudyFit(
                datalong = emprawLong[[i]],
                model = currentModel,
                same_initial_times = sameInitialTimes,
                backend_args = resolved_backend_args,
                start_values = inits,
                priors = priors,
                finishsamples = finishsamples,
                cores = coresToUse,
                verbose = verbose,
                digits = digits,
                estimation = requested_estimation
              ), i, n.studies, studyList[[i]])
              results_summary <- results$resultsSummary
            } else {
            #results <- suppressMessages(ctsem::ctStanFit(
            results <- .ctmaWithInitStudyWarnings(.ctmaCtFitCompat(
              datalong = emprawLong[[i]],
              model = currentModel,
              fit=fit,
              sameInitialTimes=sameInitialTimes,
              #inits=studyList[[i]]$startValues,
              stanmodeltext = CoTiMAStanctArgs$stanmodeltext,
              iter = CoTiMAStanctArgs$iter,
              intoverstates = CoTiMAStanctArgs$intoverstates,
              binomial = CoTiMAStanctArgs$binomial,
              intoverpop = CoTiMAStanctArgs$intoverpop,
              stationary = CoTiMAStanctArgs$stationary,
              plot = CoTiMAStanctArgs$plot,
              derrind = CoTiMAStanctArgs$derrind,
              optimize = optimize,
              optimcontrol = CoTiMAStanctArgs$optimcontrol,
              nlcontrol = CoTiMAStanctArgs$nlcontrol,
              #nopriors = NA,
              priors = priors,
              forcerecompile = CoTiMAStanctArgs$forcerecompile,
              saveCompile = CoTiMAStanctArgs$saveCompile,
              savescores = CoTiMAStanctArgs$savescores,
              savesubjectmatrices = CoTiMAStanctArgs$savesubjectmatrices,
              saveComplexPars = CoTiMAStanctArgs$saveComplexPars,
              gendata = CoTiMAStanctArgs$gendata,
              control = CoTiMAStanctArgs$control,
              vb = CoTiMAStanctArgs$vb,
              #warmup=CoTiMAStanctArgs$warmup,
              verbose=verbose,
              cores=coresToUse), i, n.studies, studyList[[i]])

            #hessianWarning[[i]] <- list(warn_hessinv = results$warn_hessinv,
            #                            warnings = results$warnings,
            #                            error= results$error)

            #if (is.na(results$error)) { # CHD 20.7.2026
            #if (is.null(results$error)) {
              #results <- results$fit # to match former fitting results w/o error handling # CHD 20.7.2026
              results_summary <- .ctmaWithInitStudyWarnings(
                summary(results, digits=2*digits, parmatrices=TRUE, residualcov=FALSE),
                i, n.studies, studyList[[i]]
              )
            }
            #if (!is.null(results$error)) { # CHD 20.7.2026
              #fit <- FALSE
            #  print(paste0("#################################################################################"))
            #  print(paste0("###########  Model could not be fitted, only data and code are returned #########"))
            #  print(paste0("#################################################################################"))
            #  hessianWarning[[i]] <- list(warn_hessinv = "There was fatal fitting error - no hessian computed.",
            #                              warnings = "There was fatal fitting error - no hessian computed.",
            #                              error= results$error)
            #}
          } # end if forPar


          if (doPar > 1) {
            # parallel re-fitting of problem study

            tmp1 <- paste0(" Parallel fit attmepts requested. Screen remains silent for a while. ")
            tmp2 <- nchar(tmp1); tmp2
            tmp3 <- (81 - tmp2)/2; tmp3
            tmp4 <- strrep("#", round(tmp3 + 0.45, 0)); tmp4
            tmp5 <- strrep("#", round(tmp3 - 0.45, 0)); tmp5
            tmp6 <- paste0(tmp4, tmp1, tmp5); tmp6

            Msg <- paste0("################################################################################# \n", tmp6, "\n#################################################################################")
            message(Msg)

            #Msg <- "Parallel fit attmepts requested. Screen remains silent for a while.\n"
            #message(Msg)

            allfits <- foreach::foreach(p=1:doPar) %dopar% {

              hessianWarning[[i]] <- FALSE

              # CHD changed 7 Oct 2022
              if (any(is.na(studyList[[i]]$startValues))) inits <- NULL else inits <- studyList[[i]]$startValues
            if (determinedDiffusions) {
              if (!is.null(inits)) .ctma_print_message("Stored primary-study starts use old raw coordinates; new starts will be generated.")
              inits <- NULL
            }

              fits <- suppressMessages(run_ctStanFit_logged(.ctmaCtFitCompat(
                datalong = emprawLong[[i]],
                model = currentModel,
                sameInitialTimes=sameInitialTimes,
                inits=inits,
                savesubjectmatrices=CoTiMAStanctArgs$savesubjectmatrices,
                stanmodeltext=CoTiMAStanctArgs$stanmodeltext,
                iter=CoTiMAStanctArgs$iter,
                intoverstates=CoTiMAStanctArgs$intoverstates,
                binomial=CoTiMAStanctArgs$binomial,
                fit=CoTiMAStanctArgs$fit,
                intoverpop=CoTiMAStanctArgs$intoverpop,
                stationary=CoTiMAStanctArgs$stationary,
                plot=CoTiMAStanctArgs$plot,
                derrind=CoTiMAStanctArgs$derrind,
                optimize=CoTiMAStanctArgs$optimize,
                optimcontrol=CoTiMAStanctArgs$optimcontrol,
                nlcontrol=CoTiMAStanctArgs$nlcontrol,
                priors=CoTiMAStanctArgs$priors,
                chains=CoTiMAStanctArgs$chains,
                forcerecompile=CoTiMAStanctArgs$forcerecompile,
                savescores=CoTiMAStanctArgs$savescores,
                gendata=CoTiMAStanctArgs$gendata,
                control=CoTiMAStanctArgs$control,
                verbose=verbose,
                warmup=.ctmaStanWarmup(CoTiMAStanctArgs),
                cores=1) ))

              hessianWarning[[i]] <- list(warn_hessinv = fits$warn_hessinv,
                                          warnings = fits$warnings,
                                          error= fits$error)

              if (is.na(fits$error)) {
                fits <- fits$fit # to match former fitting results w/o error handling
                #results_summary <- summary(results, digits=2*digits, parmatrices=TRUE, residualcov=FALSE)
              } else {
                #fit <- FALSE
                print(paste0("#################################################################################"))
                print(paste0("###########  Model could not be fitted, only data and code are returned #########"))
                print(paste0("#################################################################################"))
                hessianWarning[[i]] <- list(warn_hessinv = "There was fatal fitting error - no hessian computed.",
                                            warnings = "There was fatal fitting error - no hessian computed.",
                                            error= results$error)
              }

              return(fits)
            }
            all_loglik <- unlist(lapply(allfits, function(x) x$stanfit$optimfit$value)); all_loglik
            bestFit <- which(abs(all_loglik) == max(abs(all_loglik)))[1]; bestFit
            #
            results <- allfits[[bestFit]]
          }

          if (!studyCompute) studyFit[[i]] <- results
        }
      }
      if (!studyCompute) eval(collectStudy)
    }

    if (studyCompute) {
      fitted <- .ctmaRunStudyCompute(studyRequests, compute, studyList)
      for (i in which(!vapply(studyRequests, is.null, logical(1))))
        studyFit[[i]] <- fitted[[i]]
    }

    if (studyCompute) {
      for (i in seq_len(n.studies)) eval(collectStudy)
    }


  }

  #######################################################################################################################
  ################################ Combine summary information and fit statistics #######################################
  #######################################################################################################################

  if (fit == TRUE) {
    allStudies_Minus2LogLikelihood <- sum(unlist(studyFit_Minus2LogLikelihood)); allStudies_Minus2LogLikelihood
    allStudies_estimatedParameters <- sum(unlist(studyFit_estimatedParameters)); allStudies_estimatedParameters
    allStudies_df <- "deprecated"
    #
    tmp1 <- driftFullNames; tmp1
    tmp2 <- rep("SE", length(tmp1)); tmp2
    targetNames1 <- c(rbind(tmp1, tmp2)); targetNames1
    #
    targetNames2 <- c()
    for (i in 1:n.latent) {
      for (j in 1:n.latent) {
        if (i != j) {
          targetNames2 <- c(targetNames2, paste0("_eta", i, "_eta", j))
        } else {
          targetNames2 <- c(targetNames2, paste0("_eta", i))
        }
      }
    }
    targetNames2  <- paste0("diff", targetNames2); targetNames2
    #targetNames2  <- c(rbind(targetNames2, rep("SE", 4))); targetNames2
    targetNames2  <- c(rbind(targetNames2, rep("SE", n.latent^2))); targetNames2
    #
    targetNames3 <- targetNames2
    targetNames3 <- gsub("eta", "V", targetNames3)
    targetNames3 <- gsub("diff", "T0var", targetNames3); targetNames3

    allStudiesDRIFT_effects <- matrix(t(cbind(unlist(model_Drift_Coef), unlist(model_Drift_SE)) ), n.studies, 2*n.latent^2, byrow=T)
    colnames(allStudiesDRIFT_effects) <- targetNames1; allStudiesDRIFT_effects

    allStudiesDIFF_effects <- matrix(t(cbind(unlist(model_Diffusion_Coef), unlist(model_Diffusion_SE)) ), n.studies, 2*n.latent^2, byrow=T)
    colnames(allStudiesDIFF_effects) <- targetNames2; allStudiesDIFF_effects

    allStudiesT0VAR_effects <- matrix(t(cbind(unlist(model_T0var_Coef), unlist(model_T0var_SE)) ), n.studies, 2*n.latent^2, byrow=T)
    colnames(allStudiesT0VAR_effects) <- targetNames3; allStudiesT0VAR_effects

    #
    if (!(is.null(scaleTime))) {
      allStudiesDRIFT_effects_original_time_scale <- matrix(t(cbind(unlist(model_Drift_Coef) * scaleTime,
                                                                    unlist(model_Drift_SE) * scaleTime) ), n.studies, 2*n.latent^2, byrow=T)
      colnames(allStudiesDRIFT_effects_original_time_scale) <- targetNames1
      allStudiesDRIFT_effects_original_time_scale <- round(allStudiesDRIFT_effects_original_time_scale, digits)
      #
      allStudiesDIFF_effects_original_time_scale <- matrix(t(cbind(unlist(model_Diffusion_Coef) * scaleTime,
                                                                   unlist(model_Diffusion_SE) * scaleTime) ), n.studies, 2*n.latent^2, byrow=T)
      colnames(allStudiesDIFF_effects_original_time_scale) <- targetNames2
      allStudiesDIFF_effects_original_time_scale <- round(allStudiesDIFF_effects_original_time_scale, digits)
    } else {
      allStudiesDRIFT_effects_original_time_scale <- NULL
      allStudiesDIFFUSION_effects_original_time_scale <- NULL
    }

    source <- lapply(primaryStudies$source, function(extract) paste(extract, collapse=", ")); source
    if (length(source) > n.studies) source[n.studies +1] <- NULL # new 6. July< 2022
    for (l in 1:length(source)) if ( source[[l]] == "NA") source[[l]] <- "Reference not provided"
    #
    allStudiesDRIFT_effects_ext <- cbind(unlist(source), allStudiesDRIFT_effects)
    tmp <- allStudiesDRIFT_effects_ext
    tmp[, 2:(ncol(tmp))] <- round(as.numeric(tmp[, 2:(ncol(tmp))]), digits)
    allStudiesDRIFT_effects_ext <- tmp; allStudiesDRIFT_effects_ext
    #
    allStudiesDIFF_effects_ext <- cbind(unlist(source), allStudiesDIFF_effects)
    tmp <- allStudiesDIFF_effects_ext
    tmp[, 2:(ncol(tmp))] <- round(as.numeric(tmp[, 2:(ncol(tmp))]), digits)
    allStudiesDIFF_effects_ext <- tmp; allStudiesDIFF_effects_ext
    #
    allStudiesT0VAR_effects_ext <- cbind(unlist(source), allStudiesT0VAR_effects)
    tmp <- allStudiesT0VAR_effects_ext
    tmp[, 2:(ncol(tmp))] <- round(as.numeric(tmp[, 2:(ncol(tmp))]), digits)
    allStudiesT0VAR_effects_ext <- tmp; allStudiesT0VAR_effects_ext

    if (!(is.null(allStudiesDRIFT_effects_original_time_scale))) {
      #
      allStudiesDRIFT_effects_original_time_scale_ext <- cbind(unlist(source), allStudiesDRIFT_effects_original_time_scale)
      tmp <- allStudiesDRIFT_effects_original_time_scale_ext
      tmp[, 2:(ncol(tmp))] <- round(as.numeric(tmp[, 2:(ncol(tmp))]), digits)
      allStudiesDRIFT_effects_original_time_scale_ext <- tmp; allStudiesDRIFT_effects_original_time_scale_ext
      #allStudiesDRIFT_effects_ext
      #
      allStudiesDIFF_effects_original_time_scale_ext <- cbind(unlist(source), allStudiesDIFF_effects_original_time_scale)
      tmp <- allStudiesDIFF_effects_original_time_scale_ext
      tmp[, 2:(ncol(tmp))] <- round(as.numeric(tmp[, 2:(ncol(tmp))]), digits)
      allStudiesDIFF_effects_original_time_scale_ext <- tmp; allStudiesDIFF_effects_original_time_scale_ext
      #allStudiesDIFF_effects_ext
    } else {
      allStudiesDRIFT_effects_original_time_scale_ext <- allStudiesDRIFT_effects_ext # new 8.7.2022
      allStudiesDIFF_effects_original_time_scale_ext <- allStudiesDIFF_effects_ext # new 8.7.2022
    }

    # confidence intervals
    allStudiesDriftCI <- matrix(unlist(model_Drift_CI), nrow=n.studies, byrow=TRUE)
    colnames(allStudiesDriftCI) <- names(model_Drift_CI[[1]]); allStudiesDriftCI
    allStudiesDiffusionCI <- matrix(unlist(model_Diffusion_CI), nrow=n.studies, byrow=TRUE)
    colnames(allStudiesDiffusionCI) <- names(model_Diffusion_CI[[1]]); allStudiesDiffusionCI
    allStudiesT0varCI <- matrix(unlist(model_T0var_CI), nrow=n.studies, byrow=TRUE)
    colnames(allStudiesT0varCI) <- names(model_T0var_CI[[1]]); allStudiesT0varCI

    if (!(is.null(scaleTime))) {
      #
      allStudiesDriftCI_original_time_scale <- matrix(round(unlist(model_Drift_CI) * scaleTime, digits),
                                                      nrow=n.studies, byrow=TRUE)
      colnames(allStudiesDriftCI_original_time_scale) <- names(model_Drift_CI[[1]]); allStudiesDriftCI_original_time_scale
      #
      allStudiesDiffCI_original_time_scale <- matrix(round(unlist(model_Diffusion_CI) * scaleTime, digits),
                                                     nrow=n.studies, byrow=TRUE)
      colnames(allStudiesDiffCI_original_time_scale) <- names(model_Diffusion_CI[[1]]); allStudiesDiffCI_original_time_scale
    } else {
      allStudiesDriftCI_original_time_scale <- allStudiesDriftCI  # new 8.7.2022
      allStudiesDiffCI_original_time_scale <- allStudiesDiffusionCI  # new 8.7.2022
    }


    allStudiesCI <- cbind(allStudiesDriftCI, allStudiesDiffusionCI, allStudiesT0varCI); allStudiesCI
    allStudiesCI_ext <- cbind(allStudiesDRIFT_effects_ext[,1], allStudiesCI); allStudiesCI_ext

    allStudiesCI_original_time_scale <- allStudiesCI_ext
    tmp1 <- grep("T0", colnames(allStudiesCI_original_time_scale)); tmp1
    tmp2 <- allStudiesCI_original_time_scale[, -c(1, tmp1)]; tmp2
    if (!(is.null(scaleTime))) scaleTime2 <- scaleTime else scaleTime2 <- 1
    if (is.null(ncol(tmp2))) tmp2b <- length(tmp2) else tmp2b <- ncol(tmp2)
    tmp3b <- matrix(round(as.numeric(tmp2) * scaleTime2, digits), ncol=tmp2b); tmp3b

    if (!(is.matrix(allStudiesCI_original_time_scale[, 1]))) {
      tmp3a <- matrix(allStudiesCI_original_time_scale[, 1], nrow=length(allStudiesCI_original_time_scale[, 1]))
    } else {
      tmp3a <- allStudiesCI_original_time_scale[, 1]
    }
    if (!(is.matrix(allStudiesCI_original_time_scale[, tmp1]))) {
      tmp3c <- matrix(allStudiesCI_original_time_scale[, tmp1], nrow=length(allStudiesCI_original_time_scale[, tmp1]))
    } else {
      tmp3c <- allStudiesCI_original_time_scale[, tmp1]
    }

    if (dim(tmp3c)[2] == 1) tmp3c <- c(tmp3c) # new 6.7.2022
    if (is.null(dim(tmp3c))) tmp3 <- cbind(tmp3a, tmp3b, t(tmp3c)) else tmp3 <- cbind(tmp3a, tmp3b,   tmp3c)
    colnames(tmp3) <- colnames(allStudiesCI_original_time_scale); tmp3
    allStudiesCI_original_time_scale <- tmp3

    # Label summary table
    rownames(allStudiesDRIFT_effects) <- paste0("Study No ", primaryStudies$studyNumbers)
    rownames(allStudiesDRIFT_effects_ext) <- paste0("Study No ", primaryStudies$studyNumbers)
    rownames(allStudiesDIFF_effects) <- paste0("Study No ", primaryStudies$studyNumbers)
    rownames(allStudiesDIFF_effects_ext) <- paste0("Study No ", primaryStudies$studyNumbers)

    if (!(is.null(scaleTime))) {
      #
      rownames(allStudiesDRIFT_effects_original_time_scale) <- paste0("Study No ", primaryStudies$studyNumbers)
      rownames(allStudiesDRIFT_effects_original_time_scale_ext) <- paste0("Study No ", primaryStudies$studyNumbers)
      #
      rownames(allStudiesDIFF_effects_original_time_scale) <- paste0("Study No ", primaryStudies$studyNumbers)
      rownames(allStudiesDIFF_effects_original_time_scale_ext) <- paste0("Study No ", primaryStudies$studyNumbers)
    }

    allStudiesDRIFT_effects_ext_dt <- allStudiesDRIFT_effects_ext
    tmp1 <- grep("toV", colnames(allStudiesDRIFT_effects_ext_dt))
    tmp2 <- (allStudiesDRIFT_effects_ext[, tmp1, drop=FALSE]); tmp2
    if (!(is.null(dim(tmp2)))) {
      for (l in 1:dim(tmp2)[1]) {
        tmp3 <- matrix(as.numeric(tmp2[l, ]), n.latent, n.latent, byrow=TRUE)
        if (is.null(scaleTime)) tmp4 <- c(t(OpenMx::expm(tmp3))) else tmp4 <- c(t(OpenMx::expm(tmp3 * scaleTime)))
        allStudiesDRIFT_effects_ext_dt[l, tmp1] <- round(tmp4, digits)
      }
    } else {
      tmp3 <- matrix(as.numeric(tmp2), n.latent, n.latent, byrow=TRUE); tmp3
      if (is.null(scaleTime)) tmp4 <- c(t(OpenMx::expm(tmp3))) else tmp4 <- c(t(OpenMx::expm(tmp3 * scaleTime)))
      allStudiesDRIFT_effects_ext_dt[tmp1]
      allStudiesDRIFT_effects_ext_dt[tmp1] <- round(tmp4, digits)
    }
    tmp1 <- grep("SE", colnames(allStudiesDRIFT_effects_ext_dt))
    allStudiesDRIFT_effects_ext_dt <- allStudiesDRIFT_effects_ext_dt[, -tmp1]
    if (!(is.null(dim(allStudiesDRIFT_effects_ext_dt)))) {
      colnames(allStudiesDRIFT_effects_ext_dt) <- paste0(colnames(allStudiesDRIFT_effects_ext_dt), " discrete time")
    } else {
      names(allStudiesDRIFT_effects_ext_dt)  <- paste0(names(allStudiesDRIFT_effects_ext_dt), " discrete time")
    }

    # check single study results
    if (checkSingleStudyResults == TRUE) {
      print(allStudiesDRIFT_effects_ext)
      print(allStudiesDRIFT_effects_ext_dt)
      print(allStudies_Minus2LogLikelihood)
      #

      cat(crayon::blue(" Press 'q' to quit and prevent saving or any other key to continue. Press ENTER afterwards."))
      char <- readline(" ")
      if (char == 'q' | char == 'Q') stop("Good luck for the next try!")
    }

    DRIFTCoeff <- matrix(unlist(model_Drift_Coef), n.studies, n.latent^2, byrow=TRUE); DRIFTCoeff
    DRIFTSE <- matrix(unlist(model_Drift_SE), n.studies, n.latent^2, byrow=TRUE); DRIFTSE
    DIFFCoeff <- matrix(unlist(model_Diffusion_Coef), n.studies, n.latent^2, byrow=TRUE); DIFFCoeff
    DIFFSE <- matrix(unlist(model_Diffusion_SE), n.studies, n.latent^2, byrow=TRUE); DIFFSE
    T0VARCoeff <- matrix(unlist(model_T0var_Coef), n.studies, n.latent^2, byrow=TRUE); T0VARCoeff
    T0VARSE <- matrix(unlist(model_T0var_SE), n.studies, n.latent^2, byrow=TRUE); T0VARSE

    if (n.studies < 2) {

      Msg <- "Only a single primary study was handed over to ctmaInitFit. No further (meta-) analyses can be conducted. \nI guess this stop is intended!"
      message(Msg)
    }
  } #end if (fit == TRUE)


  ##############################################################################################################
  if (fit == TRUE) {
    message1 <- c()



    if (!(all(is.na(unlist(primaryStudies$deltas))))) {
      maxDeltas <- max(unlist(primaryStudies$deltas), na.rm=TRUE)

      if (!(is.null(scaleTime))) maxDeltas <- maxDeltas * scaleTime
      largeDelta <- which(unlist(primaryStudies$deltas) >= maxDeltas); largeDelta

      if (!(is.null(scaleTime))) {
        tmp1 <- unlist(primaryStudies$deltas)
        tmp1 <- tmp1[!(is.na(tmp1))]
        largeDelta <- which( (tmp1 * scaleTime) >= maxDeltas)
      }
      tmp1 <- table(unlist(primaryStudies$deltas)[largeDelta]); tmp1
      if (!(is.null(scaleTime))) {
        table((tmp1 * scaleTime)[largeDelta])
      }
      tmp2 <- which(tmp1 == (max(tmp1))); tmp2
      suggestedScaleTime <- as.numeric(names(tmp1[tmp2])); suggestedScaleTime
      # maxDeltas hast already been time scaled above
      if (maxDeltas > 6) {
        tmp2 <- paste0("Maximum time interval was ", round(maxDeltas * scaleTime2, digits), "."); tmp2
        if (is.null(scaleTime)) scaleTime3 <- scaleTime2 else scaleTime3 <- scaleTime / scaleTime2
        tmp3 <- paste0("timeScale=1/", round(suggestedScaleTime * scaleTime3, digits)); tmp3 # CHD 19.11.2023
        if (suggestedScaleTime != scaleTime2) {
          tmp4 <- paste0("It is recommended to fit the model again using the arguments ", tmp3, ". "); tmp4
        } else {
          tmp4 <- paste0(""); tmp4
        }
        message1 <- paste(tmp2, tmp4, "If the model fit (-2ll) is better (lower), continue using, e.g.,", tmp3, "in all subsequent models.", collapse="\n"); message1
      }
    } else {
      maxDeltas <- NA
    }; maxDeltas

    if (is.null(scaleTime)) scaleTime2 <- 1 else scaleTime2 <- scaleTime

    if (!(is.null(scaleTime))) {
      model_Drift_Coef_original_time_scale <- lapply(model_Drift_Coef, function(x) x * scaleTime)
      model_Diffusion_Coef_original_time_scale <- lapply(model_Diffusion_Coef, function(x) x * scaleTime)
    } else {
      model_Drift_Coef_original_time_scale <- model_Drift_Coef
      model_Diffusion_Coef_original_time_scale <- model_Diffusion_Coef
    }

    if (randomIntercepts == TRUE) {
      randomIntercepts <- list(popsd=model_popsd,
                               popcov_mean=model_popcov_m,
                               popcov_sd=model_popcov_sd)
      # move T0Means of phantom latents to cints
      # needs still to be done
      skip <- 0
      if (skip == 1) {
        tmp1 <- grep("T0MEAN", rownames(invariantDrift_Coeff)); tmp1
        tmp1 <- tmp1[(n.latent+1):(n.latent^2)]; tmp1
        tmp1 <- invariantDrift_Coeff[tmp1, ]; tmp1
        tmp2 <- grep("CINT", rownames(invariantDrift_Coeff)); tmp2
        invariantDrift_Coeff[tmp2[1:n.latent],] <- tmp1
        # delete irrelevant rows
        toDelete <- c()
        for (i in (n.latent+1):(n.latent^2)) {
          toDelete <- c(toDelete, grep(i, rownames(invariantDrift_Coeff)))
        }
        invariantDrift_Coeff <- invariantDrift_Coeff[-toDelete, ]
        toDelete <- c()
        for (i in (n.latent+1):(n.latent^2)) {
          toDelete <- c(toDelete, grep(i, rownames(estimates_original_time_scale)))
        }
        estimates_original_time_scale <- estimates_original_time_scale[-toDelete, ]
      } # end skip
    } else {
      if ( (indVarying == 'CINT') | (indVarying == "MANIFEST") ){
        randomIntercepts <- list(note1="Covariances are time-scaled, correlations are unaffected by time scale.",
                                 note2="Undo time-scaling by multiplying the LR parts by (1/scaleTime)^2, and the LL and UR part by (1/scaleTime).",
                                 popsd=model_popsd,
                                 popcov_mean=model_popcov_m, model_popcov_sd=model_popcov_sd,
                                 model_popcov_T=model_popcov_T, model_popcov_025=model_popcov_025,
                                 model_popcov_50=model_popcov_50, model_popcov_975=model_popcov_975,
                                 popcor_mean=model_popcor_m, model_popcor_sd=model_popcor_sd,
                                 model_popcor_T=model_popcor_T, model_popcor_025=model_popcor_025,
                                 model_popcor_50=model_popcor_50, model_popcor_975=model_popcor_975)
      } else {
        randomIntercepts <- "Not estimated."
      }

    }

    results <- list(activeDirectory=activeDirectory,
                    backend=backend,
                    backendArgs=resolved_backend_args,
                    plot.type="drift", model.type=type,
                    coresToUse=coresToUse, n.studies=n.studies,
                    n.latent=n.latent,
                    n.manifest=n.manifest,
                    lambda=lambdaParams,
                    primaryStudyList=primaryStudies,
                    studyList=studyList, studyFitList=studyFit,
                    emprawList=empraw,
                    statisticsList=statisticsList,
                    ind.mod.List=empraw.ind.mod,
                    argumentList=list(ctmaInitFit=deparse(substitute(primaryStudies)),

                                      activeDirectory=activeDirectory,
                                      chains=chains,
                                      checkSingleStudyResults=checkSingleStudyResults,
                                      cint=cint,
                                      coresToUse=coresToUse,
                                      customPar=customPar,
                                      diff=diff,
                                      digits=digits,
                                      doPar=doPar,
                                      drift=drift,

                                      finishsamples=finishsamples,
                                      indVarying=indVarying,
                                      indVaryingT0=indVaryingT0,
                                      iter=iter,
                                      lambda=lambda,
                                      loadSingleStudyModelFit=loadSingleStudyModelFit,
                                      manifestMeans=manifestMeans,
                                      manifestVars=manifestVars,
                                      n.latent=n.latent,
                                      n.manifest=n.manifest,
                                      #nopriors=nopriors,
                                      optimize=optimize,
                                      primaryStudies=primaryStudies,
                                      priors=priors,
                                      sameInitialTimes=sameInitialTimes,
                                      saveRawData=saveRawData,
                                      saveSingleStudyModelFit=saveSingleStudyModelFit,
                                      scaleTI=scaleTI,
                                      scaleTime=scaleTime,
                                      silentOverwrite=silentOverwrite,
                                      T0means=T0means,
                                      T0var=T0var,
                                      useSV=useSV,
                                      verbose=verbose,
                                      randomIntercepts=randomInterceptsSettings,
                                      CoTiMAStanctArgs=CoTiMAStanctArgs,
                                      backend=backend,
                                      backendArgs=resolved_backend_args),
                    modelResults=list(DRIFT=model_Drift_Coef, DIFFUSION=model_Diffusion_Coef, T0VAR=model_T0var_Coef, CINT=model_Cint_Coef,
                                      DRIFToriginal_time_scale=model_Drift_Coef_original_time_scale,
                                      DIFFUSIONoriginal_time_scale=model_Diffusion_Coef_original_time_scale),
                    ctModel = currentModel,
                    #ProblemWithHessianEstimation = hessianWarning,
                    parameterNames=list(DRIFT=names(model_Drift_Coef[[1]]), DIFFUSION=names(model_Diffusion_Coef[[1]]), T0VAR=names(model_T0var_Coef[[1]])),
                    summary=(list(model="all drift free (het. model)",
                                  estimates=allStudiesDRIFT_effects_ext, #allStudiesDRIFT_effects_ext, = estimates that would be obtained without the scaleTime argument
                                  estimates_discrete=allStudiesDRIFT_effects_ext_dt,
                                  scaleTime=scaleTime2,
                                  # changed 17. Aug. 2022
                                  #randomIntercepts=list(popsd=model_popsd, popcov=model_popcov, popcor=model_popcor),
                                  randomIntercepts=randomIntercepts,
                                  confidenceIntervals=allStudiesCI,
                                  minus2ll= round(allStudies_Minus2LogLikelihood, digits),
                                  n.parameters = round(allStudies_estimatedParameters, digits),
                                  drift_estimates_original_time_scale =allStudiesDRIFT_effects_original_time_scale_ext,
                                  drift_CI_original_time_scale=allStudiesDriftCI_original_time_scale,
                                  diff_estimates_original_time_scale=allStudiesDIFF_effects_original_time_scale_ext,
                                  diff_CI_original_time_scale=allStudiesDiffCI_original_time_scale,
                                  estimationProblems=sapply(estProb, paste, sep="/n"),
                                  note=message1)))
    # excel workbook is added later
  } # end if (fit == TRUE)

  if (fit == FALSE) {
    #if (hessianWarning == FALSE) {
    results <- list(summary=c(
                      "No model was fitted, only data and code were generated. See $data, $ctModel, and $ctModelText.",
                      "To reconstruct the model from the generated R code, use eval(parse(text = your_result$ctModelText))."
                    ),
                    data = empraw,
                    ctModel = currentModel,
                    ctModelText = .ctmaModelText(currentModel),
                    backend = backend,
                    backendArgs = resolved_backend_args)
    #} else {
    #  results <- list(summary = list(message=c("No model was fitted, only data and code were generated. See $data & $ctModel section."),
    #                                 error=hessianWarning),
    #                  data = empraw,
    #                  ctModel = currentModel)
    #}
  }
  class(results) <- "CoTiMAFit"

  ### prepare Excel Workbook with several sheets ################################################################
  if (fit == TRUE) {
    {
      wb <- openxlsx::createWorkbook()
      sheet1 <- openxlsx::addWorksheet(wb, sheetName="model")
      sheet2 <- openxlsx::addWorksheet(wb, sheetName="modelResults")
      sheet3 <- openxlsx::addWorksheet(wb, sheetName="estimates")
      sheet4 <- openxlsx::addWorksheet(wb, sheetName="confidenceIntervals")
      sheet5 <- openxlsx::addWorksheet(wb, sheetName="randomIntercepts")
      sheet6 <- openxlsx::addWorksheet(wb, sheetName="stats")
      openxlsx::writeData(wb, sheet1, results$summary$model)

      ### modelResults
      # DRIFT
      startCol <- 2; startCol
      startRow <- 1; startRow
      openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, matrix(c(t(matrix(driftNames, nrow=n.latent))), nrow=1), colNames = FALSE)
      startCol <- 2; startCol
      startRow <- 2; startRow
      openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, colNames = FALSE,
                          matrix(unlist(results$modelResults$DRIFT),
                                 nrow=n.studies, ncol=n.latent^2, byrow=TRUE))
      startCol <- 1; startCol
      startRow <- 2; startRow
      openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, colNames = FALSE,
                          matrix(unlist(primaryStudies$studyNumbers), ncol=1))
      # DIFFUSION
      offset <-  n.studies + 1
      startCol <- 2; startCol
      startRow <- 2 + offset; startRow
      openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow,
                          matrix(names(results$modelResults$DIFFUSION[[1]]), nrow=1), colNames = FALSE)
      startCol <- 2; startCol
      startRow <- 2 + offset + 1# offset; startRow
      openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, colNames = FALSE,
                          matrix(unlist(results$modelResults$DIFFUSION),
                                 nrow=n.studies, byrow=TRUE))
      startCol <- 1; startCol
      startRow <- 2 + offset + 1; startRow
      openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, colNames = FALSE,
                          matrix(unlist(primaryStudies$studyNumbers), ncol=1))
      # T0Var
      offset <-  offset + n.studies + 2
      startCol <- 2; startCol
      startRow <- 2 + offset; startRow
      openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow,
                          matrix(names(results$modelResults$T0VAR[[1]]), nrow=1), colNames = FALSE)
      startCol <- 2; startCol
      startRow <- 2 + offset + 1# offset; startRow
      openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, colNames = FALSE,
                          matrix(unlist(results$modelResults$T0VAR),
                                 nrow=n.studies, byrow=TRUE))
      startCol <- 1; startCol
      startRow <- 2 + offset + 1; startRow
      openxlsx::writeData(wb, sheet2, startCol=startCol, startRow = startRow, colNames = FALSE,
                          matrix(unlist(primaryStudies$studyNumbers), ncol=1))
      ### estimates
      startCol <- 2; startCol
      startRow <- 1; startRow
      openxlsx::writeData(wb, sheet3, startCol=startCol, startRow = startRow,
                          t(colnames(results$summary$estimates)), colNames = FALSE)
      openxlsx::writeData(wb, sheet3, startCol=startCol, startRow = startRow + 1, results$summary$estimates, colNames = FALSE)
      ### confidence Intervals
      startCol <- 2; startCol
      startRow <- 1; startRow
      openxlsx::writeData(wb, sheet4, startCol=startCol, startRow = startRow,
                          t(colnames(results$summary$confidenceIntervals)), colNames = FALSE)
      openxlsx::writeData(wb, sheet4, startCol=startCol, startRow = startRow + 1, results$summary$confidenceIntervals, colNames = FALSE)
      ### random Intercepts
      startCol <- 2; startCol
      startRow <- 1; startRow
      openxlsx::writeData(wb, sheet5, startCol=startCol, startRow = startRow, "random intercept results are no longer reported in excel sheets", colNames = FALSE)
      ### stats
      startCol <- 2; startCol
      startRow <- 1; startRow
      tmp <- cbind("-2ll = ", results$summary$minus2ll, "Number of Parameters = ", results$summary$n.parameters)
      openxlsx::writeData(wb, sheet6, startCol=startCol, startRow = startRow, t(tmp), colNames = FALSE)
    }
    results$excelSheets <- wb
  }

  results <- .ctmaAttachStudyDesign(results, primaryStudies, emprawLong,
    currentModel, randomInterceptsSettings)
  if (isTRUE(fit)) {
    results <- .ctmaRefreshInitUncertainty(results)
    results <- .ctmaStoreOriginalTime(results, "initialization")
  }
  invisible(.ctmaStationaryResult(results, stationary, determinedDiffusions, stationaryVariance))

}   # end function definition
