#' ctmaOptimizeInit
#'
#' @description Initial fitting (i.e., applies \code{\link{ctmaInit}}) to a primary study reFit times to capitalize on chance for obtaining
#' a hard-to-find optimal fit. This could be very helpful if a primary yields out-of-range estimates, which could happen if the fitting
#' algorithm unfortunately used random start values that resulted in a locally but not globally optimal fit. Essentially, using
#' ctmaOptimizeInit is like gambling, hoping that at leas one set of starting values (the number is tries is specified in the reFits argument)
#' eneables finding the global optimal fit. On unix-like machines (e.g. MacOS), this could be done in parallel mode if coresToUse > 1.
#'
#' @param primaryStudies list of primary study information created with \code{\link{ctmaPrep}} or \code{\link{ctmaFitToPrep}}
#' @param activeDirectory activeDirectory
#' @param problemStudy number (position in list) where the problem study in primaryStudies is found
#' @param reFits how many reFits should be done
#' @param n.latent number of latent variables of the model (hast to be specified)!
#' @param coresToUse if neg., the value is subtracted from available cores, else value = cores to use
#' @param finishsamples number of samples to draw (either from hessian based covariance or posterior distribution) for final results computation (default = 1000).
#' @param indVarying control for unobserved heterogeneity by having randomly (inter-individually) varying manifest means
#' @param randomScaleTime lower and upper limit of uniform distribution from which timeScale argument for ctmaInit is uniformly shuffled (integer)
#' @param customPar logical. If set TRUE (default) leverages the first pass using priors and ensure that the drift diagonal cannot easily go too negative (helps since ctsem > 3.4)
#' @param T0means Default 0 (assuming standardized variables). Can be assigned labels to estimate them freely.
#' @param manifestMeans Default 0 (assuming standardized variables). Can be assigned labels to estimate them freely.
#' @param CoTiMAStanctArgs parameters that can be set to improve model fitting of the \code{\link[ctsem]{ctStanFit}} Function
#' @param scaleTime scale time (interval) - sometimes desirable to improve fitting
#' @param manifestVars define the error variances of the manifests with a single time point using R-type lower triangular matrix with nrow=n.manifest & ncol=n.manifest.
#' @param checkSingleStudyResults displays estimates from single study 'ctsem' models and waits for user input to continue.
#' Useful to check estimates before they are saved.


#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar%
#' @importFrom stats runif
#'
#' @note Unix parallel retries use owned fork workers without registering a
#'   global foreach backend; Windows uses sequential retries. Input errors and
#'   retry failures leave the caller's foreach registration unchanged.
#' @note During fitting, not output is generated. Be patient.
#'
#' @examples
#' \dontrun{
#' optimFit313 <- ctmaOptimizeInit(primaryStudies=CoTiMAstudyList_3,
#'                                 activeDirectory="/Users/tmp/",  # adapt!
#'                                 problemStudy=which(CoTiMAstudyList_3$studyNumbers == 313),
#'                                 reFits=10,
#'                                 n.latent=2)
#' summary(optimFit313)
#' }
#'
#' @export ctmaOptimizeInit
#'
#' @return returns a list with bestFit (= the best fit achieved), all_minus2ll (= all -2ll values for all fitted models), and summary, which
#' is printed if the summary function is applied to the returned object, and which shows the summary information of the ctsem model with the
#' best fit.
#'
#' @param ... Reserved. Unknown arguments are rejected with a migration hint where applicable.
ctmaOptimizeInit <- function(primaryStudies=NULL,
                             activeDirectory=NULL,
                             problemStudy=NULL,
                             reFits=NULL,
                             finishsamples=NULL,
                             n.latent=NULL,
                             coresToUse=c(1),
                             indVarying=FALSE,
                             randomScaleTime=c(1,1),

                             checkSingleStudyResults=FALSE,
                             customPar=FALSE,
                             T0means=0,
                             manifestMeans=0,
                             manifestVars=NULL,
                             CoTiMAStanctArgs=NULL,
                             scaleTime=NULL, ...)
{
  .ctmaCheckRetiredArguments(match.call(expand.dots = FALSE), "ctmaOptimizeInit")

  # Reject missing/invalid inputs before touching resources or caller state.
  if (is.null(primaryStudies) || is.null(problemStudy) || is.null(reFits) ||
      is.null(activeDirectory) || is.null(n.latent)) stop("arguments are missing")
  if (!is.numeric(reFits) || length(reFits) != 1L || !is.finite(reFits) ||
      reFits < 1 || reFits != floor(reFits)) stop("reFits must be one positive integer.")
  if (!is.numeric(randomScaleTime) || length(randomScaleTime) != 2L ||
      any(!is.finite(randomScaleTime)) || randomScaleTime[2] < randomScaleTime[1])
    stop("randomScaleTime must contain two finite ordered values.")

  #######################################################################################################################
  ################################################# Check Cores To Use ##################################################
  #######################################################################################################################

  {
    coresToUse <- .ctmaResolveCores(coresToUse, backend = "stan",
      workflow = "ctmaOptimizeInit")$resolved
    .ctma_restore_threads <- .ctmaLimitThreadEnvironment()
    on.exit(.ctma_restore_threads(), add = TRUE)

  }


  if (!is.null(CoTiMAStanctArgs)) {
    .Deprecated("backendArgs", package = "CoTiMA2",
      old = "CoTiMAStanctArgs")
  }
  CoTiMAStanctArgs <- .ctmaMergeBackendArgs(
    .ctmaStanBackendDefaults(),
    if (is.null(CoTiMAStanctArgs)) list() else CoTiMAStanctArgs,
    path = "CoTiMAStanctArgs"
  )
  if (!(is.null(finishsamples))) CoTiMAStanctArgs$optimcontrol$finishsamples <- finishsamples
  if (!(is.null(scaleTime))) CoTiMAStanctArgs$scaleTime <- scaleTime

  ########################################################################################################################

  ErrorMsg <- "arguments are missing"
  if (is.null(primaryStudies) | is.null(problemStudy) | is.null(reFits) | is.null(activeDirectory) | is.null(n.latent)  ) stop(ErrorMsg)

  if (randomScaleTime[2] < randomScaleTime[1]) {

    ErrorMsg <- "\nrandomScaleTime[1] has to be <= randomScaleTime[2]! \nGood luck for the next try!"
    stop(ErrorMsg)
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

  #
  # adaptations for dealing with raw data
  # taken out 9. Aug. 2022
  #if (!(is.null(newStudyList$rawData[[1]]$studyNumbers))) {
  #  newStudyList$rawData[[1]]$studyNumbers <- 1
  #  newStudyList$studyNumbers <- 1 # otherwise raw data will not be loaded
  #  newStudyList$deltas <- unlist(newStudyList$deltas)
  #}
  #newStudyList

  # parallel re-fitting of problem study
  allfits <- .ctmaOptimizeInitApply(seq_len(reFits), function(i) {
    scaleTime <- round(stats::runif(1, min=randomScaleTime[1], max=randomScaleTime[2]), 2)
    fits <- ctmaInit(primaryStudies=newStudyList,
                     coresToUse = 1, n.latent=n.latent,
                     indVarying = indVarying,
                     scaleTime = scaleTime,
                     activeDirectory = activeDirectory,
                     checkSingleStudyResults=checkSingleStudyResults,
                     customPar=customPar,
                     T0means=T0means,
                     manifestMeans=manifestMeans,
                     manifestVars=manifestVars)
    fits
  }, coresToUse)
  all_minus2ll <- lapply(allfits, function(x) x$summary$minus2ll)
  bestFit <- which(unlist(all_minus2ll) == min(unlist(all_minus2ll)))[1]; bestFit
  bestFit <- allfits[[bestFit]]

  results <- list(bestFit=bestFit, all_minus2ll=all_minus2ll, summary=bestFit$summary,
                  resultsSummary=bestFit$studyFitList[[1]]$resultsSummary
  )
  class(results) <- "CoTiMAFit"

  invisible(results)
}


# Use owned workers rather than replacing the user's global foreach backend.
.ctmaOptimizeInitApply <- function(X, FUN, cores) {
  worker <- function(i) tryCatch(list(ok = TRUE, value = FUN(i)),
    error = function(e) list(ok = FALSE, error = conditionMessage(e)))
  results <- if (.Platform$OS.type == "unix" && cores > 1L)
    parallel::mclapply(X, worker, mc.cores = cores, mc.cleanup = TRUE,
      mc.allow.recursive = FALSE) else lapply(X, worker)
  failed <- which(vapply(results, function(x) !is.list(x) || !isTRUE(x$ok), logical(1)))
  if (length(failed)) {
    i <- failed[1L]; result <- results[[i]]
    detail <- if (is.list(result) && is.character(result$error)) result$error else
      "worker did not return a result"
    stop("Initialization retry ", X[i], " failed: ", detail, call. = FALSE)
  }
  lapply(results, `[[`, "value")
}
