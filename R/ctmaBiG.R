#' ctmaBiG
#'
#' @description Analysis of publication bias and generalizability. The function takes a CoTiMA fit object (created with \code{\link{ctmaInit}})
#' and estimates common- and random-effects meta-analyses of individual drift coefficients, heterogeneity (Q, I square, H square, tau square),
#' PET-PEESE corrections, Egger's tests, and z-curve analysis yielding expected replication and detection rates (ERR, EDR).
#'
#'
#' @param ctmaInitFit fit object created with \code{\link{ctmaInit}} containing
#'   fitted primary-study models from the Stan or Julia backend. Julia fits
#'   must retain transformed ctsem draws.
#' @param activeDirectory the directory where to save results (if not specified, it is taken from ctmaInitFit)
#' @param PETPEESEalpha probability level (condition) below which to switch from PET to PEESE (cf. Stanley, 2017, p. 582, below Eq. 2; default p = .10)
#' @param digits rounding (default = 4)
#' @param finishsamples number of samples to draw (either from hessian based covariance or posterior distribution) for dt computations (default = 1000).
#' @param zcurve performs z-curve analysis. Could fail if too few studies (e.g. around 10) are supplied. default=FALSE
#' @param undoTimeScaling if TRUE, the original time scale is used (timeScale argument possibly used in \code{\link{ctmaInit}} is undone )
#' @param dt A scalar indicating a time interval across which discrete time effects should be estimated and then used for ctmaBiG.
#' @param excludeStudies vector with studies to be excluded from ctmaBiG.
#'
#' @importFrom stats var lm pnorm
#' @importFrom zcurve zcurve
#' @importFrom ctsem ctCollapse ctStanDiscretePars
#'
#' @export ctmaBiG
#'
#' @examples
#' \dontrun{
#' # perform analyses of publication bias and generalizability
#' CoTiMAInitFit_D_BO$activeDirectory <- "/Users/tmp/" # adapt!
#' CoTiMABiG_D_BO <- ctmaBiG(ctmaInitFit=CoTiMAInitFit_D_BO, zcurve=FALSE)
#' }
#'
#' @examples
#' # display results
#' \dontrun{
#' ctmaLoadExamples("CoTiMABiG_D_BO")
#' summary(CoTiMABiG_D_BO)
#' ctmaUnloadExamples("CoTiMABiG_D_BO")
#' }
#'
#' @examples
#' \dontrun{
#' # get funnel & forest plots
#' CoTiMABiG_D_BO$activeDirectory <- "/Users/tmp/" # adapt!
#' plot(CoTiMABiG_D_BO)
#' }
#' @return ctmaBiG returns a list containing some arguments supplied, the results of analyses of publication bias and generalizability,
#' model type, and the type of plot that could be performed with the returned object. The arguments in the returned object are activeDirectory,
#' and coresToUse. Further arguments, which are just copied from the init-fit object supplied, are, n.studies, n.latent, studyList,
#' statisticsList, modelResults (all parameter estimates and their standard error), and parameter names. All new results are returned
#' as the list element "summary", which is printed if the summary function is applied to the returned object. The summary list element
#' comprises a title (model='Analysis of Publication Bias & Generalizability') and "estimates", which is another list comprising
#' "Fixed Effects of Drift Coefficients", "Heterogeneity", "Random Effects of Drift Coefficients", "PET-PEESE corrections",
#' "Egger's tests" (constant of the WLS regression of drift coefficients on their standard errors (SE) with 1/SE^2 as weights),
#' "Egger's tests Alt. Version" (constant of the OLS regression of the standard normal deviates of the drift coefficients on their
#' precision), and "Z-Curve 2.0 Results". Plot type is plot.type=c("funnel", "forest") and model.type="BiG".
#'
#' @param ... Reserved. Unknown arguments are rejected with a migration hint where applicable.
ctmaBiG <- function(
    ctmaInitFit=NULL,
    activeDirectory=NULL,
    PETPEESEalpha=.10,

    digits=4,
    finishsamples=1000,
    zcurve=FALSE,
    undoTimeScaling=TRUE,
    dt=NULL,   # try BiG with dt effects. Specify the Time for which dT effects should analyzed. Always does it with original time scale
    excludeStudies=c(), ...)


{  # begin function definition (until end of file)
  .ctmaCheckRetiredArguments(match.call(expand.dots = FALSE), "ctmaBiG")
  .ctmaVersionNote()
  .ctma_timer <- .ctmaRegisterExecutionTimer("ctmaBiG")
  on.exit(.ctma_timer(), add = TRUE)

  ### check if fit object is specified
  if (is.null(ctmaInitFit)){

    ErrorMsg <- "\nA fitted CoTiMA (\"ctmaInitFit\") object has to be supplied to analyse something. \nGood luck for the next try!"
    stop(ErrorMsg)
  }

  if (ctmaInitFit$model.type == "mx") {

    Msg <- "I found old OpenMx-fitted Init fits. Will try ctmaBiGOMX! \n"
    message(Msg)
    result <- ctmaBiGOMX(ctmaInitFit=ctmaInitFit,
                         activeDirectory=activeDirectory,
                         PETPEESEalpha=PETPEESEalpha,

                         digits=digits)

    class(results) <- "CoTiMAFit"

    return(result)
    #stop("Please check results carefully")
  }

  if (length(ctmaInitFit$summary$estimationProblems) > 0L) {
    Msg <- "There were estimation problems while fitting with ctmaInit(), which probably prevent BiG analysis! \n"
    message(Msg)
    message(paste(ctmaInitFit$summary$estimationProblems, "\n"))
    Msg <- "Consider re-estimating or excluding studies from BiG analysis with the \"excludeStudies\" argument! \n"
    message(Msg)
  }

  #######################################################################################################################
  ############# Extracting Parameters from Fitted Primary Studies created with CoTiMAprep Function  #####################
  #######################################################################################################################

  start.time <- Sys.time(); start.time

  if (ctmaInitFit$model.type != "mx") {

    { # start extracting
      n.latent <- ctmaInitFit$n.latent; n.latent
      if (is.null(activeDirectory)) activeDirectory <- ctmaInitFit$activeDirectory; activeDirectory
      extracted_studies <- .ctmaBiGExtractStudies(
        ctmaInitFit, exclude_studies = excludeStudies)
      all_Coeff <- extracted_studies$estimates
      all_SE <- extracted_studies$standard_errors
      n.studies <- nrow(all_Coeff); n.studies
      names1 <- colnames(all_Coeff)[seq_len(n.latent^2)]
      excludeStudies <- extracted_studies$excluded_studies
      backend_counts <- table(
        extracted_studies$backends[extracted_studies$included_studies])
      message(
        "Extracted transformed study estimates from ",
        paste(paste(names(backend_counts), backend_counts, sep = "="),
          collapse = ", "),
        " backend fit(s)."
      )
      allSampleSizes <- ctmaInitFit$statisticsList$allSampleSizes; allSampleSizes
    } # end extracting
    #all_Coeff; all_SE


    # undo time scaling
    all_Coeff_timeScaled <- all_Coeff
    all_SE_timeScaled <- all_SE
    if (undoTimeScaling) {
      if(!(is.null(ctmaInitFit$summary$scaleTime))) {
        all_Coeff <- all_Coeff * ctmaInitFit$summary$scaleTime
        all_SE <- all_SE * ctmaInitFit$summary$scaleTime
      }
    }

    # CHD 24.2.2023
    .ctmaBiGCheckDriftSE(all_SE, n.latent, excludeStudies)

    # CHD 24.2.2023
    if (!(is.null(dt))) {
      #nsamples <- 1000
      drift_Coeff_dt <- matrix(NA, ncol=(1*(n.latent^2)), nrow=n.studies); drift_Coeff_dt
      drift_SE_dt <- matrix(NA, ncol=(1*(n.latent^2)), nrow=n.studies); drift_SE_dt
      #
      tmpTimeScale <- ctmaInitFit$summary$scaleTime; tmpTimeScale
      if (is.null(tmpTimeScale)) tmpTimeScale <- 1 # in case old fit files are used
      tmpTimeScale <- tmpTimeScale * dt; tmpTimeScale
      for (i in 1:n.studies) {
        #i <- 1
        #(!(i %in% excludeStudies))
        if (!(i %in% excludeStudies)) {
          tmpFit <- ctmaInitFit$studyFitList[[i]]
          tmpDrift_dt <- ctsem::ctStanDiscretePars(tmpFit,
                                                   subjects = "popmean",
                                                   times = tmpTimeScale,
                                                   nsamples = finishsamples,
                                                   plot=FALSE,
                                                   indices="ALL",
                                                   cores='maxneeded')
          #tmp <- cbind(tmpDrift_dt[ , 1, 1, , 1], tmpDrift_dt[ , 1, 1, , 2]) # columnwise
          tmp <- (tmpDrift_dt[,,,1:n.latent,1:n.latent])
          tmp <- cbind(tmp[,1,], tmp[,2,])
          dimnames(tmp)[[2]] <- c(matrix(names1, 2, 2, byrow=TRUE))
          drift_Coeff_dt[i, ] <- apply(tmp, 2, mean)
          drift_SE_dt[i, ] <- apply(tmp, 2, sd)
        }
      }
      colnames(drift_Coeff_dt) <- colnames(drift_SE_dt) <- names1
    }

    #######################################################################################################################
    ##################################### Analyses of Publication Bias ####################################################
    #######################################################################################################################


    print(paste0("#################################################################################"))
    print(paste0("################################# Egger's test ##################################"))
    print(paste0("#################################################################################"))

    DRIFTCoeff <- DIFFUSIONCoeff <- T0VARCoeff <- DRIFTSE <- DIFFUSIONSE <- T0VARSE <- matrix(NA, n.studies, n.latent^2)

    DRIFTCoeff <- all_Coeff[, grep("toV", colnames(all_Coeff))]; DRIFTCoeff
    DRIFTSE <- all_SE[, grep("toV", colnames(all_SE))]; DRIFTSE
    DIFFUSIONCoeff <- all_Coeff[, grep("diff", colnames(all_Coeff))]; DIFFUSIONCoeff
    DIFFUSIONSE <- all_SE[, grep("diff", colnames(all_SE))]; DIFFUSIONSE
    tmp <- grep("T0var", colnames(all_Coeff)); tmp
    if (length(tmp) == 0) tmp <- grep("T0cov", colnames(all_Coeff)); tmp
    T0VARCoeff <- all_Coeff[, tmp]; T0VARCoeff
    tmp <- grep("T0var", colnames(all_SE)); tmp
    if (length(tmp) == 0) tmp <- grep("T0cov", colnames(all_Coeff)); tmp
    T0VARSE <- all_SE[, tmp]; T0VARSE

    # just added in case time scaled funnel and forest plots should be done with ctmaPlot
    DRIFTCoeff_timeScaled <- all_Coeff_timeScaled[, grep("toV", colnames(all_Coeff))]; DRIFTCoeff_timeScaled
    DRIFTSE_timeScaled <- all_SE_timeScaled[, grep("toV", colnames(all_SE))]; DRIFTSE_timeScaled
    DIFFUSIONCoeff_timeScaled <- all_Coeff_timeScaled[, grep("diff", colnames(all_Coeff))]; DIFFUSIONCoeff_timeScaled
    DIFFUSIONSE_timeScaled <- all_SE_timeScaled[, grep("diff", colnames(all_SE))]; DIFFUSIONSE_timeScaled


    DRIFTCoeffSND <- DRIFTCoeff / DRIFTSE; DRIFTCoeffSND
    DRIFTPrecision <- c(rep(1, n.latent^2))/(DRIFTSE); DRIFTPrecision
    colnames(DRIFTPrecision) <- colnames(DRIFTCoeffSND); DRIFTPrecision
    if (!(is.null(dt))) {
      DRIFTPrecision_dt <- c(rep(1, n.latent^2))/(drift_SE_dt); DRIFTPrecision_dt
      colnames(DRIFTPrecision_dt) <- colnames(DRIFTPrecision_dt); DRIFTPrecision_dt
    }


    message1 <- "The pos. & sign. intercept indicates that SMALLER studies produced more positive (or less negative) effects"
    message2 <- "The neg. & sign. intercept indicates that LARGER studies produced more positive (or less negative) effects"
    tmp <- c()

    eggerDrift <- list()
    for (j in 1:(n.latent^2)) {
      #j <- 1
      tmp1 <- stats::lm(DRIFTCoeffSND[,j]~DRIFTPrecision[,j]) # This is identical to a weighted regression of drift on se ...
      tmp2 <- summary(tmp1)
      eggerDrift[[j]] <- list()
      eggerDrift[[j]]$message <- "No sign. evidence for publication bias."
      eggerDrift[[j]]$summary <- tmp2[[4]]
      #if (summary(eggerDrift[[j]])$coefficients[1,1] > 0 & summary(eggerDrift[[j]])$coefficients[1,4] < .05) {
      if (tmp2$coefficients[1,1] > 0 & tmp2$coefficients[1,4] < .05) {
        eggerDrift[[j]]$message <- message1
      }
      #if (summary(eggerDrift[[j]])$coefficients[1,1] < 0 & summary(eggerDrift[[j]])$coefficients[1,4] < .05) {
      if (tmp2$coefficients[1,1] < 0 & tmp2$coefficients[1,4] < .05) {
        eggerDrift[[j]]$message <- message2
      }
    }

    FREAResults <- list()
    FREAResults[[1]] <- "############# Eggers Test for DRIFT Parameter Estimates  ###############################"
    FREACounter <- 1
    for (j in 1:(n.latent^2)) {
      #j <- 1
      FREACounter <- FREACounter + 1
      FREAResults[[FREACounter]] <- paste0("-------------------------------- Eggers Test for ",
                                           colnames(DRIFTCoeff)[j], "--------------------------------")
      FREACounter <- FREACounter + 1
      FREAResults[[FREACounter]] <- eggerDrift[[j]]$message
      FREACounter <- FREACounter + 1
      #FREAResults[[FREACounter]] <- summary(eggerDrift[[j]])
      FREAResults[[FREACounter]] <- eggerDrift[[j]]$summary
    }


    ################################### Fixed & Random Effects Analyses ###################################################

    print(paste0("#################################################################################"))
    print(paste0("########## Fixed effect analysis of each drift coefficient separately  ##########"))
    print(paste0("#################################################################################"))

    meta_results <- .ctmaBiGMetaAnalysis(DRIFTCoeff, DRIFTSE)
    fixedEffectDriftResults <- meta_results$common
    heterogeneity <- meta_results$heterogeneity
    RandomEffectDriftResults <- meta_results$random
    FixedEffect_Drift <- fixedEffectDriftResults["FixedEffect_Drift", ]
    fixedEffectDriftMessage <- c()
    tau2DriftMessage <- c()
    if (any(meta_results$k < 3L)) {
      tau2DriftMessage <- paste0(
        "At least one parameter has fewer than three valid studies; ",
        "heterogeneity intervals or prediction intervals are unavailable."
      )
    }
    if (!(is.null(dt))) {
      meta_results_dt <- .ctmaBiGMetaAnalysis(drift_Coeff_dt, drift_SE_dt)
      fixedEffectDriftResults_dt <- meta_results_dt$common
      heterogeneity_dt <- meta_results_dt$heterogeneity
      RandomEffectDriftResults_dt <- meta_results_dt$random
      FixedEffect_Drift_dt <- fixedEffectDriftResults_dt["FixedEffect_Drift", ]
    }

    # RANDOM EFFECTS ANALYSIS ###############################################################################

    print(paste0("#################################################################################"))
    print(paste0("########## Random effect analysis of each drift coefficient separately  #########"))
    print(paste0("#################################################################################"))

    # Results were computed jointly above so all statistics use the same
    # parameter-specific set of valid studies.

    ### PET, PEESE & WLS approaches to correct for bias

    print(paste0("#################################################################################"))
    print(paste0("################ PET, PEESE & WLS approaches to correct for bias  ###############"))
    print(paste0("#################################################################################"))

    PETDrift_fit <- list()
    PEESEDrift_fit <- list()
    WLSDrift_fit <- list()
    PET_PEESEDrift_fit <- list()
    Egger2Drift_fit <- list()

    sampleSizes <- unlist(allSampleSizes); sampleSizes
    for (ii in 1:(n.latent^2)) {
      #ii <- 1
      driftCoeff <- DRIFTCoeff[ , ii]; driftCoeff
      driftSE <- DRIFTSE[, ii]; driftSE

      # PET; Egger's test = constant of the weigthed least squares regression of drift coefficients on their standard errors (SE) with 1/SE^2 as weights
      IV <- driftSE; IV
      DV <- driftCoeff; DV
      currentWeigths <- (1/(driftSE^2)); currentWeigths
      PETDrift_fit[[ii]] <- stats::lm(DV ~ IV, weights=currentWeigths); PETDrift_fit[[ii]]

      # Egger's Test (alternative but algebraically identical model)
      Egger2Drift_fit[[ii]] <- t(c(summary(PETDrift_fit[[ii]])$coefficients[2,1:4]))

      # PEESE
      IV <- driftSE^2; IV
      DV <- driftCoeff
      currentWeigths <- (1/(driftSE^2)); currentWeigths

      PEESEDrift_fit[[ii]] <- stats::lm(DV ~ IV, weights=currentWeigths); PEESEDrift_fit[[ii]]

      # PET-PEESE
      if ( (summary(PETDrift_fit[[ii]]))$coefficients[1,4] > PETPEESEalpha) {
        PET_PEESEDrift_fit[[ii]] <- PETDrift_fit[[ii]]
      }
      if ( (summary(PETDrift_fit[[ii]]))$coefficients[1,4] <= PETPEESEalpha) {
        PET_PEESEDrift_fit[[ii]] <- PEESEDrift_fit[[ii]]
      }

      # WLS
      DV <- driftCoeff/driftSE
      IV <- 1/driftSE
      WLSDrift_fit[[ii]] <- stats::lm(DV ~ IV + 0); WLSDrift_fit[[ii]]; FixedEffect_Drift[[ii]] # should be identical to FixedEffect_Drift
      #WLSDriftSE_fit <- summary(WLSDrift_fit[[ii]])$coefficients[2]; WLSDriftSE_fit; FixedEffect_DriftSE[[ii]] # should outperform FixedEffect_DriftSE
    }
    # same for dt
    if (!(is.null(dt))) {
      PETDrift_fit_dt <- list()
      PEESEDrift_fit_dt <- list()
      WLSDrift_fit_dt <- list()
      PET_PEESEDrift_fit_dt <- list()
      Egger2Drift_fit_dt <- list()

      #sampleSizes <- unlist(allSampleSizes); sampleSizes
      for (ii in 1:(n.latent^2)) {
        #ii <- 1
        driftCoeff <- drift_Coeff_dt[ , ii]; driftCoeff
        driftSE <- drift_SE_dt[, ii]; driftSE
        #driftCoeff <- driftCoeff[-excludeStudies]
        #driftSE <- driftSE[-excludeStudies]

        # PET; Egger's test = constant of the weigthed least squares regression of drift coefficients on their standard errors (SE) with 1/SE^2 as weights
        IV <- driftSE; IV
        DV <- driftCoeff; DV
        currentWeigths <- (1/(driftSE^2)); currentWeigths
        PETDrift_fit_dt[[ii]] <- stats::lm(DV ~ IV, weights=currentWeigths); PETDrift_fit_dt[[ii]]

        # Egger's Test (alternative but algebraically identical model)
        Egger2Drift_fit_dt[[ii]] <- t(c(summary(PETDrift_fit_dt[[ii]])$coefficients[2,1:4]))

        # PEESE
        IV <- driftSE^2; IV
        DV <- driftCoeff
        currentWeigths <- (1/(driftSE^2)); currentWeigths

        PEESEDrift_fit_dt[[ii]] <- stats::lm(DV ~ IV, weights=currentWeigths); PEESEDrift_fit_dt[[ii]]

        # PET-PEESE
        if ( (summary(PETDrift_fit_dt[[ii]]))$coefficients[1,4] > PETPEESEalpha) {
          PET_PEESEDrift_fit_dt[[ii]] <- PETDrift_fit_dt[[ii]]
        }
        if ( (summary(PETDrift_fit_dt[[ii]]))$coefficients[1,4] <= PETPEESEalpha) {
          PET_PEESEDrift_fit_dt[[ii]] <- PEESEDrift_fit_dt[[ii]]
        }

        # WLS
        DV <- driftCoeff/driftSE
        IV <- 1/driftSE
        WLSDrift_fit_dt[[ii]] <- stats::lm(DV ~ IV + 0); WLSDrift_fit_dt[[ii]]; FixedEffect_Drift_dt[[ii]] # should be identical to FixedEffect_Drift
        #WLSDriftSE_fit_dt <- summary(WLSDrift_fit_dt[[ii]])$coefficients[2]; WLSDriftSE_fit_dt; FixedEffect_DriftSE_dt[[ii]] # should outperform FixedEffect_DriftSE
      }
    }


    ############################################## zcurve Analysis ###################################################
    zFit <- zFit_dt <- "Set zcurve=TRUE for z-curve analysis results"
    if (zcurve == TRUE) {
      print(paste0("#################################################################################"))
      print(paste0("############################## Z-curve 2.0 analysis #############################"))
      print(paste0("#################################################################################"))
      print(paste0("                                                                                 "))
      print(paste0("Note: The fitting range is from |1.96| - |5|, as z-statistics > |5| have almost "))
      print(paste0("guaranteed > 99% power (the model treats them separately and includes them in   "))
      print(paste0("the EDR/ERR calculations, but does not fit the mixture of them, so too few      "))
      print(paste0("usable significant effects could perhaps be present in your data. In that case, "))
      print(paste0("you can be quite confident that the EDR and ERR will be close 1, but you have   "))
      print(paste0("to set the argument zcurve=FALSE!                                               "))

      zFit <- list()
      for (i in 1: dim(DRIFTCoeffSND)[2]) {
        tmp1 <- abs(DRIFTCoeffSND[, i]); tmp1
        tmp1
        zcurve::zcurve(z=tmp1)
        zFit[[i]] <- summary(zcurve::zcurve(z=tmp1))
      }
      names(zFit) <- paste0("Z-Curve 2.0 analysis of ", colnames(DRIFTCoeffSND)); zFit

      # same for dt
      if (!(is.null(dt))) {
        DRIFTCoeffSND_dt <- drift_Coeff_dt / drift_SE_dt; DRIFTCoeffSND_dt
        zFit_dt <- list()
        for (i in seq_len(ncol(DRIFTCoeffSND_dt))) {
          tmp1 <- abs(DRIFTCoeffSND_dt[, i]); tmp1
          zFit_dt[[i]] <- summary(zcurve::zcurve(z=tmp1))
        }
        names(zFit_dt) <- paste0(
          "Z-Curve 2.0 analysis of ", colnames(DRIFTCoeffSND_dt)); zFit_dt
      }

      ## format results
      #z.CurveResults <- matrix(NA, ncol=length(zFit), nrow=19)
      #for (h in 1:length(zFit)) {
      #  stopRow = 0
      #  for (j in 2:length(zFit[[i]])) {
      #    startRow  <- stopRow + 1; startRow
      #    stopRow <- startRow + length(unlist((zFit[[h]][j]))) - 1; stopRow
      #    unlist((zFit[[h]][j]))
      #    if (j == 2) z.CurveResults[startRow:stopRow, h] <- round(c((matrix(unlist((zFit[[h]][j])), ncol=2, byrow=TRUE))), digits)
      #    if (j != 2) z.CurveResults[startRow:stopRow, h] <- unlist((zFit[[h]][j]))
      #  }
      #}
      #rownamesPart <- c("ERR", "ERR lower CI", "ERR upper CI", "EDR", "EDR lower CI", "EDR upper CI" )
      #rownamesPart <- c(rownamesPart, names(unlist((zFit[[1]][3]))), names(unlist((zFit[[1]][4]))), names(unlist((zFit[[1]][5]))) )
      #rownames(z.CurveResults) <- rownamesPart
    }


    ############################################## Combine Results ###################################################
    {
      PET_Drift <- unlist(lapply(PETDrift_fit, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PET_Drift
      PET_SE <- c()
      for (k in 1:(n.latent^2)) PET_SE <- c(PET_SE, summary(PETDrift_fit[[k]])$coefficients[1,2])

      PEESE_Drift <-unlist(lapply(PEESEDrift_fit, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PEESE_Drift
      PEESE_SE <- c()
      for (k in 1:(n.latent^2)) PEESE_SE <- c(PEESE_SE, summary(PEESEDrift_fit[[k]])$coefficients[1,2])

      PET_PEESE_Drift <-unlist(lapply(PET_PEESEDrift_fit, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PET_PEESE_Drift
      PET_PEESE_SE <- c()
      for (k in 1:(n.latent^2)) PET_PEESE_SE <- c(PET_PEESE_SE, summary(PET_PEESEDrift_fit[[k]])$coefficients[1,2])

      WLS_Drift <- unlist(lapply(WLSDrift_fit, function(extract) extract$coefficients)); WLS_Drift
      WLS_SE <- c()
      for (k in 1:(n.latent^2)) WLS_SE <- c(WLS_SE, summary(WLSDrift_fit[[k]])$coefficients[1,2])

      Egger2Drift_results <- matrix(unlist(Egger2Drift_fit), ncol=n.latent^2, nrow=4); Egger2Drift_results

      PET_PEESE_DRIFTresults <- rbind(PET_Drift, PET_SE,
                                      PEESE_Drift, PEESE_SE,
                                      PET_PEESE_Drift, PET_PEESE_SE,
                                      WLS_Drift, WLS_SE,
                                      Egger2Drift_results)
      colnames(PET_PEESE_DRIFTresults) <- colnames(DRIFTCoeff)
      rownames(PET_PEESE_DRIFTresults) <- c(rownames(PET_PEESE_DRIFTresults)[1:8], "Egger's b0", "SE(b0)", "T", "p")
      #PET_PEESE_DRIFTresults

      ### some corrections for the output
      eggerTest <- PET_PEESE_DRIFTresults[9:12,]; eggerTest
      PET_PEESE_DRIFTresults <- PET_PEESE_DRIFTresults[1:8,]; PET_PEESE_DRIFTresults
    }
    # same for dt
    if (!(is.null(dt))) {
      PET_Drift_dt <- unlist(lapply(PETDrift_fit_dt, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PET_Drift_dt
      PET_SE_dt <- c()
      for (k in 1:(n.latent^2)) PET_SE_dt <- c(PET_SE_dt, summary(PETDrift_fit_dt[[k]])$coefficients[1,2])

      PEESE_Drift_dt <-unlist(lapply(PEESEDrift_fit_dt, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PEESE_Drift_dt
      PEESE_SE_dt <- c()
      for (k in 1:(n.latent^2)) PEESE_SE_dt <- c(PEESE_SE_dt, summary(PEESEDrift_fit_dt[[k]])$coefficients[1,2])

      PET_PEESE_Drift_dt <-unlist(lapply(PET_PEESEDrift_fit_dt, function(extract) extract$coefficients))[seq(1, 2*n.latent^2, 2)]; PET_PEESE_Drift_dt
      PET_PEESE_SE_dt <- c()
      for (k in 1:(n.latent^2)) PET_PEESE_SE_dt <- c(PET_PEESE_SE_dt, summary(PET_PEESEDrift_fit_dt[[k]])$coefficients[1,2])

      WLS_Drift_dt <- unlist(lapply(WLSDrift_fit_dt, function(extract) extract$coefficients)); WLS_Drift_dt
      WLS_SE_dt <- c()
      for (k in 1:(n.latent^2)) WLS_SE_dt <- c(WLS_SE_dt, summary(WLSDrift_fit_dt[[k]])$coefficients[1,2])

      Egger2Drift_results_dt <- matrix(unlist(Egger2Drift_fit_dt), ncol=n.latent^2, nrow=4); Egger2Drift_results_dt

      PET_PEESE_DRIFTresults_dt <- rbind(PET_Drift_dt, PET_SE_dt,
                                         PEESE_Drift_dt, PEESE_SE_dt,
                                         PET_PEESE_Drift_dt, PET_PEESE_SE_dt,
                                         WLS_Drift_dt, WLS_SE_dt,
                                         Egger2Drift_results_dt)
      colnames(PET_PEESE_DRIFTresults_dt) <- colnames(DRIFTCoeff)
      rownames(PET_PEESE_DRIFTresults_dt) <- c(rownames(PET_PEESE_DRIFTresults_dt)[1:8], "Egger's b0", "SE(b0)", "T", "p")
      ### some corrections for the output
      eggerTest_dt <- PET_PEESE_DRIFTresults_dt[9:12,]; eggerTest_dt
      PET_PEESE_DRIFTresults_dt <- PET_PEESE_DRIFTresults_dt[1:8,]; PET_PEESE_DRIFTresults_dt
    }

    # } # End Analysis of Publication Bias
    if (is.null(dt)) {
      modelResultsList <- list(DRIFT=DRIFTCoeff, DIFFUSION=DIFFUSIONCoeff, T0VAR=T0VARCoeff, CINT=NULL,
                               DRIFTSE=DRIFTSE, DIFFUSIONSE=DIFFUSIONSE, T0VARSE=T0VARSE,
                               DRIFT_timeScaled=DRIFTCoeff_timeScaled, DIFFUSION_timeScaled=DIFFUSIONCoeff_timeScaled,
                               DRIFTSE_timeScaled=DRIFTSE_timeScaled, DIFFUSIONSE_timeScaled=DIFFUSIONSE_timeScaled)
      summaryList <- list(model="Analysis of Publication Bias & Generalizability",
                          estimates=list("Fixed Effects of Drift Coefficients"=round(fixedEffectDriftResults, digits),
                                         "Heterogeneity"=round(heterogeneity, digits),
                                         "I2 message" = fixedEffectDriftMessage,
                                         "Tau2 message" = tau2DriftMessage,
                                         "Random Effects of Drift Coefficients"=round(RandomEffectDriftResults, digits),
                                         "PET-PEESE corrections"=round(PET_PEESE_DRIFTresults, digits),
                                         "Egger's tests"=round(eggerTest, digits),
                                         #"Egger's tests Alt. Version"= FREAResults,
                                         "Z-Curve 2.0 Results:"=zFit))
    }
    if (!(is.null(dt))) {
      modelResultsList <- list(DRIFT=DRIFTCoeff, DIFFUSION=DIFFUSIONCoeff, T0VAR=T0VARCoeff, CINT=NULL,
                               DRIFTSE=DRIFTSE, DIFFUSIONSE=DIFFUSIONSE, T0VARSE=T0VARSE,
                               DRIFT_timeScaled=DRIFTCoeff_timeScaled, DIFFUSION_timeScaled=DIFFUSIONCoeff_timeScaled,
                               DRIFTSE_timeScaled=DRIFTSE_timeScaled, DIFFUSIONSE_timeScaled=DIFFUSIONSE_timeScaled,
                               DRIFT_dt=drift_Coeff_dt,
                               DRIFT_dt_SE=drift_SE_dt)
      summaryList <- list(model="Analysis of Publication Bias & Generalizability",
                          estimates=list("Fixed Effects of Drift Coefficients"=round(fixedEffectDriftResults, digits),
                                         "Heterogeneity"=round(heterogeneity, digits),
                                         "I2 message" = fixedEffectDriftMessage,
                                         "Tau2 message" = tau2DriftMessage,
                                         "Random Effects of Drift Coefficients"=round(RandomEffectDriftResults, digits),
                                         "PET-PEESE corrections"=round(PET_PEESE_DRIFTresults, digits),
                                         "Egger's tests"=round(eggerTest, digits),
                                         #"Egger's tests Alt. Version"= FREAResults,
                                         "Z-Curve 2.0 Results:"=zFit,
                                         "discreteTimeTimeScale" = dt,
                                         "Fixed Effects of DISCRETE TIME Drift Coefficients"=round(fixedEffectDriftResults_dt, digits),
                                         "Heterogeneity in DISCRETE TIME "=round(heterogeneity_dt, digits),
                                         "Random Effects of DISCRETE TIME Drift Coefficients"=round(RandomEffectDriftResults_dt, digits),
                                         "PET-PEESE corrections IN DISCRETE TIME"=round(PET_PEESE_DRIFTresults_dt, digits),
                                         "Egger's tests"=round(eggerTest_dt, digits),
                                         "Z-Curve 2.0 Results in DISCRTE TIME:"=zFit_dt))
    }

    summaryList <- c(summaryList, "Supplied Studies"=n.studies)
    if (length(excludeStudies) > 0) summaryList <- c(summaryList, "Excluded Studies"=paste(excludeStudies, collapse=","), "Analysed Studies"=(n.studies-length(excludeStudies)))

    results <- list(activeDirectory=activeDirectory,
                    plot.type=c("funnel", "forest"), model.type="BiG",
                    backend=if (length(unique(stats::na.omit(
                      extracted_studies$backends))) == 1L) {
                      unique(stats::na.omit(extracted_studies$backends))
                    } else "mixed",
                    studyBackends=extracted_studies$backends,
                    coresToUse=NULL,
                    supplied.studies=n.studies, excludeStudies=excludeStudies, n.studies=(n.studies-length(excludeStudies)),
                    n.latent=n.latent,
                    studyList=ctmaInitFit$studyList, studyFitList=NULL, # , homDRIFTallFitCI),
                    emprawList=NULL,
                    statisticsList=ctmaInitFit$statisticsList,
                    modelResults=modelResultsList,
                    parameterNames=ctmaInitFit$parameterNames,
                    summary=summaryList)

    class(results) <- "CoTiMAFit"

    invisible(results)
  }
} ### END function definition
