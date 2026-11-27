#' Summary method for CoTiMA fit objects
#'
#' @description
#' Provides the standard summary for conventional CoTiMA fit objects and
#' redirects moment-based and hierarchical study-level random-effects fits to
#' their corresponding summary methods.
#'
#' @param pointEstimate NULL inherits the fit's policy. auto selects optimized mode for ML/MAP where available; mean selects earlier draw-mean reporting; mode requires optimized values.
#' @param object A fitted object inheriting from class \code{"CoTiMAFit"}.
#' @param probs Numeric vector containing the lower, median, and upper
#'   probabilities used for random-effects posterior summaries.
#' @param ... Further arguments passed to the applicable summary method.
#'
#' @return
#' For a conventional CoTiMA fit, the stored \code{object$summary} is printed
#' and returned, with fitted- and original-time aliases and explicit labels.
#' The original-scale fields retain their established meanings; drift and
#' diffusion covariance rates, SEs and intervals use the recorded time factor.
#' For a hierarchical study-level random-effects fit, a list
#' containing population, study-specific, raw-scale and transformed-scale
#' heterogeneity, correlation, and Time-0 parameter summaries is returned.
#' The hierarchical summary distinguishes the transformed raw-scale population
#' centre from the equal-study average calculated after transforming each
#' study-specific parameter within every posterior draw. An interpretation
#' note is printed and retained in \code{estimand_interpretation}.
#' Julia primary-study fits additionally print separate stored fit checks;
#' use \code{ctmaFitDiagnostics(object)} to retrieve these diagnostics.
#' Identification notes are displayed as a compact study table and a single
#' explanation per flag. Full messages remain available with
#' \code{as.character(object$summary$identification_notes)} or
#' \code{print(object$summary$identification_notes, details = TRUE)}.
#'
#' @method summary CoTiMAFit
#' @export
#'
summary.CoTiMAFit <- function(
    object,
    ...,
    probs = c(.025, .5, .975), pointEstimate = NULL) {

  if (!inherits(object, "CoTiMAFit")) {
    stop("Not a CoTiMAFit object!", call. = FALSE)
  }

  # Post-hoc empirical moment estimates returned by ctmaReMoments().
  if (inherits(object, "CoTiMAReMoments")) {
    return(.ctmaReMomentsSummary(
      object,
      probs = object$probs,
      print_note = TRUE
    ))
  }

  # Hierarchical study-level random-effects fits returned by ctmaReFit()
  # or ctmaReMLFit().
  if (inherits(object, "CoTiMAStudyREFit")) {

    if (!exists(".ctmaReSummary", mode = "function")) {
      stop(
        "The internal CoTiMA random-effects summary function ",
        "'.ctmaReSummary()' is not available.",
        call. = FALSE
      )
    }

    return(
      .ctmaReSummary(
        object,
        probs = probs
      )
    )
  }

  if(is.null(pointEstimate)) pointEstimate <- object$argumentList$pointEstimate
  if(is.null(pointEstimate)) pointEstimate <- "auto"
  pointEstimate <- .ctmaPointChoice(pointEstimate)
  object <- .ctmaPointInitSummary(object,pointEstimate)
  object$summary <- .ctmaPointSummary(object$summary,object,pointEstimate)
  # Refresh reporting from stored native tables only: no engine or refit.
  object <- .ctmaRefreshInitUncertainty(object)
  # Keep established fields; add consistently labelled time-scale aliases.
  kind <- if (is.list(object$studyFitList) && length(object$n.studies) == 1L &&
    length(object$studyFitList) == object$n.studies &&
    !.ctmaIsNativeFit(object$studyFitList)) "initialization" else "conventional"
  object$summary <- .ctmaAddSummaryTimeScales(object$summary, object, kind)
  print(object$summary)
  if (.ctmaHasJuliaFit(object)) {
    ctmaFitDiagnostics(object)
  }
  if (!is.null(object$summary$opt.lag.orig.time)) {
    message(
      "Note: The reported optimal time lag (OTL) uses the stored point-estimate ",
      "drift matrix. For draw-based median OTLs and ",
      "uncertainty intervals, use ctmaOTL(fit, CI = TRUE)."
    )
  }
  invisible(object$summary)
}
