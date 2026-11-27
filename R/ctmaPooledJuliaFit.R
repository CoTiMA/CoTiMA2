# Keep the engine fit unchanged inside the established pooled CoTiMA container.
.ctmaWrapPooledJuliaFit <- function(fit, setup) {
  if (!is.list(fit) || !is.list(setup)) stop("A native fit and matching setup are required.", call. = FALSE)
  out <- fit
  out$studyFitList <- fit
  for (name in c("data", "ctModel", "ctModelText", "studyList", "emprawList",
      "statisticsList", "moderatorCoding", "n.latent", "activeDirectory", "model.type", "plot.type")) {
    value <- setup[[name, exact = TRUE]]
    if (!is.null(value)) out[[name]] <- value
  }
  if (is.null(out$n.latent)) out$n.latent <- setup$ctModel$n.latent
  if (is.null(out$plot.type)) out$plot.type <- "drift"
  if (is.null(out$activeDirectory)) out$activeDirectory <- setup$argumentList$activeDirectory
  out$n.studies <- 1L
  if (!is.list(out$summary)) {
    ll <- fit$estimate$loglik
    if (is.null(ll)) ll <- fit$optimization$loglik
    if (is.null(ll)) ll <- fit$likelihood$julia
    if (is.null(ll)) ll <- NA_real_
    raw <- fit$estimate$raw
    if (is.null(raw)) raw <- fit$coefficients$raw
    out$summary <- list(loglik = as.numeric(ll)[1L],
      minus2ll = -2 * as.numeric(ll)[1L],
      n.parameters = if (is.null(raw)) NA_integer_ else length(raw))
  }
  # The wrapper must not inherit an engine class: dispatch and extraction must
  # distinguish the pooled container from its native studyFitList object.
  class(out) <- c("CoTiMAPooledJuliaFit", "CoTiMAFit")
  if (.ctmaSummaryTimeFactor(out) != 1) {
    native_summary <- tryCatch(summary(fit, digits = 8), error = function(e) {
      warning("Could not store original-time results: ", conditionMessage(e),
        call. = FALSE)
      NULL
    })
    if (!is.null(native_summary)) out <- .ctmaStoreOriginalTime(out, "native",
      native_summary)
  }
  out
}

.ctmaPooledJuliaNative <- function(object) {
  native <- object$studyFitList
  if (!inherits(native, c("ctJuliaFit", "ctmaJuliaFit"))) stop(
    "This pooled object does not contain a supported native Julia fit.", call. = FALSE)
  native
}

#' Methods for pooled Julia CoTiMA fits
#'
#' Pooled Julia results returned by ctmaFit() retain the full native fit in
#' studyFitList. Summaries retain the native fields and additionally report
#' transformed DRIFT and DIFFUSIONcov estimates on the original time scale.
#' Raw parameters and diagnostics are not rescaled. Other methods delegate
#' to the unchanged engine object.
#' @param pointEstimate NULL inherits the fit's policy; auto, mean or mode select the point estimator.
#' @param x,object A pooled Julia CoTiMA fit.
#' @param ... Arguments forwarded to the corresponding native-fit method.
#' @return The native summary plus \code{estimates_fitted_time_scale},
#' \code{estimates_original_time_scale} and \code{time_scale_info}; or the native
#' coefficients or covariance matrix. Native DT intervals retain their fitted
#' time units; \code{time_scale_info} records the corresponding original intervals. Printing
#'   returns the pooled object invisibly.
#' @name pooledJuliaFit
#' @seealso \code{\link{ctmaFit}}, \code{\link{ctmaFitDiagnostics}}
NULL

#' @rdname pooledJuliaFit
#' @export
print.CoTiMAPooledJuliaFit <- function(x, ...) {
  cat("CoTiMA pooled fit; native fit stored in $studyFitList\n")
  print(.ctmaPooledJuliaNative(x), ...)
  invisible(x)
}

#' @rdname pooledJuliaFit
#' @export
summary.CoTiMAPooledJuliaFit <- function(object, ..., pointEstimate = NULL) {
  .ctmaSummaryTimeFactor(object) # Validate metadata before delegating.
  native <- .ctmaPooledJuliaNative(object)
  if (is.null(pointEstimate)) pointEstimate <- object$argumentList$pointEstimate
  if (is.null(pointEstimate)) pointEstimate <- "auto"
  answer <- summary(native, ...)
  reporting <- native
  reporting$argumentList <- object$argumentList
  answer <- .ctmaPointSummary(answer,reporting,pointEstimate)
  args <- list(...)
  interval <- if (inherits(native, "ctJuliaFit")) {
    if (is.null(args$timeinterval)) 1 else args$timeinterval
  } else NULL
  .ctmaAddSummaryTimeScales(answer, object, "native", interval)
}

#' @rdname pooledJuliaFit
#' @export
coef.CoTiMAPooledJuliaFit <- function(object, ...) {
  stats::coef(.ctmaPooledJuliaNative(object), ...)
}

#' @rdname pooledJuliaFit
#' @export
vcov.CoTiMAPooledJuliaFit <- function(object, ...) {
  stats::vcov(.ctmaPooledJuliaNative(object), ...)
}
