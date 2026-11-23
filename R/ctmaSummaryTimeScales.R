# Summary-only conversion: never transform raw coordinates or stored fit draws.
.ctmaSummaryTimeFactor <- function(object) {
  stored <- if (is.list(object$summary)) object$summary else list()
  args <- if (is.list(object$argumentList)) object$argumentList else list()
  source <- if (is.list(object$source_object)) object$source_object else list()
  candidates <- list(args[["scaleTime", exact = TRUE]],
    stored[["scaledTime", exact = TRUE]], stored[["scaleTime", exact = TRUE]],
    object$CoTiMAStanctArgs[["scaleTime", exact = TRUE]],
    source$argumentList[["scaleTime", exact = TRUE]])
  candidates <- Filter(Negate(is.null), candidates)
  if (!length(candidates)) return(1)
  for (s in candidates) if (!is.numeric(s) || length(s) != 1L || !is.finite(s) || s <= 0)
    stop("The stored scaleTime must be one finite positive number.", call. = FALSE)
  if (any(vapply(candidates, function(s) !isTRUE(all.equal(as.numeric(s),
      as.numeric(candidates[[1L]]))), logical(1))))
    stop("Conflicting stored time scales; the original-time conversion is ambiguous.", call. = FALSE)
  as.numeric(candidates[[1L]])
}

.ctmaOriginalRateTable <- function(table, scale) {
  if (is.matrix(table) && is.numeric(table)) table <- as.data.frame(table)
  if (!is.data.frame(table) || !"matrix" %in% names(table)) return(NULL)
  # DIFFUSION alone may be a factor/coordinate matrix, not the covariance Q.
  keep <- tolower(as.character(table$matrix)) %in% c("drift", "diffusioncov")
  if (!any(keep)) return(NULL)
  out <- table[keep, , drop = FALSE]
  n <- tolower(names(out))
  location <- n %in% c("mean", "estimate", "est", "median", "mode", "sd", "se",
    "lower", "upper", "lwr", "upr", "ll", "ul", "ci.lower", "ci.upper",
    "mcse", "mcse_mean", "mcse_median") |
    grepl("^[0-9]+([.][0-9]+)?%?$", n)
  variance <- n %in% c("variance", "var")
  for (j in which(location | variance)) {
    if (!is.numeric(out[[j]])) stop("Rate-summary statistics must be numeric: ",
      names(out)[j], call. = FALSE)
    out[[j]] <- out[[j]] * if (variance[j]) scale^2 else scale
  }
  out
}

# Only scale columns that carry units. Diagnostic columns retain their values.
.ctmaScaleDriftSDTable <- function(table, scale) {
  if (!is.matrix(table) && !is.data.frame(table)) return(NULL)
  columns <- grepl("^(mean|sd|se|standard_error|median|mode|lower|upper|mcse|[0-9]+([.][0-9]+)?%?)$",
    colnames(table), ignore.case = TRUE)
  table[, columns] <- table[, columns, drop = FALSE] * scale
  table
}

.ctmaOriginalStudyDrift <- function(answer, native, scale) {
  if (!inherits(native, "ctJuliaFit")) return(NULL)
  model <- native$model
  if (!is.data.frame(model$pars)) model <- native$ctstanmodelbase
  pars <- model$pars
  if (!is.data.frame(pars) || !"indvarying_study" %in% names(pars)) return(NULL)
  labels <- unique(as.character(pars$param[
    pars$matrix %in% "DRIFT" & pars$indvarying_study %in% TRUE]))
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (!length(labels)) return(NULL)
  table <- answer[["popsd.study"]]
  if (is.null(table)) {
    # The unqualified table is safe only when DRIFT does not also vary by person.
    if (any(pars$indvarying %in% TRUE)) return(NULL)
    table <- answer[["popsd"]]
  }
  keep <- intersect(labels, rownames(table))
  if (!length(keep)) return(NULL)
  sd <- .ctmaScaleDriftSDTable(table[keep, , drop = FALSE], scale)
  out <- list(sd = sd)
  constrain <- get0(".ctBackendConstrained", envir = asNamespace("ctsem"),
    inherits = FALSE)
  if (is.function(constrain) && !is.null(native$sample)) {
    draws <- tryCatch({
      levels <- constrain(native)$randomeffectlevels
      study <- Filter(function(x) identical(x$level, "study"), levels)
      if (length(study) == 1L) study[[1L]]$popsd else NULL
    }, error = function(e) NULL)
    if (is.matrix(draws) && all(keep %in% colnames(draws)) && nrow(draws) > 1L) {
      variance <- t(vapply(keep, function(label) {
        values <- (draws[, label] * scale)^2
        if (!all(is.finite(values))) return(rep(NA_real_, 5L))
        c(mean(values), stats::sd(values),
          stats::quantile(values, c(.025, .5, .975), names = FALSE))
      }, numeric(5L)))
      colnames(variance) <- c("mean", "sd", "2.5%", "50%", "97.5%")
      out$variance <- variance
    }
  }
  if (is.null(out$variance)) out$note <-
    "Posterior variance summaries require saved study-level SD draws."
  out
}

.ctmaOriginalTime <- function(answer, object, scale, rates, interval) {
  if (scale == 1) return(NULL)
  fields <- c("drift_estimates_original_time_scale", "drift_CI_original_time_scale",
    "diff_estimates_original_time_scale", "diff_CI_original_time_scale",
    "mod_effects_original_time_scale", "clus_effects_original_time_scale",
    "WEC_estimates_original_time_scale", "opt.lag.orig.time")
  out <- list(scaleTime = scale, estimates = rates)
  for (field in intersect(fields, names(answer))) out[[field]] <- answer[[field]]
  if (!is.null(interval) && is.numeric(interval) && all(is.finite(interval))) out$timeinterval <- interval / scale
  re <- .ctmaOriginalStudyDrift(answer, object$studyFitList, scale)
  if (!is.null(re)) out$studyDrift <- re
  if (!any(vapply(out[-1L], Negate(is.null), logical(1)))) return(NULL)
  out
}

.ctmaStoreOriginalTime <- function(object, kind = "conventional", native_summary = NULL) {
  answer <- if (is.null(native_summary)) object$summary else native_summary
  converted <- .ctmaAddSummaryTimeScales(answer, object, kind)
  if (is.list(converted) && !is.null(converted$originalTime))
    object$summary$originalTime <- converted$originalTime
  object
}

.ctmaAddSummaryTimeScales <- function(answer, object, kind = "native", interval = NULL) {
  if (!is.list(answer) || is.data.frame(answer)) return(answer)
  scale <- .ctmaSummaryTimeFactor(object)
  old <- attr(answer, "ctma_time_reporting")
  if (!is.null(old)) {
    if (!isTRUE(all.equal(old$factor, scale))) stop(
      "This summary is already labelled with a different time scale.", call. = FALSE)
    return(answer)
  }
  before <- names(answer)
  fitted <- answer[["parmatrices", exact = TRUE]]
  if (is.null(fitted)) fitted <- answer[["estimates", exact = TRUE]]
  original <- answer[["estimates_original_time_scale", exact = TRUE]]
  if (is.null(original) && identical(kind, "initialization"))
    original <- answer[["drift_estimates_original_time_scale", exact = TRUE]]
  if (is.null(original)) original <- .ctmaOriginalRateTable(fitted, scale)
  answer["estimates_fitted_time_scale"] <- list(fitted)
  answer["estimates_original_time_scale"] <- list(original)
  info <- list(factor = scale, convention = "t_fitted = scaleTime * t_original",
    original_rate_factor = scale,
    original_table = if (identical(kind, "initialization"))
      "DRIFT; diffusion covariance remains in diff_estimates_original_time_scale" else
      "DRIFT and DIFFUSIONcov (covariance, not a factor matrix)",
    unchanged = "Other stored sections, raw parameters, likelihood and diagnostics are not rescaled.")
  if (!is.null(interval) && is.numeric(interval) && all(is.finite(interval))) {
    info$interval_fitted <- interval
    info$interval_original <- interval / scale
  }
  answer$time_scale_info <- info
  answer$originalTime <- .ctmaOriginalTime(answer, object, scale, original, interval)
  attr(answer, "ctma_time_reporting") <- list(factor = scale, kind = kind,
    added = setdiff(names(answer), before))
  class(answer) <- unique(c("summary.ctmaTimeScales", class(answer)))
  answer
}

#' @export
print.summary.ctmaTimeScales <- function(x, ...) {
  reporting <- attr(x, "ctma_time_reporting")
  s <- reporting$factor
  cat("Time scales: t_fitted = ", format(s, digits = 8), " * t_original.\n", sep = "")
  cat("CT rate estimates: fitted units; *_original_time_scale: original units.\n")
  view <- x
  hidden <- c(reporting$added, "originalTime",
    grep("_original_time_scale$|^opt\\.lag\\.orig\\.time$", names(view), value = TRUE))
  for (name in hidden) view[[name]] <- NULL
  class(view) <- setdiff(class(view), "summary.ctmaTimeScales")
  attr(view, "ctma_time_reporting") <- NULL
  if (reporting$kind != "native") {
    # Use explicit printed labels while preserving the original field aliases.
    names(view)[names(view) == "estimates"] <- "estimates_fitted_time_scale"
    names(view)[names(view) == "confidenceIntervals"] <- "confidenceIntervals_fitted_time_scale"
  }
  print(view, ...)
  if (!is.null(x$originalTime)) cat(
    "Original-time results are available in fit$summary$originalTime (new fits) ",
    "or in the summary() result under $originalTime.\n", sep = "")
  invisible(x)
}
