#' Reduction of study-level random-effect heterogeneity
#'
#' @description
#' Compares an unconditional and a moderator-adjusted \code{ctmaReFit()} model
#' and estimates the proportional reduction in each study-level random-effect
#' variance. The conventional estimate is
#' \deqn{1 - \tau^2_{moderated} / \tau^2_{unconditional}.}
#' For Bayesian fits, uncertainty is described by independently combining
#' draws from the two fitted posterior distributions. This is a descriptive
#' pseudo-R-squared comparison, not a joint posterior or a causal variance
#' decomposition.
#'
#' @param unconditional A fitted \code{"ctmaStudyREFit"} without the focal
#'   study-level moderators.
#' @param moderated A fitted \code{"ctmaStudyREFit"} with the focal
#'   study-level moderators.
#' @param probs Lower, median, and upper probabilities used for the Monte
#'   Carlo interval. Defaults to \code{c(.025, .5, .975)}.
#' @param comparison_draws Number of independent cross-model draw combinations
#'   for Bayesian uncertainty summaries. The default is 10000.
#' @param seed Seed used for reproducible cross-model draw combinations.
#' @param practical_threshold Optional proportional reduction threshold. The
#'   returned probability describes how often the draw-wise reduction exceeds
#'   this value; the default is zero.
#'
#' @return An object of class \code{"ctmaReRedHet"} containing parameter-wise
#'   and overall raw-variance reductions, model metadata, moderator names, and
#'   the Monte Carlo comparison draws.
#' @export
#' @seealso \code{\link{ctmaRePrep}}, \code{\link{ctmaReFit}}
#'
#' @examples
#' \dontrun{
#' reduction <- ctmaReRedHet(
#'   unconditional = reFitWithoutModerators,
#'   moderated = reFitWithModerators
#' )
#' summary(reduction)
#' }
ctmaReRedHet <- function(
    unconditional,
    moderated,
    probs = c(.025, .5, .975),
    comparison_draws = 10000L,
    seed = 1234L,
    practical_threshold = 0) {
  .ctmaVersionNote()
  .ctma_timer <- .ctmaRegisterExecutionTimer("ctmaReRedHet")
  on.exit(.ctma_timer(), add = TRUE)

  comparison <- .ctmaReRedHetValidate(unconditional, moderated)
  probs <- as.numeric(probs)
  if (length(probs) != 3L || any(!is.finite(probs)) ||
      any(probs <= 0 | probs >= 1) || is.unsorted(probs, strictly = TRUE)) {
    stop("'probs' must contain three strictly increasing probabilities ",
      "between zero and one.", call. = FALSE)
  }
  comparison_draws <- as.integer(comparison_draws)[1L]
  if (is.na(comparison_draws) || comparison_draws < 1L) {
    stop("'comparison_draws' must be a positive integer.", call. = FALSE)
  }
  practical_threshold <- as.numeric(practical_threshold)[1L]
  if (!is.finite(practical_threshold)) {
    stop("'practical_threshold' must be finite.", call. = FALSE)
  }

  null_sd <- as.matrix(unconditional$draws$studyre_sd)[, comparison$null_order,
    drop = FALSE]
  moderated_sd <- as.matrix(moderated$draws$studyre_sd)[,
    comparison$moderated_order, drop = FALSE]
  if (any(!is.finite(null_sd)) || any(!is.finite(moderated_sd)) ||
      any(null_sd < 0) || any(moderated_sd < 0)) {
    stop("Both fits must contain finite non-negative studyre_sd draws.",
      call. = FALSE)
  }
  null_variance <- null_sd^2
  moderated_variance <- moderated_sd^2
  if (any(colMeans(null_variance) <= 0)) {
    stop("The unconditional fit estimates zero variance for at least one ",
      "random effect; proportional reduction is undefined.", call. = FALSE)
  }

  bayesian <- nrow(null_variance) > 1L && nrow(moderated_variance) > 1L
  if (bayesian) {
    if (!is.null(seed)) set.seed(as.integer(seed)[1L])
    null_index <- sample.int(nrow(null_variance), comparison_draws,
      replace = TRUE)
    moderated_index <- sample.int(nrow(moderated_variance), comparison_draws,
      replace = TRUE)
  } else {
    comparison_draws <- 1L
    null_index <- moderated_index <- 1L
  }
  reduction_draws <- 1 - moderated_variance[moderated_index, , drop = FALSE] /
    null_variance[null_index, , drop = FALSE]
  colnames(reduction_draws) <- comparison$keys

  summarise_reduction <- function(null_values, moderated_values,
      reduction_values) {
    interval <- if (length(reduction_values) > 1L) {
      stats::quantile(reduction_values, probs, na.rm = TRUE, names = FALSE)
    } else rep(reduction_values, 3L)
    c(
      unconditional_variance = mean(null_values),
      moderated_variance = mean(moderated_values),
      proportional_reduction = 1 - mean(moderated_values) /
        mean(null_values),
      lower = interval[1L],
      median = interval[2L],
      upper = interval[3L],
      probability_above_threshold = if (length(reduction_values) > 1L) {
        mean(reduction_values > practical_threshold, na.rm = TRUE)
      } else NA_real_
    )
  }

  spec <- comparison$specification
  parameter_results <- do.call(rbind, lapply(seq_len(nrow(spec)), function(i) {
    statistics <- summarise_reduction(
      null_variance[, i], moderated_variance[, i], reduction_draws[, i])
    data.frame(
      re = spec$re[i], matrix = spec$matrix[i], row = spec$row[i],
      col = spec$col[i], parameter = spec$parameter[i], block = spec$block[i],
      scale = "raw study-effect variance", t(statistics),
      check.names = FALSE, row.names = NULL
    )
  }))

  null_total <- rowSums(null_variance)
  moderated_total <- rowSums(moderated_variance)
  total_reduction_draws <- 1 -
    moderated_total[moderated_index] / null_total[null_index]
  overall <- as.data.frame(as.list(summarise_reduction(
    null_total, moderated_total, total_reduction_draws)),
    check.names = FALSE)
  overall$measure <- "sum of marginal raw study-effect variances"
  overall <- overall[, c("measure", setdiff(names(overall), "measure"))]

  structure(list(
    parameter_reduction = parameter_results,
    overall_reduction = overall,
    practical_threshold = practical_threshold,
    probs = probs,
    comparison_draws = if (bayesian) reduction_draws else NULL,
    overall_comparison_draws = if (bayesian) total_reduction_draws else NULL,
    unconditional_moderators = .ctmaReRedHetModeratorNames(unconditional),
    moderated_moderators = .ctmaReRedHetModeratorNames(moderated),
    backends = c(
      unconditional = if (is.null(unconditional$backend)) "stan" else
        unconditional$backend,
      moderated = if (is.null(moderated$backend)) "stan" else
        moderated$backend
    ),
    estimations = c(unconditional = unconditional$estimation,
      moderated = moderated$estimation),
    study_labels = unconditional$prepared$study_labels,
    specification = spec,
    call = match.call()
  ), class = "ctmaReRedHet")
}

.ctmaReRedHetKey <- function(specification) {
  do.call(paste, c(specification[c(
    "matrix", "row", "col", "parameter", "type", "transform", "block"
  )], sep = "\r"))
}

.ctmaReRedHetValidate <- function(unconditional, moderated) {
  if (!inherits(unconditional, "ctmaStudyREFit") ||
      !inherits(moderated, "ctmaStudyREFit")) {
    stop("Both inputs must be fitted objects returned by ctmaReFit().",
      call. = FALSE)
  }
  if (is.null(unconditional$draws$studyre_sd) ||
      is.null(moderated$draws$studyre_sd)) {
    stop("Both fits must contain studyre_sd estimates.", call. = FALSE)
  }
  null_draws <- nrow(as.matrix(unconditional$draws$studyre_sd))
  moderated_draws <- nrow(as.matrix(moderated$draws$studyre_sd))
  if (xor(null_draws > 1L, moderated_draws > 1L)) {
    stop("Compare either two posterior fits or two point-estimate fits; one ",
      "input contains posterior draws and the other does not.", call. = FALSE)
  }
  if (!identical(as.character(unconditional$prepared$study_labels),
      as.character(moderated$prepared$study_labels))) {
    stop("The two fits must contain the same studies in the same order.",
      call. = FALSE)
  }
  null_spec <- unconditional$prepared$specification
  moderated_spec <- moderated$prepared$specification
  null_keys <- .ctmaReRedHetKey(null_spec)
  moderated_keys <- .ctmaReRedHetKey(moderated_spec)
  if (anyDuplicated(null_keys) || anyDuplicated(moderated_keys) ||
      !setequal(null_keys, moderated_keys)) {
    stop("The two fits must use the same random-effect specification, ",
      "transformations, and covariance blocks.", call. = FALSE)
  }
  moderated_order <- match(null_keys, moderated_keys)
  if (!identical(
      unconditional$prepared$covariance_blocks$structure,
      moderated$prepared$covariance_blocks$structure)) {
    stop("The two fits must use the same random-effect covariance structure.",
      call. = FALSE)
  }
  if (!isTRUE(all.equal(unconditional$prepared$priors,
      moderated$prepared$priors))) {
    warning("The two fits use different random-effect priors; interpret the ",
      "variance-reduction comparison cautiously.", call. = FALSE)
  }
  list(
    specification = null_spec,
    keys = null_keys,
    null_order = seq_along(null_keys),
    moderated_order = moderated_order
  )
}

.ctmaReRedHetModeratorNames <- function(fit) {
  prepared <- fit$prepared
  setup <- prepared$standata$TIPREDEFFECTsetup
  ti_names <- prepared$original$ctstanmodelbase$TIpredNames
  spec <- prepared$specification
  rows <- unique(spec$index[spec$type == 1L])
  if (is.null(setup) || !length(rows) || !length(ti_names)) {
    return(character())
  }
  active_columns <- which(colSums(setup[rows, , drop = FALSE] > 0L) > 0L)
  setdiff(as.character(ti_names[active_columns]), prepared$study_dummy_names)
}

#' @export
print.ctmaReRedHet <- function(x, digits = 4L, ...) {
  cat("ctmaReRedHet: reduction in study-level random-effect variance\n\n")
  print(x$parameter_reduction, digits = digits, row.names = FALSE)
  cat("\nOverall\n")
  print(x$overall_reduction, digits = digits, row.names = FALSE)
  invisible(x)
}

#' @export
summary.ctmaReRedHet <- function(object, ...) object
