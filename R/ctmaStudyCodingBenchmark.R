#' Controlled Stan comparison of study coding schemes
#'
#' @description Runs matched short Bayesian Stan pilots with dummy and weighted
#' study coding. Both fits use the same model arguments, iteration controls,
#' and random seed. The returned diagnostics compare runtime and NUTS geometry;
#' short pilots are not intended for substantive posterior inference. When the
#' input is simulated with identical true effects in every study, both coding
#' schemes should describe the same population model; this is especially useful
#' for isolating computational effects of the parameterization. Such a scenario
#' does not evaluate performance under genuine between-study heterogeneity.
#'
#' @param ctmaInitFit A fitted \code{ctmaInit} object supplied to
#'   \code{\link{ctmaFit}}.
#' @param studyVarying Parameter families allowed to vary between studies.
#' @param fitArgs Named list of additional model arguments passed unchanged to
#'   both calls to \code{ctmaFit}, for example \code{binaries},
#'   \code{randomIntercepts}, and \code{scaleTime}.
#' @param iter Total Stan iterations per chain.
#' @param warmup Stan warm-up iterations per chain. Must be smaller than
#'   \code{iter}.
#' @param chains Number of matched Stan chains.
#' @param coresToUse Number of parallel chains.
#' @param seed Integer seed reset before each coding condition.
#' @param adapt_delta Stan target acceptance probability.
#' @param max_treedepth Maximum NUTS tree depth.
#' @param backendArgs Additional Stan backend controls. The function supplies
#'   \code{warmup}, \code{adapt_delta}, and \code{max_treedepth}; do not repeat
#'   those entries here.
#' @param retain_fits If \code{TRUE}, retain both potentially large fitted
#'   objects in the result. The default keeps diagnostics and summaries only.
#' @param verbose Verbosity passed to \code{ctmaFit}.
#'
#' @return An object of class \code{ctmaStudyCodingBenchmark} containing the
#'   matched controls, a comparison table, parameter diagnostics, timings, and
#'   optionally the fitted objects.
#' @export ctmaStudyCodingBenchmark
#'
#' @examples
#' \dontrun{
#' coding_benchmark <- ctmaStudyCodingBenchmark(
#'   ctmaInitFit = initFit,
#'   studyVarying = c("DIFFUSION", "T0VAR"),
#'   fitArgs = list(
#'     binaries = c(1, 0),
#'     randomIntercepts = "CINT",
#'     scaleTime = 1 / 12
#'   )
#' )
#' coding_benchmark
#' }
ctmaStudyCodingBenchmark <- function(
    ctmaInitFit,
    studyVarying = c("DIFFUSION", "T0VAR"),
    fitArgs = list(),
    iter = 200L,
    warmup = 150L,
    chains = 4L,
    coresToUse = min(chains, 4L),
    seed = 20260907L,
    adapt_delta = 0.80,
    max_treedepth = 10L,
    backendArgs = list(),
    retain_fits = FALSE,
    verbose = 0L
) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("ctmaStudyCodingBenchmark requires the suggested package 'rstan'.",
         call. = FALSE)
  }
  if (missing(ctmaInitFit) || is.null(ctmaInitFit)) {
    stop("'ctmaInitFit' must be supplied.", call. = FALSE)
  }
  if (!is.list(fitArgs) || (length(fitArgs) &&
      (is.null(names(fitArgs)) || any(!nzchar(names(fitArgs)))))) {
    stop("'fitArgs' must be a fully named list.", call. = FALSE)
  }
  if (!is.list(backendArgs)) {
    stop("'backendArgs' must be a named list.", call. = FALSE)
  }

  integer_control <- function(x, name, minimum = 1L) {
    value <- suppressWarnings(as.integer(x)[1L])
    if (is.na(value) || value < minimum) {
      stop("'", name, "' must be an integer >= ", minimum, ".",
           call. = FALSE)
    }
    value
  }
  iter <- integer_control(iter, "iter")
  warmup <- integer_control(warmup, "warmup", minimum = 0L)
  chains <- integer_control(chains, "chains")
  coresToUse <- integer_control(coresToUse, "coresToUse")
  seed <- integer_control(seed, "seed", minimum = 0L)
  max_treedepth <- integer_control(max_treedepth, "max_treedepth")
  if (warmup >= iter) {
    stop("'warmup' must be smaller than 'iter'.", call. = FALSE)
  }
  adapt_delta <- suppressWarnings(as.numeric(adapt_delta)[1L])
  if (!is.finite(adapt_delta) || adapt_delta <= 0 || adapt_delta >= 1) {
    stop("'adapt_delta' must lie strictly between 0 and 1.", call. = FALSE)
  }

  protected_fit <- c(
    "ctmaInitFit", "backend", "optimize", "priors", "studyCoding",
    "studyVarying", "chains", "iter", "coresToUse", "verbose",
    "backendArgs", "CoTiMAStanctArgs"
  )
  overlap <- intersect(names(fitArgs), protected_fit)
  if (length(overlap)) {
    stop(
      "Controlled arguments must not be supplied in 'fitArgs': ",
      paste(overlap, collapse = ", "), ".",
      call. = FALSE
    )
  }
  protected_backend <- c("warmup")
  backend_overlap <- intersect(names(backendArgs), protected_backend)
  if (length(backend_overlap)) {
    stop("Do not repeat controlled backend arguments: ",
         paste(backend_overlap, collapse = ", "), ".", call. = FALSE)
  }
  if (!is.null(backendArgs$control) && !is.list(backendArgs$control)) {
    stop("'backendArgs$control' must be a list.", call. = FALSE)
  }
  controlled_names <- c("adapt_delta", "max_treedepth")
  repeated_control <- intersect(names(backendArgs$control), controlled_names)
  if (length(repeated_control)) {
    stop("Do not repeat controlled Stan controls: ",
         paste(repeated_control, collapse = ", "), ".", call. = FALSE)
  }
  stan_args <- backendArgs
  stan_args$warmup <- warmup
  stan_args$control <- utils::modifyList(
    if (is.null(backendArgs$control)) list() else backendArgs$control,
    list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  )

  run_one <- function(coding) {
    set.seed(seed)
    call <- c(
      list(
        ctmaInitFit = ctmaInitFit,
        backend = "stan",
        optimize = FALSE,
        priors = TRUE,
        studyCoding = coding,
        studyVarying = studyVarying,
        chains = chains,
        iter = iter,
        coresToUse = coresToUse,
        verbose = verbose,
        backendArgs = stan_args
      ),
      fitArgs
    )
    elapsed <- system.time(fit <- do.call(ctmaFit, call))
    stanfit <- .ctmaStudyCodingStanfit(fit)
    summary_matrix <- rstan::summary(stanfit)$summary
    sampler_params <- rstan::get_sampler_params(stanfit, inc_warmup = FALSE)
    diagnostics <- .ctmaSummarizeStanPilot(
      sampler_params = sampler_params,
      summary_matrix = summary_matrix,
      max_treedepth = max_treedepth,
      elapsed_seconds = unname(elapsed[["elapsed"]])
    )
    list(
      coding = coding,
      diagnostics = diagnostics,
      parameter_diagnostics = .ctmaSubstantiveStanSummary(summary_matrix),
      elapsed = elapsed,
      fit = if (isTRUE(retain_fits)) fit else NULL
    )
  }

  results <- lapply(c("dummy", "weighted"), run_one)
  names(results) <- c("dummy", "weighted")
  comparison <- do.call(rbind, lapply(results, function(x) {
    data.frame(coding = x$coding, x$diagnostics, check.names = FALSE)
  }))
  rownames(comparison) <- NULL
  weighted_speedup <- comparison$seconds_per_postwarmup_iteration[
    comparison$coding == "dummy"
  ] / comparison$seconds_per_postwarmup_iteration[
    comparison$coding == "weighted"
  ]

  structure(
    list(
      call = match.call(),
      controls = list(
        studyVarying = studyVarying, fitArgs = fitArgs, iter = iter,
        warmup = warmup, postwarmup = iter - warmup, chains = chains,
        coresToUse = coresToUse, seed = seed, adapt_delta = adapt_delta,
        max_treedepth = max_treedepth, backendArgs = stan_args
      ),
      interpretation = paste(
        "This is a computational comparison, not a substantive posterior",
        "analysis. If simulated studies share identical true effects, coding",
        "differences should not create systematic differences in the fitted",
        "population effects; genuine heterogeneity requires a separate test."
      ),
      comparison = comparison,
      weighted_speedup_vs_dummy = unname(weighted_speedup),
      results = results
    ),
    class = "ctmaStudyCodingBenchmark"
  )
}

.ctmaStudyCodingStanfit <- function(fit) {
  candidate <- tryCatch(fit$studyFitList$stanfit$stanfit,
                        error = function(e) NULL)
  if (is.null(candidate) || !methods::is(candidate, "stanfit")) {
    stop("The fitted CoTiMA object does not contain the expected Stan fit.",
         call. = FALSE)
  }
  candidate
}

.ctmaSubstantiveStanSummary <- function(summary_matrix) {
  keep <- !grepl("^(llrow\\[|lp__$)", rownames(summary_matrix))
  summary_matrix[keep, , drop = FALSE]
}

.ctmaSummarizeStanPilot <- function(
    sampler_params,
    summary_matrix,
    max_treedepth,
    elapsed_seconds
) {
  combined <- do.call(rbind, sampler_params)
  substantive <- .ctmaSubstantiveStanSummary(summary_matrix)
  rhat <- substantive[, "Rhat"]
  neff <- substantive[, "n_eff"]
  n_draws <- nrow(combined)
  postwarmup_per_chain <- if (length(sampler_params)) {
    nrow(sampler_params[[1L]])
  } else {
    0L
  }
  list(
    elapsed_seconds = elapsed_seconds,
    chains = length(sampler_params),
    postwarmup_per_chain = postwarmup_per_chain,
    seconds_per_postwarmup_iteration =
      elapsed_seconds / max(1L, postwarmup_per_chain),
    divergences = sum(combined[, "divergent__"]),
    treedepth_saturated = sum(combined[, "treedepth__"] >= max_treedepth),
    treedepth_saturation_fraction =
      mean(combined[, "treedepth__"] >= max_treedepth),
    mean_treedepth = mean(combined[, "treedepth__"]),
    mean_leapfrog = mean(combined[, "n_leapfrog__"]),
    median_leapfrog = stats::median(combined[, "n_leapfrog__"]),
    mean_acceptance = mean(combined[, "accept_stat__"]),
    largest_substantive_rhat = max(rhat, na.rm = TRUE),
    median_substantive_rhat = stats::median(rhat, na.rm = TRUE),
    smallest_substantive_neff = min(neff, na.rm = TRUE),
    median_substantive_neff = stats::median(neff, na.rm = TRUE),
    substantive_parameters = nrow(substantive),
    total_postwarmup_draws = n_draws
  )
}

#' @export
print.ctmaStudyCodingBenchmark <- function(x, ...) {
  cat("Controlled Stan study-coding benchmark\n")
  cat("======================================\n")
  print(x$comparison, row.names = FALSE)
  cat("\nWeighted speedup versus dummy coding: ",
      format(x$weighted_speedup_vs_dummy, digits = 4L), "x\n", sep = "")
  cat(paste0(
    "Interpret short pilots as computational diagnostics only. ",
    "High R-hat, low effective sample size, or frequent maximum-tree-depth ",
    "events require a longer or better-parameterized fit.\n"
  ))
  cat(x$interpretation, "\n")
  invisible(x)
}
