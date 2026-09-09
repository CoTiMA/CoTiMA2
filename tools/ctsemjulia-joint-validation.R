#!/usr/bin/env Rscript

# Joint recovery test for backend = "ctsemjulia". Independent study effects
# are simulated for four parameter families, but ctsem estimates their full
# raw-scale correlation matrix.

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
.script_dir <- function() {
  x <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(x)) return(normalizePath("tools", mustWork = TRUE))
  dirname(normalizePath(sub("^--file=", "", x[[1L]]), mustWork = TRUE))
}
source(file.path(.script_dir(), "ctsemjulia-validation-suite.R"))

.joint_options <- function(args = commandArgs(trailingOnly = TRUE)) {
  x <- list(
    project = normalizePath(getwd(), mustWork = TRUE),
    output = path.expand("~/CoTiMA_Benchmark/ctsemjulia-joint"),
    development_library = "/private/tmp/ctsem-juliaFit-lib",
    julia_bindir = Sys.getenv("JULIA_BINDIR"), cores = 4L,
    mode = "quick", seed = 20260909L, study_effects = "random",
    retain_fit = FALSE, simulate_only = FALSE
  )
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) stop(
      "Arguments must use --name=value syntax.", call. = FALSE)
    z <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    name <- gsub("-", "_", z[[1L]], fixed = TRUE)
    if (!name %in% names(x)) stop("Unknown option: --", z[[1L]], call. = FALSE)
    x[[name]] <- paste(z[-1L], collapse = "=")
  }
  flag <- function(v) tolower(as.character(v)) %in% c("true", "t", "1", "yes")
  x$cores <- as.integer(x$cores); x$seed <- as.integer(x$seed)
  x$mode <- match.arg(tolower(x$mode), c("quick", "full"))
  x$study_effects <- match.arg(tolower(x$study_effects), c("random", "quantiles"))
  x$retain_fit <- flag(x$retain_fit); x$simulate_only <- flag(x$simulate_only)
  x
}

.joint_effects <- function(nstudy, nparameters, mode) {
  if (mode == "random") return(scale(
    matrix(stats::rnorm(nstudy * nparameters), nstudy, nparameters), TRUE, TRUE))
  if (nstudy <= nparameters) stop("Too few studies for orthogonal effects.")
  scale(stats::contr.helmert(nstudy)[, seq_len(nparameters), drop = FALSE],
    TRUE, TRUE)
}

.joint_simulate <- function(nstudy, subjects_per_study, nobs, seed,
    study_effects = c("random", "quantiles")) {
  set.seed(seed); study_effects <- match.arg(study_effects)
  spec <- data.frame(
    component = c("CINT", "DRIFT", "DIFFUSION", "T0VAR"),
    parameter = c("cint_eta1", "drift_eta1", "diff_eta1", "T0var_eta1"),
    raw_location = c(.01, .35, -1.70, -.80),
    raw_study_sd = c(.03, .30, .22, .25), stringsAsFactors = FALSE)
  transforms <- list(CINT = function(x) 10 * x, DRIFT = .ctma_drift_transform,
    DIFFUSION = .ctma_diffusion_transform, T0VAR = .ctma_t0var_transform)
  spec$true_reported_sd <- vapply(seq_len(nrow(spec)), function(i)
    .ctma_transformed_normal_sd(transforms[[spec$component[[i]]]],
      spec$raw_location[[i]], spec$raw_study_sd[[i]]), numeric(1L))
  z <- .joint_effects(nstudy, nrow(spec), study_effects)
  colnames(z) <- spec$component
  raw <- sweep(sweep(z, 2L, spec$raw_study_sd, `*`),
    2L, spec$raw_location, `+`)
  values <- cbind(CINT = 10 * raw[, "CINT"],
    DRIFT = .ctma_drift_transform(raw[, "DRIFT"]),
    DIFFUSION = .ctma_diffusion_transform(raw[, "DIFFUSION"]),
    T0VAR = .ctma_t0var_transform(raw[, "T0VAR"]))
  subject_raw_cint_sd <- .04; manifest_sd <- .10
  rows <- vector("list", nstudy * subjects_per_study); person <- 0L
  for (study in seq_len(nstudy)) {
    transition <- exp(values[study, "DRIFT"])
    q <- values[study, "DIFFUSION"]^2 / (-2 * values[study, "DRIFT"]) *
      (1 - transition^2)
    for (within in seq_len(subjects_per_study)) {
      person <- person + 1L
      cint <- values[study, "CINT"] + 10 * stats::rnorm(1L, sd = subject_raw_cint_sd)
      intercept <- cint / values[study, "DRIFT"] * (transition - 1)
      state <- stats::rnorm(1L, sd = values[study, "T0VAR"])
      y <- numeric(nobs)
      for (occasion in seq_len(nobs)) {
        if (occasion > 1L) state <- transition * state + intercept +
          stats::rnorm(1L, sd = sqrt(q))
        y[[occasion]] <- state + stats::rnorm(1L, sd = manifest_sd)
      }
      rows[[person]] <- data.frame(id = person, study = study,
        time = 0:(nobs - 1L), Y1 = y)
    }
  }
  list(data = do.call(rbind, rows), specification = spec,
    raw_study_values = raw, realised_raw_correlations = stats::cor(raw),
    study_effects = study_effects,
    true_subject_cint_sd = 10 * subject_raw_cint_sd)
}

.joint_run <- function(options) {
  dimensions <- if (options$mode == "quick")
    list(nstudy = 8L, subjects_per_study = 20L, nobs = 12L) else
    list(nstudy = 16L, subjects_per_study = 35L, nobs = 20L)
  design <- .joint_simulate(dimensions$nstudy, dimensions$subjects_per_study,
    dimensions$nobs, options$seed, options$study_effects)
  dir.create(options$output, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(design$specification,
    file.path(options$output, "joint-truth.csv"), row.names = FALSE)
  if (options$simulate_only) {
    saveRDS(design, file.path(options$output, "joint-simulation.rds"))
    message("Simulation written; fitting was skipped."); return(invisible(design))
  }
  if (dir.exists(options$development_library))
    .libPaths(c(options$development_library, .libPaths()))
  if (nzchar(options$julia_bindir)) Sys.setenv(JULIA_BINDIR = options$julia_bindir)
  Sys.setenv(JULIA_NUM_THREADS = as.character(options$cores))
  devtools::load_all(options$project, quiet = TRUE)
  status <- CoTiMA2:::.ctmaCtsemJuliaLaplaceStatus()
  if (!isTRUE(status$available)) stop(status$reason, call. = FALSE)
  ctsem::ctJuliaSetup(threads = options$cores, force = TRUE)
  model <- .ctma_validation_model(); model$pars$indvarying <- FALSE
  setup <- structure(list(data = design$data[, c("id", "time", "Y1")],
    ctModel = model, argumentList = list()), class = "CoTiMAFit")
  args <- CoTiMA2::ctmaBackendDefaults("ctsemjulia")
  args$study_id <- design$data$study
  args$study_parameters <- design$specification$parameter
  args$subject_parameters <- "cint_eta1"
  fit_time <- system.time(fit <- CoTiMA2:::.ctmaFitCtsemJuliaBackend(
    setup, args, options$cores, cotima_call = match.call(), verbose = 0))
  check_time <- system.time(check <- ctsem::ctLaplaceCheck(fit, nodes = 5L))
  s <- summary(fit)
  sd_table <- s[["popsd.study"]] %||% s[["popsd"]]
  cor_table <- s[["rawpopcorr.study"]] %||% s[["rawpopcorr"]]
  recovery <- design$specification
  recovery$estimated_reported_sd <- vapply(recovery$parameter, function(p)
    if (!is.null(sd_table) && p %in% rownames(sd_table))
      .ctma_scalar(sd_table[p, "mean"]) else NA_real_, numeric(1L))
  recovery$relative_error <- abs(recovery$estimated_reported_sd -
    recovery$true_reported_sd) / recovery$true_reported_sd
  interval <- t(vapply(recovery$parameter,
    function(p) .ctma_interval(sd_table, p), numeric(2L)))
  recovery$ci_lower <- interval[, "lower"]; recovery$ci_upper <- interval[, "upper"]
  recovery$covers_truth <- recovery$true_reported_sd >= recovery$ci_lower &
    recovery$true_reported_sd <= recovery$ci_upper
  correlations <- if (is.null(cor_table)) data.frame() else data.frame(
    pair = rownames(cor_table), estimated_raw_correlation = cor_table[, "mean"],
    row.names = NULL)
  subject_table <- s[["popsd.id"]]
  subject_estimate <- if (!is.null(subject_table) &&
    "cint_eta1" %in% rownames(subject_table))
      .ctma_scalar(subject_table["cint_eta1", "mean"]) else NA_real_
  result <- list(options = options, dimensions = dimensions,
    converged = isTRUE(fit$estimate$converged),
    gradient_norm = .ctma_scalar(fit$estimate$gradient_norm),
    loglik = .ctma_scalar(s$loglik), fit_seconds = unname(fit_time[["elapsed"]]),
    check_seconds = unname(check_time[["elapsed"]]),
    gap_per_subject = .ctma_scalar(check$gap_per_subject), recovery = recovery,
    correlations = correlations,
    realised_raw_correlations = design$realised_raw_correlations,
    true_subject_cint_sd = design$true_subject_cint_sd,
    estimated_subject_cint_sd = subject_estimate,
    fit = if (options$retain_fit) fit else NULL, check = check)
  utils::write.csv(recovery, file.path(options$output, "joint-recovery.csv"),
    row.names = FALSE)
  utils::write.csv(correlations,
    file.path(options$output, "joint-correlations.csv"), row.names = FALSE)
  saveRDS(result, file.path(options$output, "joint-validation-results.rds"))
  print(recovery, row.names = FALSE); print(correlations, row.names = FALSE)
  invisible(result)
}

if (sys.nframe() == 0L) .joint_run(.joint_options())
