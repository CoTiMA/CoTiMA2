#!/usr/bin/env Rscript

# Resumable Monte Carlo controller for ctsemjulia-validation-suite.R.
# Each replication runs in an independent R/Julia process. This avoids
# forking a live Julia runtime and works on both macOS and Linux.

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

.mc_options <- function(args = commandArgs(trailingOnly = TRUE)) {
  defaults <- list(
    project = normalizePath(getwd(), mustWork = TRUE),
    output = path.expand("~/CoTiMA_Benchmark/ctsemjulia-monte-carlo"),
    development_library = "",
    julia_bindir = Sys.getenv("JULIA_BINDIR"),
    mode = "full", scenarios = "", seed = 20260909L,
    replications = 10L, replication_start = 1L,
    workers = 1L, cores_per_fit = 4L, retain_fits = FALSE,
    aggregate_only = FALSE
  )
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) {
      stop("Arguments must use --name=value syntax.", call. = FALSE)
    }
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    name <- gsub("-", "_", parts[[1L]], fixed = TRUE)
    if (!name %in% names(defaults)) stop("Unknown option: --", parts[[1L]], call. = FALSE)
    defaults[[name]] <- paste(parts[-1L], collapse = "=")
  }
  integers <- c("seed", "replications", "replication_start", "workers", "cores_per_fit")
  defaults[integers] <- lapply(defaults[integers], as.integer)
  logicals <- c("retain_fits", "aggregate_only")
  defaults[logicals] <- lapply(defaults[logicals], function(x) {
    tolower(as.character(x)) %in% c("true", "t", "1", "yes")
  })
  defaults$mode <- match.arg(tolower(defaults$mode), c("quick", "full"))
  if (any(vapply(defaults[integers[-1L]], function(x) is.na(x) || x < 1L, logical(1L)))) {
    stop("Replication and worker counts must be positive integers.", call. = FALSE)
  }
  defaults
}

.mc_quote_arg <- function(name, value) shQuote(paste0("--", name, "=", value))

.mc_run_replication <- function(replication, options, runner, rscript) {
  replication_dir <- file.path(options$output, "replications", sprintf("rep-%04d", replication))
  log_file <- file.path(options$output, "logs", sprintf("rep-%04d.log", replication))
  dir.create(replication_dir, recursive = TRUE, showWarnings = FALSE)
  seed <- as.integer(options$seed + replication * 1000L)
  args <- c(
    shQuote(runner),
    .mc_quote_arg("project", options$project),
    .mc_quote_arg("output", replication_dir),
    .mc_quote_arg("cores", options$cores_per_fit),
    .mc_quote_arg("mode", options$mode),
    .mc_quote_arg("seed", seed),
    .mc_quote_arg("retain-fits", tolower(options$retain_fits)),
    .mc_quote_arg("resume", "true")
  )
  if (nzchar(options$scenarios)) args <- c(args, .mc_quote_arg("scenarios", options$scenarios))
  if (nzchar(options$development_library)) {
    args <- c(args, .mc_quote_arg("development-library", options$development_library))
  }
  if (nzchar(options$julia_bindir)) {
    args <- c(args, .mc_quote_arg("julia-bindir", options$julia_bindir))
  }
  started <- Sys.time()
  status <- tryCatch(system2(rscript, args, stdout = log_file, stderr = log_file),
    error = function(e) structure(1L, error = conditionMessage(e)))
  data.frame(
    replication = replication, seed = seed, status = as.integer(status),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    log = log_file, stringsAsFactors = FALSE
  )
}

.mc_read_results <- function(options) {
  files <- list.files(file.path(options$output, "replications"), pattern = "^[a-z0-9]+[.]rds$",
    recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("validation-results[.]rds$", files)]
  rows <- lapply(files, function(file) {
    x <- tryCatch(readRDS(file), error = function(e) NULL)
    if (is.null(x) || is.null(x$scenario)) return(NULL)
    rep_name <- basename(dirname(file))
    replication <- suppressWarnings(as.integer(sub("^rep-", "", rep_name)))
    fields <- c("scenario", "status", "seed", "converged", "fit_seconds", "check_seconds",
      "gap_per_subject", "true_raw_study_sd", "true_study_sd",
      "estimated_study_sd", "study_sd_relative_error",
      "study_ci_covers_truth", "true_subject_cint_sd", "estimated_subject_cint_sd",
      "subject_cint_ci_covers_truth")
    out <- setNames(rep(list(NA), length(fields)), fields)
    for (name in intersect(fields, names(x))) out[[name]] <- x[[name]]
    data.frame(replication = replication, out, stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

.mc_mean <- function(x) if (any(is.finite(x))) mean(x, na.rm = TRUE) else NA_real_
.mc_median <- function(x) if (any(is.finite(x))) median(x, na.rm = TRUE) else NA_real_

.mc_summarize <- function(raw, requested_replications) {
  if (!nrow(raw)) return(data.frame())
  split_rows <- split(raw, raw$scenario)
  do.call(rbind, lapply(split_rows, function(x) {
    ok <- x$status == "completed" & is.finite(x$estimated_study_sd)
    error <- x$estimated_study_sd - x$true_study_sd
    truth <- unique(x$true_study_sd[is.finite(x$true_study_sd)])
    subject_error <- x$estimated_subject_cint_sd - x$true_subject_cint_sd
    data.frame(
      scenario = x$scenario[[1L]], requested = requested_replications,
      observed = nrow(x), completed = sum(x$status == "completed"),
      converged = sum(x$converged %in% TRUE),
      convergence_rate = sum(x$converged %in% TRUE) / requested_replications,
      true_raw_study_sd = {
        value <- unique(x$true_raw_study_sd[is.finite(x$true_raw_study_sd)])
        if (length(value)) value[[1L]] else NA_real_
      },
      true_study_sd = if (length(truth)) truth[[1L]] else NA_real_,
      mean_estimated_study_sd = .mc_mean(x$estimated_study_sd),
      bias = .mc_mean(error), relative_bias = .mc_mean(error / x$true_study_sd),
      rmse = sqrt(.mc_mean(error^2)), mae = .mc_mean(abs(error)),
      median_relative_error = .mc_median(x$study_sd_relative_error),
      coverage_95 = if (any(!is.na(x$study_ci_covers_truth))) mean(x$study_ci_covers_truth, na.rm = TRUE) else NA_real_,
      subject_cint_bias = .mc_mean(subject_error),
      subject_cint_rmse = sqrt(.mc_mean(subject_error^2)),
      subject_cint_coverage_95 = if (any(!is.na(x$subject_cint_ci_covers_truth))) mean(x$subject_cint_ci_covers_truth, na.rm = TRUE) else NA_real_,
      mean_abs_gap_per_subject = .mc_mean(abs(x$gap_per_subject)),
      maximum_abs_gap_per_subject = if (any(is.finite(x$gap_per_subject))) max(abs(x$gap_per_subject), na.rm = TRUE) else NA_real_,
      mean_fit_seconds = .mc_mean(x$fit_seconds), mean_check_seconds = .mc_mean(x$check_seconds),
      usable = sum(ok), stringsAsFactors = FALSE
    )
  }))
}

.mc_write_outputs <- function(options, manifest = NULL) {
  raw <- .mc_read_results(options)
  summary <- .mc_summarize(raw, options$replications)
  utils::write.csv(raw, file.path(options$output, "replication-summary.csv"), row.names = FALSE)
  utils::write.csv(summary, file.path(options$output, "monte-carlo-summary.csv"), row.names = FALSE)
  lines <- c("# CoTiMA2 ctsemjulia Monte Carlo validation", "",
    paste0("Generated: ", format(Sys.time(), tz = "UTC", usetz = TRUE)), "",
    paste0("- Requested replications: ", options$replications),
    paste0("- Mode: ", options$mode), paste0("- Workers: ", options$workers),
    paste0("- Cores per fit: ", options$cores_per_fit), "", "## Aggregate results", "")
  if (!nrow(summary)) {
    lines <- c(lines, "No completed scenario files were found.")
  } else {
    lines <- c(lines,
      "| Scenario | Completed | Convergence | Bias | RMSE | 95% coverage | Mean fit s |",
      "|---|---:|---:|---:|---:|---:|---:|")
    for (i in seq_len(nrow(summary))) lines <- c(lines, sprintf(
      "| %s | %d/%d | %.3f | %.4f | %.4f | %s | %.1f |",
      summary$scenario[i], summary$completed[i], summary$requested[i],
      summary$convergence_rate[i], summary$bias[i], summary$rmse[i],
      ifelse(is.na(summary$coverage_95[i]), "--", sprintf("%.3f", summary$coverage_95[i])),
      summary$mean_fit_seconds[i]))
  }
  lines <- c(lines, "", "Coverage is reported only when ctsem exposes finite 95% limits for the random-effect SD.",
    "Interrupted or failed replications can be resumed by rerunning the same command.")
  writeLines(lines, file.path(options$output, "MONTE_CARLO_REPORT.md"))
  saveRDS(list(options = options, manifest = manifest, raw = raw, summary = summary),
    file.path(options$output, "monte-carlo-results.rds"))
  summary
}

.mc_main <- function(options) {
  runner <- file.path(options$project, "tools", "ctsemjulia-validation-suite.R")
  if (!file.exists(runner)) stop("Validation runner not found: ", runner, call. = FALSE)
  dir.create(file.path(options$output, "logs"), recursive = TRUE, showWarnings = FALSE)
  detected <- parallel::detectCores(logical = FALSE) %||% parallel::detectCores()
  requested <- options$workers * options$cores_per_fit
  if (is.finite(detected) && requested > max(1L, detected - 1L)) stop(
    "Requested ", requested, " fit cores, but only ", detected,
    " physical cores were detected. Leave at least one core free.", call. = FALSE)
  replications <- seq.int(options$replication_start,
    length.out = options$replications)
  manifest <- list(created = Sys.time(), options = options, replications = replications,
    detected_physical_cores = detected)
  saveRDS(manifest, file.path(options$output, "manifest.rds"))
  if (!isTRUE(options$aggregate_only)) {
    worker <- function(replication) .mc_run_replication(
      replication, options, runner, file.path(R.home("bin"), "Rscript"))
    if (options$workers == 1L) {
      runs <- lapply(replications, worker)
    } else {
      cluster <- parallel::makeCluster(options$workers)
      on.exit(parallel::stopCluster(cluster), add = TRUE)
      parallel::clusterExport(cluster,
        c("options", "runner", ".mc_run_replication", ".mc_quote_arg"), envir = environment())
      runs <- parallel::parLapply(cluster, replications, function(replication) {
        .mc_run_replication(replication, options, runner, file.path(R.home("bin"), "Rscript"))
      })
    }
    utils::write.csv(do.call(rbind, runs), file.path(options$output, "run-status.csv"), row.names = FALSE)
  }
  summary <- .mc_write_outputs(options, manifest)
  print(summary, row.names = FALSE)
  message("Monte Carlo report: ", file.path(options$output, "MONTE_CARLO_REPORT.md"))
  invisible(summary)
}

if (sys.nframe() == 0L) .mc_main(.mc_options())
