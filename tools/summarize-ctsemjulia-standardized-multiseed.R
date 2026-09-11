#!/usr/bin/env Rscript

# Aggregate independent runs from ctsemjulia-standardized-validation.R.
# Usage:
# Rscript tools/summarize-ctsemjulia-standardized-multiseed.R \
#   --input=/path/to/stationary-correlation-multiseed

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

.options <- function(args = commandArgs(trailingOnly = TRUE)) {
  x <- list(input = path.expand(
    "~/CoTiMA2-Benchmark/stationary-correlation-multiseed"
  ))
  for (a in args) {
    if (!startsWith(a, "--input=")) stop(
      "Use --input=/path/to/results.", call. = FALSE
    )
    x$input <- sub("^--input=", "", a)
  }
  x$input <- normalizePath(x$input, mustWork = TRUE)
  x
}

.result_rows <- function(file) {
  stored <- readRDS(file)
  seed <- as.integer(stored$options$seed %||% NA_integer_)
  rows <- lapply(stored$results, function(x) {
    if (!identical(x$status, "completed")) return(data.frame(
      seed = seed, scenario = x$scenario %||% NA_character_,
      status = x$status %||% "failed", error = x$error %||% "",
      converged = FALSE, gradient_norm = NA_real_, fit_seconds = NA_real_,
      laplace_gap = NA_real_, gap_per_subject = NA_real_,
      nweak = NA_integer_, identifiability_condition = NA_real_,
      weak_directions = NA_character_,
      parameter = NA_character_, truth = NA_real_, estimate = NA_real_,
      error_sd = NA_real_, relative_error = NA_real_,
      stringsAsFactors = FALSE
    ))
    truth <- x$truth
    estimate <- x$estimated_raw_study_sd[names(truth)]
    data.frame(
      seed = seed, scenario = x$scenario, status = x$status, error = "",
      converged = isTRUE(x$converged), gradient_norm = x$gradient_norm,
      fit_seconds = x$fit_seconds, laplace_gap = x$laplace_gap,
      gap_per_subject = x$gap_per_subject,
      nweak = x$identifiability$nweak %||% NA_integer_,
      identifiability_condition = x$identifiability$condition %||% NA_real_,
      weak_directions = paste(
        x$identifiability$parameters %||% character(), collapse = "; "
      ),
      parameter = names(truth), truth = as.numeric(truth),
      estimate = as.numeric(estimate),
      error_sd = as.numeric(estimate - truth),
      relative_error = as.numeric(abs(estimate - truth) / abs(truth)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.main <- function(opts) {
  files <- Sys.glob(file.path(
    opts$input, "seed-*", "standardized-validation-results.rds"
  ))
  if (!length(files)) stop("No validation result files found.", call. = FALSE)
  detail <- do.call(rbind, lapply(files, .result_rows))
  completed <- detail[detail$status == "completed", , drop = FALSE]
  summary_rows <- lapply(split(
    completed, interaction(completed$scenario, completed$parameter, drop = TRUE)
  ), function(x) {
    data.frame(
      scenario = x$scenario[[1L]], parameter = x$parameter[[1L]],
      runs = nrow(x),
      true_sd = x$truth[[1L]], mean_estimate = mean(x$estimate),
      median_estimate = stats::median(x$estimate),
      empirical_sd = stats::sd(x$estimate),
      mean_error = mean(x$error_sd),
      rmse = sqrt(mean(x$error_sd^2)),
      mean_absolute_relative_error = mean(x$relative_error),
      stringsAsFactors = FALSE
    )
  })
  aggregate <- do.call(rbind, summary_rows)
  rownames(aggregate) <- NULL
  run_directions <- unique(completed[c(
    "seed", "scenario", "weak_directions"
  )])
  direction_values <- unlist(strsplit(
    run_directions$weak_directions[nzchar(run_directions$weak_directions)],
    "; ", fixed = TRUE
  ), use.names = FALSE)
  direction_summary <- if (length(direction_values)) {
    counts <- sort(table(direction_values), decreasing = TRUE)
    data.frame(
      direction = names(counts), seeds_flagged = as.integer(counts),
      fraction_of_completed_seeds = as.integer(counts) /
        length(unique(completed$seed)),
      stringsAsFactors = FALSE
    )
  } else data.frame(
    direction = character(), seeds_flagged = integer(),
    fraction_of_completed_seeds = numeric(), stringsAsFactors = FALSE
  )
  write.csv(detail, file.path(opts$input, "MULTISEED_DETAIL.csv"),
    row.names = FALSE)
  write.csv(aggregate, file.path(opts$input, "MULTISEED_SUMMARY.csv"),
    row.names = FALSE)
  write.csv(direction_summary,
    file.path(opts$input, "MULTISEED_DIRECTIONS.csv"), row.names = FALSE)
  saveRDS(list(
    detail = detail, summary = aggregate, directions = direction_summary
  ),
    file.path(opts$input, "MULTISEED_SUMMARY.rds"))
  runs <- unique(completed[c("seed", "converged", "gradient_norm",
    "fit_seconds", "gap_per_subject", "nweak",
    "identifiability_condition")])
  lines <- c(
    "# Stationary-correlation multiseed validation", "",
    paste0("Completed seeds: ", nrow(runs)),
    paste0("Converged seeds: ", sum(runs$converged)),
    paste0("Maximum gradient norm: ", signif(max(runs$gradient_norm), 5L)),
    paste0("Maximum Laplace gap per subject: ",
      signif(max(abs(runs$gap_per_subject)), 5L)),
    paste0("Seeds with unidentified directions: ", sum(runs$nweak > 0L)),
    "", "| Scenario | Parameter | True SD | Mean estimate | Median | Empirical SD | RMSE | Mean absolute relative error |",
    "|---|---|---:|---:|---:|---:|---:|---:|"
  )
  for (i in seq_len(nrow(aggregate))) lines <- c(lines, sprintf(
    "| %s | %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |",
    aggregate$scenario[i], aggregate$parameter[i], aggregate$true_sd[i],
    aggregate$mean_estimate[i], aggregate$median_estimate[i],
    aggregate$empirical_sd[i], aggregate$rmse[i],
    aggregate$mean_absolute_relative_error[i]
  ))
  lines <- c(lines, "", "## Weak identifiability directions", "")
  if (!nrow(direction_summary)) {
    lines <- c(lines, "No weak directions were returned.")
  } else {
    lines <- c(lines,
      "| Direction | Seeds flagged | Fraction of completed seeds |",
      "|---|---:|---:|"
    )
    for (i in seq_len(nrow(direction_summary))) lines <- c(lines, sprintf(
      "| `%s` | %d | %.3f |", direction_summary$direction[i],
      direction_summary$seeds_flagged[i],
      direction_summary$fraction_of_completed_seeds[i]
    ))
  }
  writeLines(lines, file.path(opts$input, "MULTISEED_REPORT.md"))
  print(aggregate)
  invisible(list(
    detail = detail, summary = aggregate, directions = direction_summary
  ))
}

if (sys.nframe() == 0L) .main(.options())
