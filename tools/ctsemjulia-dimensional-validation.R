#!/usr/bin/env Rscript

# Sequential, resumable localization of multivariate Laplace error. Each
# configuration is fitted in a fresh R/Julia process so a failed case cannot
# corrupt the remaining comparisons.

.dim_script_directory <- function() {
  x <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(x)) return(normalizePath("tools", mustWork = TRUE))
  dirname(normalizePath(sub("^--file=", "", x[[1L]]), mustWork = TRUE))
}

.dim_options <- function(args = commandArgs(trailingOnly = TRUE)) {
  x <- list(
    project = normalizePath(getwd(), mustWork = TRUE),
    output = path.expand("~/CoTiMA_Benchmark/ctsemjulia-dimensional"),
    development_library = "/private/tmp/ctsem-juliaFit-lib",
    julia_bindir = Sys.getenv("JULIA_BINDIR"), cores = 4L,
    mode = "full", seed = 20260909L, check_nodes = 5L,
    study_effects = "quantiles", retain_fit = FALSE, resume = TRUE,
    all_four_result = ""
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
  x$check_nodes <- as.integer(x$check_nodes)
  x$mode <- match.arg(tolower(x$mode), c("quick", "full"))
  x$study_effects <- match.arg(tolower(x$study_effects), c("quantiles", "random"))
  x$retain_fit <- flag(x$retain_fit); x$resume <- flag(x$resume)
  x
}

.dim_configurations <- list(
  cint_drift = c("CINT", "DRIFT"),
  cint_diffusion = c("CINT", "DIFFUSION"),
  cint_t0var = c("CINT", "T0VAR"),
  drift_diffusion = c("DRIFT", "DIFFUSION"),
  drift_t0var = c("DRIFT", "T0VAR"),
  diffusion_t0var = c("DIFFUSION", "T0VAR"),
  cint_drift_diffusion = c("CINT", "DRIFT", "DIFFUSION"),
  all_four = c("CINT", "DRIFT", "DIFFUSION", "T0VAR")
)

.dim_reusable <- function(file, options, components) {
  if (!file.exists(file)) return(FALSE)
  x <- tryCatch(readRDS(file), error = function(e) NULL)
  !is.null(x) && identical(x$options$mode, options$mode) &&
    identical(as.integer(x$options$seed), options$seed) &&
    identical(as.integer(x$options$check_nodes), options$check_nodes) &&
    identical(x$options$study_effects, options$study_effects) &&
    identical(x$options$components, components)
}

.dim_external_all_four <- function(file, options) {
  x <- readRDS(file)
  expected <- if (options$mode == "quick")
    c(nstudy = 8L, subjects_per_study = 20L, nobs = 12L) else
    c(nstudy = 16L, subjects_per_study = 35L, nobs = 20L)
  observed <- unlist(x$dimensions[names(expected)], use.names = TRUE)
  compatible <- identical(as.integer(observed), as.integer(expected)) &&
    identical(as.integer(x$options$seed), options$seed) &&
    identical(x$options$study_effects, options$study_effects)
  if (!compatible) stop(
    "all_four_result does not match mode, seed, or study-effect generation.",
    call. = FALSE
  )
  # Results created before --check-nodes was exposed used the hard-coded value 5.
  prior_nodes <- x$options$check_nodes
  if (is.null(prior_nodes)) prior_nodes <- 5L
  if (!identical(as.integer(prior_nodes), options$check_nodes)) stop(
    "all_four_result used a different number of quadrature nodes.",
    call. = FALSE
  )
  TRUE
}

.dim_row <- function(name, components, file, exit_status = 0L) {
  if (!file.exists(file)) return(data.frame(
    configuration = name, components = paste(components, collapse = "+"),
    n_study_effects = length(components), status = "failed",
    exit_status = exit_status, converged = FALSE, fit_seconds = NA_real_,
    check_seconds = NA_real_, gap_per_subject = NA_real_,
    maximum_sd_relative_error = NA_real_, maximum_abs_correlation = NA_real_
  ))
  x <- readRDS(file)
  data.frame(
    configuration = name, components = paste(components, collapse = "+"),
    n_study_effects = length(components), status = "completed",
    exit_status = exit_status, converged = isTRUE(x$converged),
    fit_seconds = x$fit_seconds, check_seconds = x$check_seconds,
    gap_per_subject = x$gap_per_subject,
    maximum_sd_relative_error = max(x$recovery$relative_error, na.rm = TRUE),
    maximum_abs_correlation = if (nrow(x$correlations))
      max(abs(x$correlations$estimated_raw_correlation), na.rm = TRUE) else 0
  )
}

.dim_write_report <- function(table, options) {
  utils::write.csv(table, file.path(options$output, "dimensional-summary.csv"),
    row.names = FALSE)
  fmt <- function(x) ifelse(is.finite(x), formatC(x, digits = 4, format = "fg"), "--")
  lines <- c(
    "# ctsemjulia random-effect dimensionality validation", "",
    paste0("- Mode: ", options$mode),
    paste0("- Quadrature nodes: ", options$check_nodes),
    paste0("- Study effects: ", options$study_effects), "",
    "| Configuration | k | Converged | Fit s | Check s | Gap/person | Max SD error | Max |correlation| |",
    "|---|---:|---:|---:|---:|---:|---:|---:|"
  )
  for (i in seq_len(nrow(table))) lines <- c(lines, paste0(
    "| ", table$components[[i]], " | ", table$n_study_effects[[i]], " | ",
    table$converged[[i]], " | ", fmt(table$fit_seconds[[i]]), " | ",
    fmt(table$check_seconds[[i]]), " | ", fmt(table$gap_per_subject[[i]]),
    " | ", fmt(table$maximum_sd_relative_error[[i]]), " | ",
    fmt(table$maximum_abs_correlation[[i]]), " |"
  ))
  lines <- c(lines, "", paste(
    "The sequence localizes whether the Laplace-check gap increases mainly",
    "with random-effect dimension or with a particular parameter combination."
  ))
  writeLines(lines, file.path(options$output, "DIMENSIONAL_REPORT.md"))
}

.dim_run <- function(options) {
  dir.create(options$output, recursive = TRUE, showWarnings = FALSE)
  joint_script <- file.path(.dim_script_directory(), "ctsemjulia-joint-validation.R")
  rows <- vector("list", length(.dim_configurations))
  for (i in seq_along(.dim_configurations)) {
    name <- names(.dim_configurations)[[i]]
    components <- .dim_configurations[[i]]
    case_dir <- file.path(options$output, name)
    result_file <- file.path(case_dir, "joint-validation-results.rds")
    dir.create(case_dir, recursive = TRUE, showWarnings = FALSE)
    status <- 0L
    external_all_four <- name == "all_four" &&
      nzchar(options$all_four_result) && file.exists(options$all_four_result)
    if (external_all_four) {
      .dim_external_all_four(options$all_four_result, options)
      message("Reusing external all-four result: ", options$all_four_result)
      result_file <- options$all_four_result
    } else if (!(options$resume &&
        .dim_reusable(result_file, options, components))) {
      message("Running dimensional validation: ", name)
      args <- c(joint_script,
        paste0("--project=", shQuote(options$project)),
        paste0("--output=", shQuote(case_dir)),
        paste0("--development-library=", shQuote(options$development_library)),
        paste0("--julia-bindir=", shQuote(options$julia_bindir)),
        paste0("--cores=", options$cores), paste0("--mode=", options$mode),
        paste0("--seed=", options$seed),
        paste0("--study-effects=", options$study_effects),
        paste0("--components=", paste(components, collapse = ",")),
        paste0("--check-nodes=", options$check_nodes),
        paste0("--retain-fit=", tolower(options$retain_fit)))
      status <- system2(file.path(R.home("bin"), "Rscript"), args,
        stdout = file.path(case_dir, "run.log"),
        stderr = file.path(case_dir, "error.log"))
    } else message("Reusing completed dimensional validation: ", name)
    rows[[i]] <- .dim_row(name, components, result_file, status)
    .dim_write_report(do.call(rbind, rows[seq_len(i)]), options)
  }
  table <- do.call(rbind, rows)
  saveRDS(list(options = options, summary = table),
    file.path(options$output, "dimensional-validation-results.rds"))
  print(table, row.names = FALSE); invisible(table)
}

if (sys.nframe() == 0L) .dim_run(.dim_options())
