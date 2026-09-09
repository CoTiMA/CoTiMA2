#!/usr/bin/env Rscript

# Reproducible parameter-recovery validation for CoTiMA2's experimental
# backend = "ctsemjulia". The simulations use ctsem's exact raw-parameter
# transforms and therefore compare estimated random-effect SDs with known
# raw-scale targets.

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
.ctma_validation_schema <- 2L

.ctma_validation_options <- function(args = commandArgs(trailingOnly = TRUE)) {
  defaults <- list(
    project = normalizePath(getwd(), mustWork = TRUE),
    output = path.expand("~/CoTiMA_Benchmark/ctsemjulia-validation"),
    development_library = "/private/tmp/ctsem-juliaFit-lib",
    julia_bindir = Sys.getenv("JULIA_BINDIR"),
    cores = 4L,
    mode = "quick",
    scenarios = "",
    seed = 20260909L,
    retain_fits = FALSE,
    resume = TRUE
  )
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) stop(
      "Arguments must use --name=value syntax.", call. = FALSE
    )
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    name <- gsub("-", "_", parts[[1L]], fixed = TRUE)
    if (!name %in% names(defaults)) stop("Unknown option: --", parts[[1L]],
      call. = FALSE)
    defaults[[name]] <- paste(parts[-1L], collapse = "=")
  }
  defaults$cores <- as.integer(defaults$cores)
  defaults$seed <- as.integer(defaults$seed)
  defaults$retain_fits <- tolower(as.character(defaults$retain_fits)) %in%
    c("true", "t", "1", "yes")
  defaults$resume <- tolower(as.character(defaults$resume)) %in%
    c("true", "t", "1", "yes")
  defaults$mode <- match.arg(tolower(defaults$mode), c("quick", "full"))
  defaults$scenarios <- trimws(strsplit(
    as.character(defaults$scenarios), ",", fixed = TRUE
  )[[1L]])
  defaults$scenarios <- defaults$scenarios[nzchar(defaults$scenarios)]
  unknown <- setdiff(
    defaults$scenarios, c("cint", "drift", "diffusion", "t0var")
  )
  if (length(unknown)) stop(
    "Unknown validation scenario(s): ", paste(unknown, collapse = ", "),
    call. = FALSE
  )
  defaults
}

.ctma_drift_transform <- function(raw) {
  -(1e-6 + 2 * log1p(exp(-(2 * raw))))
}

.ctma_diffusion_transform <- function(raw) {
  10 * log1p(exp(2 * raw)) + 1e-10
}

.ctma_t0var_transform <- function(raw) {
  1e-10 + 5 * log1p(exp(2 * raw))
}

.ctma_transformed_normal_sd <- function(transform, location, scale) {
  if (!is.finite(scale) || scale <= 0) return(0)
  first <- stats::integrate(function(z) {
    transform(location + scale * z) * stats::dnorm(z)
  }, -8, 8, subdivisions = 1000L, rel.tol = 1e-10)$value
  second <- stats::integrate(function(z) {
    transform(location + scale * z)^2 * stats::dnorm(z)
  }, -8, 8, subdivisions = 1000L, rel.tol = 1e-10)$value
  sqrt(max(0, second - first^2))
}

.ctma_standardized_quantiles <- function(n) {
  as.numeric(scale(stats::qnorm((seq_len(n) - 0.5) / n), TRUE, TRUE))
}

.ctma_simulate_univariate <- function(
    scenario, nstudy, subjects_per_study, nobs, seed) {
  set.seed(seed)
  scenario <- match.arg(
    scenario, c("cint", "drift", "diffusion", "t0var")
  )
  truth <- list(
    raw_cint_mean = 0.01,
    raw_drift_mean = 0.35,
    raw_diffusion_mean = -1.70,
    raw_t0var_mean = -0.80,
    raw_subject_cint_sd = if (scenario %in% c("cint", "drift")) 0.04 else 0,
    raw_study_sd = switch(
      scenario, cint = 0.03, drift = 0.30, diffusion = 0.22, t0var = 0.25
    ),
    manifest_sd = 0.10
  )
  transform <- switch(
    scenario, cint = function(x) 10 * x,
    drift = .ctma_drift_transform, diffusion = .ctma_diffusion_transform,
    t0var = .ctma_t0var_transform
  )
  location <- switch(
    scenario, cint = truth$raw_cint_mean, drift = truth$raw_drift_mean,
    diffusion = truth$raw_diffusion_mean, t0var = truth$raw_t0var_mean
  )
  truth$study_sd <- .ctma_transformed_normal_sd(
    transform, location, truth$raw_study_sd
  )
  truth$subject_cint_sd <- 10 * truth$raw_subject_cint_sd
  z <- .ctma_standardized_quantiles(nstudy)
  study_raw <- truth$raw_study_sd * z
  raw_drift <- rep(truth$raw_drift_mean, nstudy)
  raw_diffusion <- rep(truth$raw_diffusion_mean, nstudy)
  raw_t0var <- rep(truth$raw_t0var_mean, nstudy)
  raw_study_cint <- rep(0, nstudy)
  if (scenario == "cint") raw_study_cint <- study_raw
  if (scenario == "drift") raw_drift <- raw_drift + study_raw
  if (scenario == "diffusion") raw_diffusion <- raw_diffusion + study_raw
  if (scenario == "t0var") raw_t0var <- raw_t0var + study_raw

  drift <- .ctma_drift_transform(raw_drift)
  diffusion <- .ctma_diffusion_transform(raw_diffusion)
  t0var <- .ctma_t0var_transform(raw_t0var)
  rows <- vector("list", nstudy * subjects_per_study)
  person <- 0L
  for (study in seq_len(nstudy)) {
    transition <- exp(drift[[study]])
    process_variance <- diffusion[[study]]^2 / (-2 * drift[[study]]) *
      (1 - transition^2)
    for (within in seq_len(subjects_per_study)) {
      person <- person + 1L
      raw_subject_cint <- stats::rnorm(1L, sd = truth$raw_subject_cint_sd)
      cint <- 10 * (
        truth$raw_cint_mean + raw_study_cint[[study]] + raw_subject_cint
      )
      intercept <- (cint / drift[[study]]) * (transition - 1)
      x <- stats::rnorm(1L, mean = 0, sd = t0var[[study]])
      y <- numeric(nobs)
      for (occasion in seq_len(nobs)) {
        if (occasion > 1L) {
          x <- transition * x + intercept +
            stats::rnorm(1L, sd = sqrt(process_variance))
        }
        y[[occasion]] <- x + stats::rnorm(1L, sd = truth$manifest_sd)
      }
      rows[[person]] <- data.frame(
        id = person, study = study, time = 0:(nobs - 1L), Y1 = y
      )
    }
  }
  list(
    data = do.call(rbind, rows), truth = truth,
    realised_study_raw = study_raw,
    transformed = data.frame(
      study = seq_len(nstudy), cint = 10 * raw_study_cint, drift = drift,
      diffusion = diffusion, t0var = t0var
    )
  )
}

.ctma_validation_model <- function() {
  ctsem::ctModel(
    type = "ct", id = "id", manifestNames = "Y1", latentNames = "eta1",
    LAMBDA = matrix(1), CINT = matrix("cint_eta1"),
    DRIFT = matrix("drift_eta1"), DIFFUSION = matrix("diff_eta1"),
    T0VAR = matrix("T0var_eta1"), MANIFESTMEANS = matrix(0)
  )
}

.ctma_scalar <- function(x, default = NA_real_) {
  value <- suppressWarnings(as.numeric(x)[1L])
  if (!length(value) || is.na(value)) default else value
}

.ctma_interval <- function(table, parameter) {
  if (is.null(table) || !parameter %in% rownames(table)) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  columns <- colnames(table)
  lower_name <- intersect(c("2.5%", "2.5%ile", "l95", "lower"), columns)[1L]
  upper_name <- intersect(c("97.5%", "97.5%ile", "u95", "upper"), columns)[1L]
  c(
    lower = if (is.na(lower_name)) NA_real_ else .ctma_scalar(table[parameter, lower_name]),
    upper = if (is.na(upper_name)) NA_real_ else .ctma_scalar(table[parameter, upper_name])
  )
}

.ctma_run_recovery_case <- function(
    scenario, design, cores, output, retain_fit = FALSE, seed = NA_integer_) {
  parameter <- switch(
    scenario, cint = "cint_eta1", drift = "drift_eta1",
    diffusion = "diff_eta1", t0var = "T0var_eta1"
  )
  model <- .ctma_validation_model()
  model$pars$indvarying <- FALSE
  setup <- structure(list(
    data = design$data[, c("id", "time", "Y1")],
    ctModel = model, argumentList = list()
  ), class = "CoTiMAFit")
  args <- CoTiMA2::ctmaBackendDefaults("ctsemjulia")
  args$study_id <- design$data$study
  args$study_parameters <- parameter
  if (design$truth$subject_cint_sd > 0) {
    args$subject_parameters <- "cint_eta1"
  }
  fit_time <- system.time({
    fit <- CoTiMA2:::.ctmaFitCtsemJuliaBackend(
      setup = setup, backend_args = args, cores = cores,
      cotima_call = match.call(), verbose = 0
    )
  })
  check_time <- system.time({
    check <- ctsem::ctLaplaceCheck(fit, nodes = 5L)
  })
  fit_summary <- summary(fit)
  study_table <- fit_summary[["popsd.study"]] %||% fit_summary[["popsd"]]
  study_estimate <- if (!is.null(study_table) &&
      parameter %in% rownames(study_table)) {
    .ctma_scalar(study_table[parameter, "mean"])
  } else NA_real_
  study_interval <- .ctma_interval(study_table, parameter)
  subject_table <- fit_summary[["popsd.id"]]
  subject_estimate <- if (design$truth$subject_cint_sd > 0 &&
      !is.null(subject_table) && "cint_eta1" %in% rownames(subject_table)) {
    .ctma_scalar(subject_table["cint_eta1", "mean"])
  } else NA_real_
  subject_interval <- .ctma_interval(subject_table, "cint_eta1")
  result <- list(
    validation_schema = .ctma_validation_schema,
    scenario = scenario,
    seed = seed,
    status = "completed",
    nstudy = length(unique(design$data$study)),
    nsubjects = length(unique(design$data$id)),
    nobs = nrow(design$data),
    converged = isTRUE(fit$estimate$converged),
    gradient_norm = .ctma_scalar(fit$estimate$gradient_norm),
    loglik = .ctma_scalar(fit_summary$loglik),
    fit_seconds = unname(fit_time[["elapsed"]]),
    check_seconds = unname(check_time[["elapsed"]]),
    laplace_gap = .ctma_scalar(check$gap),
    gap_per_subject = .ctma_scalar(check$gap_per_subject),
    true_raw_study_sd = design$truth$raw_study_sd,
    true_study_sd = design$truth$study_sd,
    estimated_study_sd = study_estimate,
    study_ci_lower = unname(study_interval[["lower"]]),
    study_ci_upper = unname(study_interval[["upper"]]),
    study_ci_covers_truth = if (all(is.finite(study_interval))) {
      design$truth$study_sd >= study_interval[["lower"]] &&
        design$truth$study_sd <= study_interval[["upper"]]
    } else NA,
    study_sd_relative_error = abs(study_estimate - design$truth$study_sd) /
      design$truth$study_sd,
    true_subject_cint_sd = design$truth$subject_cint_sd,
    estimated_subject_cint_sd = subject_estimate,
    subject_cint_ci_lower = unname(subject_interval[["lower"]]),
    subject_cint_ci_upper = unname(subject_interval[["upper"]]),
    subject_cint_ci_covers_truth = if (
        design$truth$subject_cint_sd > 0 && all(is.finite(subject_interval))) {
      design$truth$subject_cint_sd >= subject_interval[["lower"]] &&
        design$truth$subject_cint_sd <= subject_interval[["upper"]]
    } else NA,
    fit = if (retain_fit) fit else NULL,
    check = check
  )
  saveRDS(result, file.path(output, paste0(scenario, ".rds")))
  result
}

.ctma_failed_case <- function(scenario, condition) {
  list(
    validation_schema = .ctma_validation_schema,
    scenario = scenario, status = "failed",
    error = conditionMessage(condition), converged = FALSE,
    fit_seconds = NA_real_, check_seconds = NA_real_,
    gap_per_subject = NA_real_, true_study_sd = NA_real_,
    estimated_study_sd = NA_real_, study_sd_relative_error = NA_real_,
    true_subject_cint_sd = NA_real_, estimated_subject_cint_sd = NA_real_
  )
}

.ctma_validation_table <- function(results) {
  fields <- c(
    "scenario", "status", "nstudy", "nsubjects", "nobs", "converged",
    "gradient_norm", "loglik", "fit_seconds", "check_seconds",
    "laplace_gap", "gap_per_subject", "true_study_sd",
    "true_raw_study_sd",
    "estimated_study_sd", "study_sd_relative_error",
    "study_ci_lower", "study_ci_upper", "study_ci_covers_truth",
    "true_subject_cint_sd", "estimated_subject_cint_sd",
    "subject_cint_ci_lower", "subject_cint_ci_upper",
    "subject_cint_ci_covers_truth", "seed", "error"
  )
  rows <- lapply(results, function(x) {
    out <- setNames(rep(list(NA), length(fields)), fields)
    for (name in intersect(fields, names(x))) out[[name]] <- x[[name]]
    as.data.frame(out, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

.ctma_decision <- function(row) {
  if (!identical(row$status, "completed")) return("FAIL")
  if (!isTRUE(row$converged)) return("FAIL")
  if (!is.finite(row$estimated_study_sd)) return("FAIL")
  if (is.finite(row$gap_per_subject) && abs(row$gap_per_subject) > 0.01) {
    return("WARN")
  }
  if (is.finite(row$study_sd_relative_error) &&
      row$study_sd_relative_error > 0.50) return("WARN")
  "PASS"
}

.ctma_write_report <- function(table, provenance, output) {
  decisions <- vapply(seq_len(nrow(table)), function(i) {
    .ctma_decision(as.list(table[i, , drop = FALSE]))
  }, character(1L))
  table$decision <- decisions
  utils::write.csv(table, file.path(output, "summary.csv"), row.names = FALSE)
  fmt <- function(x, digits = 4L) {
    ifelse(is.na(x), "--", formatC(x, digits = digits, format = "fg"))
  }
  lines <- c(
    "# CoTiMA2 ctsemjulia parameter-recovery validation",
    "",
    paste0("Generated: ", provenance$timestamp),
    "",
    paste0("- CoTiMA2: ", provenance$CoTiMA2),
    paste0("- ctsem: ", provenance$ctsem),
    paste0("- R: ", provenance$R),
    paste0("- Platform: ", provenance$platform),
    paste0("- Julia threads requested: ", provenance$cores),
    paste0("- Validation mode: ", provenance$mode),
    "",
    "## Results",
    "",
    paste(
      "| Scenario | Decision | Converged | Fit s | Gap/person |",
      "True raw SD | True reported-scale SD | Estimated SD | Relative error |"
    ),
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|"
  )
  for (i in seq_len(nrow(table))) {
    lines <- c(lines, paste0(
      "| ", table$scenario[[i]], " | ", table$decision[[i]], " | ",
      table$converged[[i]] %||% FALSE, " | ", fmt(table$fit_seconds[[i]], 5),
      " | ", fmt(table$gap_per_subject[[i]], 4), " | ",
      fmt(table$true_raw_study_sd[[i]], 4), " | ",
      fmt(table$true_study_sd[[i]], 4), " | ",
      fmt(table$estimated_study_sd[[i]], 4), " | ",
      fmt(table$study_sd_relative_error[[i]], 4), " |"
    ))
  }
  lines <- c(
    lines, "", "## Decision rules", "",
    "- FAIL: execution failed or the optimizer did not converge.",
    paste(
      "- WARN: absolute Laplace/quadrature gap per subject exceeds 0.01,",
      "or the single-run study-SD recovery error exceeds 50%."
    ),
    "- PASS: neither condition occurred.",
    "",
    paste(
      "These thresholds are engineering diagnostics, not Monte Carlo coverage",
      "criteria. A single successful data set does not establish unbiasedness."
    ),
    "Replicated full-mode runs are required before promoting the backend from experimental status.",
    "",
    "## Recommended interpretation", "",
    paste(
      "CINT tests an identity-transformed location effect. DRIFT, DIFFUSION,",
      "and T0VAR test nonlinear ctsem transformations. Random-effect SDs are",
      "generated on ctsem's unconstrained raw scale. summary(fit)$popsd",
      "reports the induced SD after the parameter transform, so recovery is",
      "evaluated against the corresponding transformed population SD."
    ),
    "",
    paste(
      "A fixed-effects Stan or CoTiMA2-Julia result is not used as a gold",
      "standard because it targets a different likelihood and estimand.",
      "Cross-backend timing is valid only for equivalent hierarchical models."
    )
  )
  writeLines(lines, file.path(output, "REPORT.md"))
  table
}

.ctma_run_validation_suite <- function(options) {
  if (dir.exists(options$development_library)) {
    .libPaths(c(options$development_library, .libPaths()))
  }
  if (nzchar(options$julia_bindir)) Sys.setenv(JULIA_BINDIR = options$julia_bindir)
  Sys.setenv(JULIA_NUM_THREADS = as.character(options$cores))
  dir.create(options$output, recursive = TRUE, showWarnings = FALSE)
  devtools::load_all(options$project, quiet = TRUE)
  status <- CoTiMA2:::.ctmaCtsemJuliaLaplaceStatus()
  if (!isTRUE(status$available)) stop(status$message, call. = FALSE)
  ctsem::ctJuliaSetup(threads = options$cores, force = TRUE)

  scenarios <- if (length(options$scenarios)) options$scenarios else if (
      options$mode == "quick") {
    c("cint", "drift")
  } else c("cint", "drift", "diffusion", "t0var")
  dimensions <- if (options$mode == "quick") {
    list(nstudy = 8L, subjects_per_study = 20L, nobs = 12L)
  } else list(nstudy = 16L, subjects_per_study = 35L, nobs = 20L)
  results <- vector("list", length(scenarios))
  names(results) <- scenarios
  for (i in seq_along(scenarios)) {
    scenario <- scenarios[[i]]
    scenario_seed <- options$seed + i
    scenario_file <- file.path(options$output, paste0(scenario, ".rds"))
    if (isTRUE(options$resume) && file.exists(scenario_file)) {
      prior <- tryCatch(readRDS(scenario_file), error = function(e) NULL)
      reusable <- !is.null(prior) && identical(prior$status, "completed") &&
        identical(prior$validation_schema, .ctma_validation_schema) &&
        identical(prior$scenario, scenario) &&
        identical(as.integer(prior$seed), as.integer(scenario_seed)) &&
        identical(as.integer(prior$nstudy), as.integer(dimensions$nstudy)) &&
        identical(as.integer(prior$nsubjects),
          as.integer(dimensions$nstudy * dimensions$subjects_per_study)) &&
        identical(as.integer(prior$nobs), as.integer(
          dimensions$nstudy * dimensions$subjects_per_study * dimensions$nobs))
      if (isTRUE(reusable)) {
        message("Reusing completed ctsemjulia recovery scenario: ", scenario)
        results[[i]] <- prior
        next
      }
    }
    message("Running ctsemjulia recovery scenario: ", scenario)
    design <- .ctma_simulate_univariate(
      scenario = scenario, nstudy = dimensions$nstudy,
      subjects_per_study = dimensions$subjects_per_study,
      nobs = dimensions$nobs, seed = scenario_seed
    )
    results[[i]] <- tryCatch(
      .ctma_run_recovery_case(
        scenario, design, options$cores, options$output,
        options$retain_fits, seed = scenario_seed
      ),
      error = function(e) .ctma_failed_case(scenario, e)
    )
    if (!identical(results[[i]]$status, "completed")) {
      results[[i]]$seed <- scenario_seed
      saveRDS(results[[i]], scenario_file)
    }
  }
  provenance <- list(
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
    CoTiMA2 = as.character(utils::packageVersion("CoTiMA2")),
    ctsem = as.character(utils::packageVersion("ctsem")),
    R = R.version.string, platform = R.version$platform,
    cores = options$cores, mode = options$mode, seed = options$seed,
    ctsemjulia_status = status
  )
  table <- .ctma_validation_table(results)
  table <- .ctma_write_report(table, provenance, options$output)
  saveRDS(
    list(provenance = provenance, options = options, summary = table,
      results = results),
    file.path(options$output, "validation-results.rds")
  )
  print(table, row.names = FALSE)
  message("Validation report: ", file.path(options$output, "REPORT.md"))
  invisible(table)
}

if (sys.nframe() == 0L) {
  .ctma_run_validation_suite(.ctma_validation_options())
}
