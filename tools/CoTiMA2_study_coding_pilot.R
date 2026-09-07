# CoTiMA2 short study-coding diagnostic
#
# Purpose
# -------
# Compare the numerical geometry of dummy and weighted study coding without
# starting another full Bayesian fit.  The setup-only comparison takes only a
# few seconds.  An optional short Julia ML pilot can be enabled below.
#
# This file is a development diagnostic and is excluded from package builds
# because the tools/ directory is listed in .Rbuildignore.

benchmark_dir <- path.expand("~/CoTiMA_Benchmark")
init_file <- file.path(benchmark_dir, "initBayes.rds")

# The setup-only diagnostic is the safe default. To add the short Julia
# optimization comparison, run the script with the environment variable
# COTIMA_RUN_JULIA_PILOTS=true.
run_julia_pilots <- identical(
  tolower(Sys.getenv("COTIMA_RUN_JULIA_PILOTS", unset = "false")),
  "true"
)
julia_pilot_iterations <- 50L

if (!file.exists(init_file)) {
  stop("Cannot find the benchmark input: ", init_file)
}
# When the script is run from a CoTiMA2 source checkout, load that checkout so
# the diagnostic cannot accidentally use an older installed package version.
script_file <- tryCatch(
  normalizePath(sys.frame(1L)$ofile),
  error = function(e) NA_character_
)
project_root <- if (!is.na(script_file)) {
  dirname(dirname(script_file))
} else {
  normalizePath(getwd())
}
description_file <- file.path(project_root, "DESCRIPTION")
is_cotima2_checkout <- file.exists(description_file) &&
  identical(
    unname(read.dcf(description_file, fields = "Package")[1L, 1L]),
    "CoTiMA2"
  )

if (is_cotima2_checkout && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(project_root, quiet = TRUE)
} else if (!requireNamespace("CoTiMA2", quietly = TRUE)) {
  stop("Install CoTiMA2, or run this script from its source checkout.")
}

init_bayes <- readRDS(init_file)
varying_families <- c("DIFFUSION", "T0VAR")

make_setup <- function(coding) {
  CoTiMA2::ctmaFit(
    ctmaInitFit = init_bayes,
    fit = FALSE,
    studyCoding = coding,
    studyVarying = varying_families,
    verbose = 0
  )
}

study_design_diagnostics <- function(setup, coding) {
  dat <- setup$data
  study_rows <- !duplicated(dat$id)
  ti_names <- grep("^TI[0-9]+$", names(dat), value = TRUE)
  if (!length(ti_names)) {
    stop("No generated study-indicator columns were found.")
  }

  contrasts <- as.matrix(dat[study_rows, ti_names, drop = FALSE])
  design <- cbind(population = 1, contrasts)
  gram <- crossprod(design) / nrow(design)
  intercept_correlations <- gram[1L, -1L] /
    sqrt(gram[1L, 1L] * diag(gram)[-1L])

  list(
    coding = coding,
    n_subjects = nrow(design),
    unique_study_codes = unique(as.data.frame(contrasts)),
    normalized_crossproduct = gram,
    condition_number = kappa(design, exact = TRUE),
    maximum_population_contrast_correlation =
      max(abs(intercept_correlations)),
    contrast_correlation = stats::cor(contrasts)
  )
}

dummy_setup <- make_setup("dummy")
weighted_setup <- make_setup("weighted")

design_results <- list(
  dummy = study_design_diagnostics(dummy_setup, "dummy"),
  weighted = study_design_diagnostics(weighted_setup, "weighted")
)

cat("\nStudy-coding design comparison\n")
cat("==============================\n")
for (coding in names(design_results)) {
  result <- design_results[[coding]]
  cat("\n", coding, " coding\n", sep = "")
  cat("Condition number: ", format(result$condition_number, digits = 6L), "\n",
      sep = "")
  cat("Maximum population/contrast correlation: ",
      format(result$maximum_population_contrast_correlation, digits = 6L),
      "\n", sep = "")
  cat("Study codes:\n")
  print(result$unique_study_codes)
  cat("Normalized cross-product matrix:\n")
  print(result$normalized_crossproduct)
}

cat("\nInterpretation\n")
cat("--------------\n")
cat(paste(
  "A smaller condition number indicates a better-conditioned study design.",
  "A population/contrast correlation of zero means that the population",
  "coefficient is orthogonal to the study contrasts at the design level.",
  "This diagnostic concerns parameterization only; it does not establish",
  "Bayesian convergence or repair weak model identification.\n"
))

julia_results <- NULL
equivalence_result <- NULL
if (isTRUE(run_julia_pilots)) {
  run_julia_pilot <- function(coding) {
    elapsed <- system.time({
      fit <- CoTiMA2::ctmaFit(
        ctmaInitFit = init_bayes,
        backend = "julia",
        studyCoding = coding,
        studyVarying = varying_families,
        verbose = 0,
        backendArgs = list(
          maxiters = julia_pilot_iterations,
          carefulFit = FALSE,
          uncertainty = FALSE,
          transformed = FALSE,
          validate = FALSE,
          multistart = 1L
        )
      )
    })

    list(
      coding = coding,
      elapsed_seconds = unname(elapsed[["elapsed"]]),
      objective = tail(fit$optimization$stage_objectives, 1L),
      maximum_absolute_gradient =
        fit$optimization$maximum_absolute_gradient,
      iterations = fit$optimization$iterations,
      converged = isTRUE(fit$optimization$converged),
      fit = fit
    )
  }

  set.seed(20260907L)
  julia_results <- list(
    dummy = run_julia_pilot("dummy"),
    weighted = run_julia_pilot("weighted")
  )

  cat("\nShort Julia optimization pilots\n")
  cat("===============================\n")
  print(do.call(rbind, lapply(julia_results, function(x) {
    data.frame(
      coding = x$coding,
      elapsed_seconds = x$elapsed_seconds,
      objective = x$objective,
      maximum_absolute_gradient = x$maximum_absolute_gradient,
      iterations = x$iterations,
      converged = x$converged
    )
  })))

  # A fair coding comparison must use equivalent starting values. Raw values
  # from dummy coding cannot be passed unchanged to weighted coding because
  # they describe a different point in the model. Transform the nearly fitted
  # dummy coordinates so that every study-specific raw parameter is preserved.
  transform_coding <- function(raw, from_setup, to_setup) {
    pars <- from_setup$ctModel$pars
    parameter_rows <- which(!is.na(pars$param))
    base_names <- pars$param[parameter_rows]
    effect_columns <- grep("^TI[0-9]+_effect$", names(pars), value = TRUE)
    active <- which(rowSums(pars[parameter_rows, effect_columns,
      drop = FALSE]) > 0L)

    from_subjects <- !duplicated(from_setup$data$id)
    to_subjects <- !duplicated(to_setup$data$id)
    ti_names <- sub("_effect$", "", effect_columns)
    from_design <- unique(cbind(
      population = 1,
      as.matrix(from_setup$data[from_subjects, ti_names, drop = FALSE])
    ))
    to_design <- unique(cbind(
      population = 1,
      as.matrix(to_setup$data[to_subjects, ti_names, drop = FALSE])
    ))
    if (nrow(from_design) != ncol(from_design) ||
        !identical(dim(from_design), dim(to_design))) {
      stop("The study designs are not square, equivalent contrast systems.")
    }

    n_base <- length(base_names)
    n_contrasts <- length(effect_columns)
    expected <- n_base + length(active) * n_contrasts
    if (length(raw) != expected) {
      stop("Unexpected raw-vector length in the coding transformation.")
    }

    transformed <- raw
    effect_values <- matrix(
      raw[seq.int(n_base + 1L, expected)],
      nrow = length(active),
      ncol = n_contrasts,
      byrow = TRUE
    )
    for (i in seq_along(active)) {
      parameter_index <- active[[i]]
      from_coordinates <- c(raw[[parameter_index]], effect_values[i, ])
      study_values <- as.numeric(from_design %*% from_coordinates)
      to_coordinates <- as.numeric(solve(to_design, study_values))
      transformed[[parameter_index]] <- to_coordinates[[1L]]
      effect_values[i, ] <- to_coordinates[-1L]
    }
    transformed[seq.int(n_base + 1L, expected)] <- as.vector(t(effect_values))
    transformed
  }

  equivalent_weighted_inits <- transform_coding(
    julia_results$dummy$fit$coefficients$raw,
    dummy_setup,
    weighted_setup
  )
  equivalence_elapsed <- system.time({
    equivalence_fit <- CoTiMA2::ctmaFit(
      ctmaInitFit = init_bayes,
      backend = "julia",
      studyCoding = "weighted",
      studyVarying = varying_families,
      verbose = 0,
      backendArgs = list(
        inits = equivalent_weighted_inits,
        maxiters = 5L,
        optimizer_initial_step = 0.001,
        carefulFit = FALSE,
        uncertainty = FALSE,
        transformed = FALSE,
        validate = FALSE,
        multistart = 1L
      )
    )
  })
  equivalence_result <- list(
    elapsed_seconds = unname(equivalence_elapsed[["elapsed"]]),
    transformed_inits = equivalent_weighted_inits,
    objective = tail(equivalence_fit$optimization$stage_objectives, 1L),
    maximum_absolute_gradient =
      equivalence_fit$optimization$maximum_absolute_gradient,
    iterations = equivalence_fit$optimization$iterations,
    converged = isTRUE(equivalence_fit$optimization$converged),
    fit = equivalence_fit
  )
  cat("\nWeighted coding from equivalent dummy coordinates\n")
  cat("=================================================\n")
  print(equivalence_result[c(
    "elapsed_seconds", "objective", "maximum_absolute_gradient",
    "iterations", "converged"
  )])
}

output <- list(
  created = Sys.time(),
  package_version = as.character(utils::packageVersion("CoTiMA2")),
  input_file = normalizePath(init_file),
  varying_families = varying_families,
  design = design_results,
  julia = julia_results,
  equivalence = equivalence_result
)

output_file <- file.path(
  benchmark_dir,
  paste0(
    "CoTiMA2_study_coding_pilot_",
    format(Sys.Date(), "%Y-%m-%d"),
    ".rds"
  )
)
saveRDS(output, output_file)
cat("\nSaved diagnostic results to:\n", normalizePath(output_file), "\n")
