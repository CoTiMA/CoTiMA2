#!/usr/bin/env Rscript

# End-to-end validation of the internal CoTiMA2 adapter against ctsem's native
# Julia/Laplace backend. The source data contain three simulated primary
# studies with identical true effects. The validation keeps 30 subjects per
# study so that it remains a diagnostic rather than a production fit.

args <- commandArgs(trailingOnly = TRUE)
development_library <- if (length(args)) args[[1L]] else "/private/tmp/ctsem-juliaFit-lib"
julia_bindir <- if (length(args) >= 2L) args[[2L]] else Sys.getenv("JULIA_BINDIR")
project <- if (length(args) >= 3L) args[[3L]] else
  "/Users/cdormann/SynologyDrive/Synology Drive/CHRISTIAN/TEXTE/METHODEN/R/GitHub CoTiMA2"

.libPaths(c(development_library, .libPaths()))
Sys.setenv(JULIA_BINDIR = julia_bindir, JULIA_NUM_THREADS = "4")
devtools::load_all(project, quiet = TRUE)

init <- readRDS(path.expand("~/CoTiMA_Benchmark/initBayes.rds"))
setup <- ctmaFit(
  ctmaInitFit = init,
  fit = FALSE,
  cint = "auto",
  studyCoding = "dummy",
  backend = "stan",
  coresToUse = 1L,
  verbose = 0
)

# Identify the three study profiles from their two dummy columns and retain the
# first 30 distinct people in each profile. Selection is deterministic.
profile <- interaction(setup$data$TI1, setup$data$TI2, drop = TRUE)
subject_by_study <- split(setup$data$id, profile)
selected_ids <- unlist(lapply(subject_by_study, function(id) {
  head(unique(id), 30L)
}), use.names = FALSE)
setup$data <- setup$data[setup$data$id %in% selected_ids, , drop = FALSE]

cint_parameters <- unique(stats::na.omit(
  setup$ctModel$pars$param[setup$ctModel$pars$matrix == "CINT"]
))
stopifnot(length(cint_parameters) == 2L)

prepared <- .ctmaPrepareCtsemJuliaLaplace(
  setup,
  study_dummy_ti = c("TI1", "TI2"),
  subject_parameters = cint_parameters,
  study_parameters = cint_parameters,
  remove_study_fixed_effects = TRUE,
  verbose = TRUE
)

# The selected CINT rows must vary by study and their fixed dummy effects must
# no longer be active, otherwise fixed and random study effects are redundant.
cint_rows <- prepared$model$pars$matrix == "CINT"
stopifnot(
  all(prepared$model$pars$indvarying_study[cint_rows]),
  all(!prepared$model$pars[cint_rows, c("TI1_effect", "TI2_effect")]),
  length(unique(prepared$data$study)) == 3L,
  length(unique(prepared$data$id)) == 90L
)

ctsem::ctJuliaSetup(threads = 4L, force = TRUE)
fit_time <- system.time({
  fit <- .ctmaFitCtsemJuliaLaplace(prepared, cores = 4L)
})
check_time <- system.time({
  check <- ctsem::ctLaplaceCheck(fit)
})

output <- path.expand("~/CoTiMA_Benchmark/ctsemjulia_laplace_equal_effects.rds")
saveRDS(list(
  provenance = list(
    ctsem = as.character(utils::packageVersion("ctsem")),
    CoTiMA2 = as.character(utils::packageVersion("CoTiMA2")),
    R = R.version.string,
    platform = R.version$platform,
    subjects_per_study = 30L,
    true_study_heterogeneity = 0
  ),
  prepared = prepared,
  fit = fit,
  check = check,
  fit_time = fit_time,
  check_time = check_time
), output)

cat("\nCoTiMA2 Laplace validation completed.\n")
cat("Selected CINT parameters:", paste(cint_parameters, collapse = ", "), "\n")
cat("Fit elapsed seconds:", unname(fit_time[["elapsed"]]), "\n")
cat("Laplace check elapsed seconds:", unname(check_time[["elapsed"]]), "\n")
cat("Saved result:", output, "\n")
print(check)
