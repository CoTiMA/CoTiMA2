#!/usr/bin/env Rscript

# Nonlinear validation of the internal CoTiMA2 -> ctsem Julia/Laplace adapter.
# Subjects vary in CINT and studies vary in DRIFT. Unlike identity-transformed
# CINT effects, DRIFT effects make Laplace a genuine approximation.

args <- commandArgs(trailingOnly = TRUE)
development_library <- if (length(args)) args[[1L]] else "/private/tmp/ctsem-juliaFit-lib"
julia_bindir <- if (length(args) >= 2L) args[[2L]] else Sys.getenv("JULIA_BINDIR")
project <- if (length(args) >= 3L) args[[3L]] else
  "/Users/cdormann/SynologyDrive/Synology Drive/CHRISTIAN/TEXTE/METHODEN/R/GitHub CoTiMA2"

.libPaths(c(development_library, .libPaths()))
Sys.setenv(JULIA_BINDIR = julia_bindir, JULIA_NUM_THREADS = "4")
devtools::load_all(project, quiet = TRUE)

set.seed(20260909L)
nstudy <- 12L
subjects_per_study <- 30L
nobs <- 10L
true_drift <- seq(-0.25, -0.85, length.out = nstudy)
population_cint <- 0.20
subject_cint_sd <- 0.20
diffusion <- 0.45
manifest_sd <- 0.15

rows <- vector("list", nstudy * subjects_per_study)
truth <- vector("list", length(rows))
id <- 0L
for (study in seq_len(nstudy)) {
  for (person in seq_len(subjects_per_study)) {
    id <- id + 1L
    cint <- population_cint + rnorm(1L, sd = subject_cint_sd)
    x <- rnorm(1L, mean = -cint / true_drift[study], sd = 0.40)
    y <- numeric(nobs)
    for (occasion in seq_len(nobs)) {
      if (occasion > 1L) {
        transition <- exp(true_drift[study])
        intercept <- (cint / true_drift[study]) * (transition - 1)
        process_variance <- diffusion^2 / (-2 * true_drift[study]) *
          (1 - transition^2)
        x <- transition * x + intercept + rnorm(1L, sd = sqrt(process_variance))
      }
      y[occasion] <- x + rnorm(1L, sd = manifest_sd)
    }
    rows[[id]] <- data.frame(
      id = id,
      study = study,
      time = 0:(nobs - 1L),
      Y1 = y
    )
    truth[[id]] <- data.frame(id = id, study = study, cint = cint)
  }
}
dat <- do.call(rbind, rows)
truth <- do.call(rbind, truth)

model <- ctsem::ctModel(
  type = "ct",
  id = "id",
  manifestNames = "Y1",
  latentNames = "eta1",
  LAMBDA = matrix(1),
  CINT = matrix("cint_eta1"),
  DRIFT = matrix("drift_eta1"),
  DIFFUSION = matrix("diff_eta1"),
  MANIFESTMEANS = matrix(0)
)
model$pars$indvarying <- FALSE

cotima_like <- structure(list(
  data = dat[, c("id", "time", "Y1")],
  ctModel = model,
  argumentList = list()
), class = "CoTiMAFit")

prepared <- .ctmaPrepareCtsemJuliaLaplace(
  cotima_like,
  study_id = dat$study,
  subject_parameters = "cint_eta1",
  study_parameters = "drift_eta1",
  remove_study_fixed_effects = TRUE,
  verbose = FALSE
)

stopifnot(
  sum(prepared$model$pars$indvarying) == 1L,
  sum(prepared$model$pars$indvarying_study) == 1L,
  length(unique(prepared$data$study)) == nstudy
)

ctsem::ctJuliaSetup(threads = 4L, force = TRUE)
fit_time <- system.time({
  fit <- .ctmaFitCtsemJuliaLaplace(prepared, cores = 4L)
})
check_time <- system.time({
  check <- ctsem::ctLaplaceCheck(fit, nodes = 5L)
})

fit_summary <- summary(fit)
output <- path.expand("~/CoTiMA_Benchmark/ctsemjulia_laplace_drift.rds")
saveRDS(list(
  provenance = list(
    ctsem = as.character(utils::packageVersion("ctsem")),
    CoTiMA2 = as.character(utils::packageVersion("CoTiMA2")),
    R = R.version.string,
    platform = R.version$platform
  ),
  design = list(
    nstudy = nstudy,
    subjects_per_study = subjects_per_study,
    nobs = nobs,
    true_drift = true_drift,
    population_cint = population_cint,
    subject_cint_sd = subject_cint_sd,
    diffusion = diffusion,
    manifest_sd = manifest_sd
  ),
  truth = truth,
  prepared = prepared,
  fit = fit,
  summary = fit_summary,
  check = check,
  fit_time = fit_time,
  check_time = check_time
), output)

cat("\nNonlinear DRIFT Laplace validation completed.\n")
cat("Fit elapsed seconds:", unname(fit_time[["elapsed"]]), "\n")
cat("Laplace check elapsed seconds:", unname(check_time[["elapsed"]]), "\n")
cat("Converged:", isTRUE(fit$estimate$converged), "\n")
cat("Gradient norm:", fit$estimate$gradient_norm, "\n")
cat("Saved result:", output, "\n")
cat("\nStudy-level random-effect SD:\n")
print(fit_summary[["popsd.study"]])
print(check)
