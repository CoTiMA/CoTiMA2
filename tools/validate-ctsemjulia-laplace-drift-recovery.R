#!/usr/bin/env Rscript

# Parameter-recovery validation for nonlinear study-level DRIFT effects.
# DRIFT is simulated on ctsem's own unconstrained raw scale and transformed by
# the exact model transform, making the target population SD directly known.

args <- commandArgs(trailingOnly = TRUE)
development_library <- if (length(args)) args[[1L]] else "/private/tmp/ctsem-juliaFit-lib"
julia_bindir <- if (length(args) >= 2L) args[[2L]] else Sys.getenv("JULIA_BINDIR")
project <- if (length(args) >= 3L) args[[3L]] else
  "/Users/cdormann/SynologyDrive/Synology Drive/CHRISTIAN/TEXTE/METHODEN/R/GitHub CoTiMA2"

.libPaths(c(development_library, .libPaths()))
Sys.setenv(JULIA_BINDIR = julia_bindir, JULIA_NUM_THREADS = "4")
devtools::load_all(project, quiet = TRUE)

set.seed(20260910L)
nstudy <- 16L
subjects_per_study <- 35L
nobs <- 20L
raw_drift_mean <- 0.35
raw_drift_sd <- 0.30
population_cint <- 0.15
subject_cint_sd <- 0.50
diffusion <- 0.40
manifest_sd <- 0.12

# This is the transform printed in ctModel() for a free diagonal DRIFT cell.
drift_transform <- function(raw) {
  -(1e-6 + 2 * log1p(exp(-(2 * raw))))
}

# Use deterministic standard-normal quantiles. Their exact sample SD is then
# rescaled to one, so the realised between-study raw SD equals the target.
study_z <- stats::qnorm((seq_len(nstudy) - 0.5) / nstudy)
study_z <- as.numeric(scale(study_z, center = TRUE, scale = TRUE))
true_raw_drift <- raw_drift_mean + raw_drift_sd * study_z
true_drift <- drift_transform(true_raw_drift)
stopifnot(abs(stats::sd(true_raw_drift) - raw_drift_sd) < 1e-12)

rows <- vector("list", nstudy * subjects_per_study)
truth <- vector("list", length(rows))
id <- 0L
for (study in seq_len(nstudy)) {
  transition <- exp(true_drift[study])
  process_variance <- diffusion^2 / (-2 * true_drift[study]) *
    (1 - transition^2)
  for (person in seq_len(subjects_per_study)) {
    id <- id + 1L
    cint <- population_cint + rnorm(1L, sd = subject_cint_sd)
    equilibrium <- -cint / true_drift[study]
    x <- rnorm(1L, mean = equilibrium, sd = 0.45)
    y <- numeric(nobs)
    intercept <- (cint / true_drift[study]) * (transition - 1)
    for (occasion in seq_len(nobs)) {
      if (occasion > 1L) {
        x <- transition * x + intercept +
          rnorm(1L, sd = sqrt(process_variance))
      }
      y[occasion] <- x + rnorm(1L, sd = manifest_sd)
    }
    rows[[id]] <- data.frame(
      id = id, study = study, time = 0:(nobs - 1L), Y1 = y
    )
    truth[[id]] <- data.frame(
      id = id, study = study, cint = cint,
      raw_drift = true_raw_drift[study], drift = true_drift[study]
    )
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
  verbose = FALSE
)

ctsem::ctJuliaSetup(threads = 4L, force = TRUE)
fit_time <- system.time({
  fit <- .ctmaFitCtsemJuliaLaplace(prepared, cores = 4L)
})
check_time <- system.time({
  check <- ctsem::ctLaplaceCheck(fit, nodes = 5L)
})
correction_time <- system.time({
  corrected_fit <- ctsem::ctLaplaceCorrect(
    fit, nodes = 5L, draws = "normal", correct_estimate = TRUE
  )
})

fit_summary <- summary(fit)
subject_sd <- fit_summary[["popsd.id"]]["cint_eta1", , drop = FALSE]
study_sd <- fit_summary[["popsd.study"]]["drift_eta1", , drop = FALSE]
study_sd_estimate <- unname(study_sd[1L, "mean"])
study_sd_relative_error <- abs(study_sd_estimate - raw_drift_sd) / raw_drift_sd

output <- path.expand(
  "~/CoTiMA_Benchmark/ctsemjulia_laplace_drift_recovery.rds"
)
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
    raw_drift_mean = raw_drift_mean,
    raw_drift_sd = raw_drift_sd,
    true_raw_drift = true_raw_drift,
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
  corrected_fit = corrected_fit,
  recovery = list(
    subject_cint_sd = subject_sd,
    study_drift_sd = study_sd,
    study_drift_sd_relative_error = study_sd_relative_error
  ),
  fit_time = fit_time,
  check_time = check_time,
  correction_time = correction_time
), output)

cat("\nNonlinear DRIFT recovery validation completed.\n")
cat("Fit elapsed seconds:", unname(fit_time[["elapsed"]]), "\n")
cat("Laplace check elapsed seconds:", unname(check_time[["elapsed"]]), "\n")
cat("Laplace correction elapsed seconds:",
  unname(correction_time[["elapsed"]]), "\n")
cat("Converged:", isTRUE(fit$estimate$converged), "\n")
cat("Gradient norm:", fit$estimate$gradient_norm, "\n")
cat("True raw study DRIFT SD:", raw_drift_sd, "\n")
cat("Estimated raw study DRIFT SD:", study_sd_estimate, "\n")
cat("Relative recovery error:", study_sd_relative_error, "\n")
cat("Saved result:", output, "\n")
cat("\nSubject-level CINT SD:\n")
print(subject_sd)
cat("\nStudy-level DRIFT SD:\n")
print(study_sd)
print(check)
