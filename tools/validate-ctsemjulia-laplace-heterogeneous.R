#!/usr/bin/env Rscript

# Heterogeneity validation for the internal ctsem Julia/Laplace adapter. Nine
# pseudo-studies are formed from the three equal-effect benchmark studies and
# receive known equilibrium shifts. This tests recovery away from the zero-
# variance boundary; it is not intended as a substantive simulation study.

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

original_profile <- interaction(setup$data$TI1, setup$data$TI2, drop = TRUE)
original_by_subject <- tapply(
  original_profile, setup$data$id,
  function(x) as.character(x[[1L]])
)

# Select 90 people from each original study and divide each set into three
# groups. Interleaving the offsets across original studies reduces confounding
# between an original-study peculiarity and the imposed heterogeneity.
subject_groups <- split(as.numeric(names(original_by_subject)), original_by_subject)
selected <- lapply(subject_groups, function(ids) head(ids, 90L))
pseudo_map <- data.frame(id = numeric(), pseudo_study = integer())
counter <- 0L
for (block in seq_along(selected)) {
  ids <- selected[[block]]
  stopifnot(length(ids) == 90L)
  for (within in 1:3) {
    counter <- counter + 1L
    take <- ids[seq.int((within - 1L) * 30L + 1L, within * 30L)]
    pseudo_map <- rbind(
      pseudo_map,
      data.frame(id = take, pseudo_study = counter)
    )
  }
}

dat <- merge(setup$data, pseudo_map, by = "id", sort = FALSE)
dat <- dat[order(dat$id, dat$time), , drop = FALSE]

# Deterministic, centered shifts on the observed scale. CINT is identity-
# transformed, so its random-effect integral remains in the exact Laplace case.
shift_v1 <- seq(-0.8, 0.8, length.out = 9L)
shift_v2 <- rev(seq(-0.6, 0.6, length.out = 9L))
dat$V1 <- dat$V1 + shift_v1[dat$pseudo_study]
dat$V2 <- dat$V2 + shift_v2[dat$pseudo_study]
setup$data <- dat[, setdiff(names(dat), "pseudo_study"), drop = FALSE]

cint_parameters <- unique(stats::na.omit(
  setup$ctModel$pars$param[setup$ctModel$pars$matrix == "CINT"]
))
prepared <- .ctmaPrepareCtsemJuliaLaplace(
  setup,
  study_id = dat$pseudo_study,
  study_dummy_ti = c("TI1", "TI2"),
  subject_parameters = cint_parameters,
  study_parameters = cint_parameters,
  remove_study_fixed_effects = TRUE,
  verbose = FALSE
)

stopifnot(
  length(unique(prepared$data$study)) == 9L,
  length(unique(prepared$data$id)) == 270L
)

ctsem::ctJuliaSetup(threads = 4L, force = TRUE)
fit_time <- system.time({
  fit <- .ctmaFitCtsemJuliaLaplace(prepared, cores = 4L)
})
check_time <- system.time({
  check <- ctsem::ctLaplaceCheck(fit)
})

output <- path.expand("~/CoTiMA_Benchmark/ctsemjulia_laplace_heterogeneous.rds")
saveRDS(list(
  provenance = list(
    ctsem = as.character(utils::packageVersion("ctsem")),
    CoTiMA2 = as.character(utils::packageVersion("CoTiMA2")),
    R = R.version.string,
    platform = R.version$platform,
    subjects_per_study = 30L,
    imposed_shifts = data.frame(
      study = 1:9, V1 = shift_v1, V2 = shift_v2
    )
  ),
  prepared = prepared,
  fit = fit,
  check = check,
  fit_time = fit_time,
  check_time = check_time
), output)

cat("\nHeterogeneous CoTiMA2 Laplace validation completed.\n")
cat("Fit elapsed seconds:", unname(fit_time[["elapsed"]]), "\n")
cat("Laplace check elapsed seconds:", unname(check_time[["elapsed"]]), "\n")
cat("Saved result:", output, "\n")
print(fit$collapsedScales)
print(check)
