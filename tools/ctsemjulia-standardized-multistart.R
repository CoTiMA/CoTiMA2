#!/usr/bin/env Rscript

# Diagnose optimizer-basin sensitivity of the bivariate covariance-stationary
# validation model. Run the main standardized validation first.

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(name, default) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(hit)) default else sub(paste0("^--", name, "="), "", hit[[1L]])
}
project <- normalizePath(getarg("project", getwd()), mustWork = TRUE)
output <- path.expand(getarg("output",
  "~/CoTiMA_Benchmark/ctsemjulia-standardized-multistart"))
devlib <- getarg("development-library", "/private/tmp/ctsem-juliaFit-lib")
julia_bindir <- getarg("julia-bindir", Sys.getenv("JULIA_BINDIR"))
cores <- as.integer(getarg("cores", "4"))
seed <- as.integer(getarg("seed", "20260912"))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
if (dir.exists(devlib)) .libPaths(c(devlib, .libPaths()))
if (nzchar(julia_bindir)) Sys.setenv(JULIA_BINDIR = julia_bindir)
Sys.setenv(JULIA_NUM_THREADS = as.character(cores))
devtools::load_all(project, quiet = TRUE)
source(file.path(project, "tools", "ctsemjulia-standardized-validation.R"),
  local = .GlobalEnv)
ctsem::ctJuliaSetup(threads = cores, force = TRUE)

design <- .simulate("bivariate_drift", 16L, 35L, 20L, seed)
model <- .model("bivariate_drift")
configs <- data.frame(
  label = c("careful10_default", "plain_default", "careful10_001",
    "careful10_005", "careful10_020", "plain_005"),
  carefulfit = c(10, 0, 10, 10, 10, 0),
  initsd = c(NA, NA, .01, .05, .20, .05)
)

one <- function(i) {
  set.seed(seed + i)
  control <- list(estonly = TRUE, maxiter = 1000L,
    carefulfit = if (configs$carefulfit[[i]] == 0) FALSE else
      configs$carefulfit[[i]])
  if (is.finite(configs$initsd[[i]])) control$initsd <- configs$initsd[[i]]
  tm <- system.time(fit <- ctsem::ctFit(
    datalong = design$data, model = model, backend = "julia",
    intoverpop = "laplace", optimize = TRUE, priors = FALSE, cores = cores,
    optimcontrol = control, verbose = 0
  ))
  sm <- summary(fit)
  popsd <- sm$popsd
  means <- sm$popmeans
  corr <- sm$rawpopcorr
  row <- data.frame(
    label = configs$label[[i]], carefulfit = configs$carefulfit[[i]],
    initsd = configs$initsd[[i]], elapsed = unname(tm[["elapsed"]]),
    converged = isTRUE(fit$estimate$converged),
    gradient_norm = as.numeric(fit$estimate$gradient_norm)[1L],
    loglik = as.numeric(sm$loglik)[1L],
    drift1_mean = means["drift1", "mean"],
    drift2_mean = means["drift2", "mean"],
    drift1_sd = popsd["drift1", "mean"],
    drift2_sd = popsd["drift2", "mean"],
    correlation = corr["drift2__drift1", "mean"], error = NA_character_
  )
  saveRDS(list(row = row, fit = fit),
    file.path(output, paste0(configs$label[[i]], ".rds")))
  row
}

rows <- vector("list", nrow(configs))
for (i in seq_len(nrow(configs))) {
  message("Running multistart configuration: ", configs$label[[i]])
  rows[[i]] <- tryCatch(one(i), error = function(e) data.frame(
    label = configs$label[[i]], carefulfit = configs$carefulfit[[i]],
    initsd = configs$initsd[[i]], elapsed = NA, converged = FALSE,
    gradient_norm = NA, loglik = NA, drift1_mean = NA, drift2_mean = NA,
    drift1_sd = NA, drift2_sd = NA, correlation = NA,
    error = conditionMessage(e)
  ))
  result <- do.call(rbind, rows[seq_len(i)])
  saveRDS(result, file.path(output, "multistart-summary.rds"))
  utils::write.csv(result, file.path(output, "multistart-summary.csv"),
    row.names = FALSE)
  print(result, row.names = FALSE)
}
