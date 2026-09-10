#!/usr/bin/env Rscript

# Refit the selected stable basin from the standardized multistart diagnostic,
# retain Hessian uncertainty, and perform adaptive-quadrature validation.

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(name, default) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (!length(hit)) default else sub(paste0("^--", name, "="), "", hit[[1L]])
}
project <- normalizePath(getarg("project", getwd()), mustWork = TRUE)
output <- path.expand(getarg("output",
  "~/CoTiMA_Benchmark/ctsemjulia-standardized-refinement"))
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
set.seed(seed + 6L) # configuration plain_005 in the multistart controller
tm <- system.time(fit <- ctsem::ctFit(
  datalong = design$data, model = model, backend = "julia",
  intoverpop = "laplace", optimize = TRUE, priors = FALSE, cores = cores,
  optimcontrol = list(estonly = FALSE, maxiter = 1000L,
    carefulfit = FALSE, initsd = .05), verbose = 0
))
check_time <- system.time(check <- ctsem::ctLaplaceCheck(fit, nodes = 5L))
sm <- summary(fit)
result <- list(
  seed = seed, fit_seconds = unname(tm[["elapsed"]]),
  check_seconds = unname(check_time[["elapsed"]]),
  converged = isTRUE(fit$estimate$converged),
  gradient_norm = as.numeric(fit$estimate$gradient_norm)[1L],
  loglik = as.numeric(sm$loglik)[1L], popmeans = sm$popmeans,
  popsd = sm$popsd, rawpopcorr = sm$rawpopcorr,
  laplace_gap = as.numeric(check$gap)[1L],
  gap_per_subject = as.numeric(check$gap_per_subject)[1L],
  fit = fit, check = check
)
saveRDS(result, file.path(output, "refinement-results.rds"))
capture.output(
  print(result[c("seed", "fit_seconds", "check_seconds", "converged",
    "gradient_norm", "loglik", "popmeans", "popsd", "rawpopcorr",
    "laplace_gap", "gap_per_subject")]),
  file = file.path(output, "RESULTS.txt")
)
print(result[c("seed", "fit_seconds", "check_seconds", "converged",
  "gradient_norm", "loglik", "popmeans", "popsd", "rawpopcorr",
  "laplace_gap", "gap_per_subject")])
