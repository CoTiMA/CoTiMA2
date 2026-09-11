#!/usr/bin/env Rscript

# Focused smoke test for parameter-dependent T0VAR expressions on ctsem's
# native Julia/Laplace route. Run against an isolated ctsem juliaFit library.

args <- commandArgs(trailingOnly = TRUE)
development_library <- if (length(args)) args[[1L]] else {
  "/private/tmp/ctsem-juliaFit-lib"
}
.libPaths(c(development_library, .libPaths()))

suppressPackageStartupMessages(library(ctsem))
stopifnot(utils::packageVersion("ctsem") >= "3.12.0")

set.seed(20260911)
n_subjects <- 12L
n_obs <- 5L
dat <- do.call(rbind, lapply(seq_len(n_subjects), function(id) {
  person_mean <- stats::rnorm(1L, 0, 0.35)
  data.frame(
    id = id,
    time = 0:(n_obs - 1L),
    Y1 = person_mean + stats::rnorm(n_obs, sd = 0.7),
    Y2 = person_mean + stats::rnorm(n_obs, sd = 0.8)
  )
}))

model <- suppressWarnings(suppressMessages(ctModel(
  type = "ct",
  manifestNames = c("Y1", "Y2"),
  latentNames = c("eta1", "eta2"),
  LAMBDA = diag(2),
  DRIFT = diag(-0.5, 2),
  DIFFUSION = diag(0.6, 2),
  MANIFESTVAR = diag(0.2, 2),
  MANIFESTMEANS = matrix(c("mmean1", "mmean2"), 2, 1),
  T0MEANS = matrix(0, 2, 1),
  CINT = matrix(0, 2, 1),
  PARS = c("aa", "bb"),
  T0VAR = matrix(
    c("log1p_exp(aa + bb)", 0, "tanh(aa - bb)", "log1p_exp(aa)"),
    2, 2, byrow = TRUE
  )
)))

model$pars$indvarying <- FALSE
model$pars$indvarying[
  !is.na(model$pars$param) & model$pars$param %in% c("mmean1", "mmean2")
] <- TRUE

prepared <- suppressMessages(ctFit(
  datalong = dat,
  model = model,
  backend = "julia",
  intoverpop = "laplace",
  fit = FALSE,
  cores = 1L
))

t0_rows <- prepared$parameter_table$matrix == "T0VAR"
t0_table <- prepared$parameter_table[t0_rows, c(
  "matrix", "row", "col", "parnumber", "transform"
)]
stopifnot(any(grepl("param\\[", t0_table$transform)))

timing <- system.time({
  fitted <- suppressMessages(ctFit(
    datalong = dat,
    model = model,
    backend = "julia",
    intoverpop = "laplace",
    cores = 1L,
    optimcontrol = list(estonly = TRUE)
  ))
})

result <- list(
  ctsem_version = as.character(utils::packageVersion("ctsem")),
  prepared = TRUE,
  t0var_parameter_table = t0_table,
  finite_loglik = is.finite(fitted$estimate$loglik),
  converged = isTRUE(fitted$estimate$converged),
  gradient_norm = fitted$estimate$gradient_norm,
  elapsed_seconds = unname(timing[["elapsed"]])
)
print(result)
stopifnot(result$finite_loglik)
