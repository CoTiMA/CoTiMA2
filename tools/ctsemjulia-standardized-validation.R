#!/usr/bin/env Rscript

# Validate the native ctsem Julia/Laplace backend under the covariance-
# stationary parameterization appropriate for idealized standardized latent
# processes.  Unlike the generic engineering suite, DIFFUSION is not given an
# independent random effect: it is derived from DRIFT and the fixed/unit
# marginal covariance through Q = -(A P + P A').

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

.options <- function(args = commandArgs(trailingOnly = TRUE)) {
  x <- list(
    project = normalizePath(getwd(), mustWork = TRUE),
    output = path.expand("~/CoTiMA_Benchmark/ctsemjulia-standardized"),
    development_library = "/private/tmp/ctsem-juliaFit-lib",
    julia_bindir = Sys.getenv("JULIA_BINDIR"), cores = 4L,
    mode = "quick", seed = 20260910L, retain_fits = FALSE
  )
  for (a in args) {
    if (!startsWith(a, "--") || !grepl("=", a, fixed = TRUE))
      stop("Arguments must use --name=value syntax.", call. = FALSE)
    z <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1L]]
    n <- gsub("-", "_", z[[1L]], fixed = TRUE)
    if (!n %in% names(x)) stop("Unknown option: --", z[[1L]], call. = FALSE)
    x[[n]] <- paste(z[-1L], collapse = "=")
  }
  x$cores <- as.integer(x$cores)
  x$seed <- as.integer(x$seed)
  x$mode <- match.arg(tolower(x$mode), c("quick", "full"))
  x$retain_fits <- tolower(as.character(x$retain_fits)) %in%
    c("true", "t", "1", "yes")
  x
}

.drift <- function(raw) -(1e-6 + 2 * log1p(exp(-2 * raw)))
.z <- function(n) as.numeric(scale(qnorm((seq_len(n) - .5) / n)))
.scalar <- function(x) suppressWarnings(as.numeric(x)[1L])

.expm2 <- function(A) {
  if (length(A) == 1L) return(matrix(exp(A[[1L]]), 1L, 1L))
  e <- eigen(A)
  Re(e$vectors %*% diag(exp(e$values), 2L) %*% solve(e$vectors))
}

.simulate <- function(scenario, nstudy, nper, nobs, seed) {
  set.seed(seed)
  scenario <- match.arg(scenario, c("univariate_drift", "bivariate_drift"))
  zz <- .z(nstudy)
  rows <- vector("list", nstudy * nper)
  person <- 0L
  if (scenario == "univariate_drift") {
    raw <- .35 + .30 * zz
    A <- lapply(raw, function(x) matrix(.drift(x), 1L))
    P <- replicate(nstudy, matrix(1), simplify = FALSE)
    truth <- c(drift1_raw_sd = .30)
  } else {
    # Three study effects deliberately use non-collinear deterministic
    # quantiles.  The two diagonal drift effects and the initial correlation
    # vary; marginal variances remain one in every study.
    z2 <- zz[c(seq(2L, nstudy, 2L), seq(1L, nstudy, 2L))]
    r1 <- .35 + .20 * zz
    r2 <- .15 + .16 * z2
    rho_raw <- rep(atanh(.25), nstudy)
    A <- Map(function(x, y) matrix(c(.drift(x), .08, .04, .drift(y)), 2L, 2L),
      r1, r2)
    P <- lapply(tanh(rho_raw), function(r) matrix(c(1, r, r, 1), 2L))
    raw <- cbind(drift1 = r1, drift2 = r2)
    truth <- c(drift1_raw_sd = .20, drift2_raw_sd = .16)
  }
  for (s in seq_len(nstudy)) {
    F <- .expm2(A[[s]])
    V <- (P[[s]] - F %*% P[[s]] %*% t(F))
    V <- (V + t(V)) / 2
    L <- t(chol(V))
    for (j in seq_len(nper)) {
      person <- person + 1L
      x <- drop(t(chol(P[[s]])) %*% rnorm(nrow(P[[s]])))
      yy <- matrix(NA_real_, nobs, length(x))
      for (k in seq_len(nobs)) {
        if (k > 1L) x <- drop(F %*% x + L %*% rnorm(length(x)))
        yy[k, ] <- x
      }
      d <- data.frame(id = person, study = s, time = 0:(nobs - 1L), yy)
      names(d)[-(1:3)] <- paste0("Y", seq_len(ncol(yy)))
      rows[[person]] <- d
    }
  }
  list(data = do.call(rbind, rows), A = A, P = P, raw = raw, truth = truth)
}

.model <- function(scenario) {
  d1 <- "-(1e-06 + 2 * log1p_exp(-(2 * PARS[1,1])))"
  if (scenario == "univariate_drift") {
    m <- ctsem::ctModel(
      type = "ct", id = c("id", "study"), manifestNames = "Y1",
      latentNames = "eta1", LAMBDA = matrix(1),
      PARS = matrix("drift1"),
      DRIFT = matrix(d1), DIFFUSION = matrix(paste0("sqrt(-2 * (", d1, "))")),
      T0VAR = matrix(1), T0MEANS = matrix(0), CINT = matrix(0),
      MANIFESTMEANS = matrix(0), MANIFESTVAR = matrix(0)
    )
    varying <- "drift1"
  } else {
    d2 <- "-(1e-06 + 2 * log1p_exp(-(2 * PARS[2,1])))"
    rho <- "0.25"
    l11 <- paste0("sqrt(-2 * (", d1, " + 0.08 * ", rho, "))")
    q21 <- paste0("-(0.04 + (", d2, ") * ", rho, " + ", rho,
      " * (", d1, ") + 0.08)")
    l21 <- paste0("(", q21, ") / (", l11, ")")
    m <- ctsem::ctModel(
      type = "ct", id = c("id", "study"),
      manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"),
      LAMBDA = diag(2),
      PARS = matrix(c("drift1", "drift2"), 2L, 1L),
      DRIFT = matrix(c(d1, .08, .04, d2), 2L, 2L),
      T0VAR = matrix(c(1, .25, 0, sqrt(1 - .25^2)), 2L, 2L),
      DIFFUSION = matrix(c(
        l11, l21,
        0,
        paste0("sqrt(-2 * ((", d2, ") + 0.04 * ", rho,
          ") - (", l21, ")^2)")
      ), 2L, 2L),
      T0MEANS = matrix(0, 2L, 1L), CINT = matrix(0, 2L, 1L),
      MANIFESTMEANS = matrix(0, 2L, 1L), MANIFESTVAR = matrix(0, 2L, 2L)
    )
    varying <- c("drift1", "drift2")
  }
  m$pars$indvarying <- FALSE
  m$pars$transform[m$pars$matrix == "PARS"] <- "param"
  m$pars$indvarying_study <- FALSE
  m$pars$indvarying_study[m$pars$param %in% varying] <- TRUE
  m
}

.run_one <- function(scenario, dims, opts) {
  design <- .simulate(scenario, dims$nstudy, dims$nper, dims$nobs,
    opts$seed + match(scenario, c("univariate_drift", "bivariate_drift")))
  model <- .model(scenario)
  tm <- system.time({
    fit <- ctsem::ctFit(
      datalong = design$data, model = model, backend = "julia",
      intoverpop = "laplace", optimize = TRUE, priors = FALSE,
      cores = opts$cores,
      optimcontrol = list(estonly = FALSE, carefulfit = 10L), verbose = 0
    )
  })
  ck <- ctsem::ctLaplaceCheck(fit, nodes = 5L)
  sm <- summary(fit)
  tab <- sm[["popsd.study"]] %||% sm[["popsd"]]
  estimates <- if (is.null(tab)) numeric() else setNames(tab[, "mean"], rownames(tab))
  target_parameters <- sub("_raw_sd$", "", names(design$truth))
  list(
    scenario = scenario, status = "completed",
    nstudy = dims$nstudy, nsubjects = dims$nstudy * dims$nper,
    nobs = nrow(design$data), converged = isTRUE(fit$estimate$converged),
    gradient_norm = .scalar(fit$estimate$gradient_norm),
    fit_seconds = unname(tm[["elapsed"]]), laplace_gap = .scalar(ck$gap),
    gap_per_subject = .scalar(ck$gap_per_subject), truth = design$truth,
    estimated_raw_study_sd = setNames(
      estimates[target_parameters], names(design$truth)
    ),
    fit = if (opts$retain_fits) fit else NULL, check = ck
  )
}

.main <- function(opts) {
  if (dir.exists(opts$development_library))
    .libPaths(c(opts$development_library, .libPaths()))
  if (nzchar(opts$julia_bindir)) Sys.setenv(JULIA_BINDIR = opts$julia_bindir)
  Sys.setenv(JULIA_NUM_THREADS = as.character(opts$cores))
  dir.create(opts$output, recursive = TRUE, showWarnings = FALSE)
  devtools::load_all(opts$project, quiet = TRUE)
  ctsem::ctJuliaSetup(threads = opts$cores, force = TRUE)
  dims <- if (opts$mode == "quick") {
    list(nstudy = 8L, nper = 20L, nobs = 12L)
  } else list(nstudy = 16L, nper = 35L, nobs = 20L)
  scenarios <- c("univariate_drift", "bivariate_drift")
  out <- lapply(scenarios, function(s) tryCatch(
    .run_one(s, dims, opts),
    error = function(e) list(scenario = s, status = "failed",
      error = conditionMessage(e))
  ))
  names(out) <- scenarios
  saveRDS(list(options = opts, dimensions = dims, results = out),
    file.path(opts$output, "standardized-validation-results.rds"))
  lines <- c(
    "# CoTiMA2 covariance-stationary standardized validation", "",
    paste0("Generated: ", format(Sys.time(), tz = "UTC", usetz = TRUE)), "",
    paste0("- Mode: ", opts$mode), paste0("- Julia cores: ", opts$cores), "",
    paste("This validation fixes latent marginal variances at one. DIFFUSION is",
      "derived from DRIFT and P by Q = -(A P + P A'); it has no independent",
      "study random effect. The bivariate scenario permits study variation in",
      "two diagonal DRIFT parameters and fixes the T0 correlation at 0.25.",
      "The present native Julia API rejects parameter-dependent expressions in",
      "T0VAR, so a varying T0 correlation with both marginal variances fixed",
      "exactly at one remains an explicitly documented capability gap."), "",
    "| Scenario | Status | Converged | Fit s | Gap/person |", "|---|---|---:|---:|---:|"
  )
  for (x in out) lines <- c(lines, sprintf("| %s | %s | %s | %.3f | %.6g |",
    x$scenario, x$status, x$converged %||% NA, x$fit_seconds %||% NA,
    x$gap_per_subject %||% NA))
  writeLines(lines, file.path(opts$output, "REPORT.md"))
  print(out)
  message("Validation report: ", file.path(opts$output, "REPORT.md"))
  invisible(out)
}

if (sys.nframe() == 0L) .main(.options())
