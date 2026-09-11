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

.drift <- function(raw) -(.2 + 2 * log1p(exp(-2 * raw)))
.bounded_drift <- function(raw) -(.2 + .8 / (1 + exp(-raw)))
.z <- function(n) as.numeric(scale(qnorm((seq_len(n) - .5) / n)))
.scalar <- function(x) suppressWarnings(as.numeric(x)[1L])

.match_study_sd_parameters <- function(truth_names, estimate_names) {
  vapply(truth_names, function(x) {
    candidates <- unique(c(
      sub("_sd$", "", x),
      sub("_raw_sd$", "", x)
    ))
    hit <- candidates[candidates %in% estimate_names]
    if (length(hit)) hit[[1L]] else NA_character_
  }, character(1L), USE.NAMES = FALSE)
}

.expm2 <- function(A) {
  if (length(A) == 1L) return(matrix(exp(A[[1L]]), 1L, 1L))
  e <- eigen(A)
  Re(e$vectors %*% diag(exp(e$values), 2L) %*% solve(e$vectors))
}

.simulate <- function(scenario, nstudy, nper, nobs, seed) {
  set.seed(seed)
  scenario <- match.arg(scenario, c(
    "univariate_drift", "bivariate_drift", "bivariate_drift_correlation"
  ))
  zz <- .z(nstudy)
  rows <- vector("list", nstudy * nper)
  person <- 0L
  if (scenario == "univariate_drift") {
    raw <- .35 + .30 * zz
    A <- lapply(raw, function(x) matrix(.drift(x), 1L))
    P <- replicate(nstudy, matrix(1), simplify = FALSE)
    truth <- c(drift1_raw_sd = .30)
  } else if (scenario == "bivariate_drift") {
    # Three study effects deliberately use non-collinear deterministic
    # quantiles.  The two diagonal drift effects and the initial correlation
    # vary; marginal variances remain one in every study.
    # A cosine sequence is exactly uncorrelated with the symmetric ordered
    # normal quantiles, avoiding a nearly singular two-effect design.
    z2 <- as.numeric(scale(cos(2 * pi * (seq_len(nstudy) - .5) / nstudy)))
    r1 <- .35 + .20 * zz
    r2 <- .15 + .16 * z2
    rho_raw <- rep(0, nstudy)
    A <- Map(function(x, y) matrix(
      c(.drift(x), .08, .04, .drift(y)), 2L, 2L, byrow = TRUE
    ),
      r1, r2)
    P <- lapply(tanh(rho_raw), function(r) matrix(c(1, r, r, 1), 2L))
    raw <- cbind(drift1 = r1, drift2 = r2)
    truth <- c(drift1_raw_sd = .20, drift2_raw_sd = .16)
  } else {
    # Three mutually orthogonal deterministic study scores isolate recovery of
    # both drift coordinates and the raw stationary-correlation coordinate.
    z2 <- as.numeric(scale(cos(2 * pi * (seq_len(nstudy) - .5) / nstudy)))
    candidate <- sin(4 * pi * (seq_len(nstudy) - .5) / nstudy)
    z3 <- as.numeric(scale(qr.resid(qr(cbind(1, zz, z2)), candidate)))
    r1 <- .25 + .20 * zz
    r2 <- -.10 + .16 * z2
    rho_raw <- .18 * z3
    rho <- .3 * tanh(rho_raw)
    A <- Map(function(x, y) matrix(
      c(.bounded_drift(x), .08, .04, .bounded_drift(y)),
      2L, 2L, byrow = TRUE
    ), r1, r2)
    P <- lapply(rho, function(r) matrix(c(1, r, r, 1), 2L))
    raw <- cbind(drift1 = r1, drift2 = r2, rho_raw = rho_raw)
    truth <- c(
      drift1_raw_sd = .20, drift2_raw_sd = .16, rho_raw_sd = .18
    )
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
  # With positive fixed cross-effects, merely constraining each diagonal below
  # zero is insufficient to keep Q positive definite throughout optimization.
  # A -0.2 upper bound is safely inside the admissible region for this design.
  d1 <- "-(0.2 + 2 * log1p_exp(-(2 * PARS[1,1])))"
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
  } else if (scenario == "bivariate_drift") {
    d2 <- "-(0.2 + 2 * log1p_exp(-(2 * PARS[2,1])))"
    l11 <- paste0("sqrt(-2 * (", d1, "))")
    q21 <- "-(0.04 + 0.08)"
    l21 <- paste0("(", q21, ") / (", l11, ")")
    m <- ctsem::ctModel(
      type = "ct", id = c("id", "study"),
      manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"),
      LAMBDA = diag(2),
      PARS = matrix(c("drift1", "drift2"), 2L, 1L),
      DRIFT = matrix(c(d1, .08, .04, d2), 2L, 2L, byrow = TRUE),
      T0VAR = diag(2),
      DIFFUSION = matrix(c(
        l11, l21,
        0,
        paste0("sqrt(-2 * (", d2, ") - (", l21, ")^2)")
      ), 2L, 2L),
      T0MEANS = matrix(0, 2L, 1L), CINT = matrix(0, 2L, 1L),
      MANIFESTMEANS = matrix(0, 2L, 1L), MANIFESTVAR = matrix(0, 2L, 2L)
    )
    varying <- c("drift1", "drift2")
  } else {
    # Bounded diagonal drifts and |rho| < .3 make the derived Q positive
    # definite over the complete admissible parameter domain for the fixed
    # cross-effects used here, rather than merely at the starting values.
    d1 <- "-(0.2 + 0.8 / (1 + exp(-PARS[1,1])))"
    d2 <- "-(0.2 + 0.8 / (1 + exp(-PARS[2,1])))"
    rho <- "(0.3 * tanh(PARS[3,1]))"
    q11 <- paste0("-2 * (", d1, " + 0.08 * ", rho, ")")
    q22 <- paste0("-2 * (", d2, " + 0.04 * ", rho, ")")
    q21 <- paste0("-(0.08 + 0.04 + ", rho, " * (", d1, " + ", d2, "))")
    l11 <- paste0("sqrt(", q11, ")")
    l21 <- paste0("(", q21, ") / (", l11, ")")
    l22 <- paste0("sqrt(", q22, " - (", l21, ")^2)")
    m <- ctsem::ctModel(
      type = "ct", id = c("id", "study"),
      manifestNames = c("Y1", "Y2"), latentNames = c("eta1", "eta2"),
      LAMBDA = diag(2),
      PARS = matrix(c("drift1", "drift2", "rho_raw"), 3L, 1L),
      DRIFT = matrix(c(d1, .08, .04, d2), 2L, 2L, byrow = TRUE),
      T0VAR = matrix(c(1, 0, rho, 1), 2L, 2L, byrow = TRUE),
      DIFFUSION = matrix(c(l11, l21, 0, l22), 2L, 2L),
      T0MEANS = matrix(0, 2L, 1L), CINT = matrix(0, 2L, 1L),
      MANIFESTMEANS = matrix(0, 2L, 1L), MANIFESTVAR = matrix(0, 2L, 2L)
    )
    varying <- c("drift1", "drift2", "rho_raw")
  }
  m$pars$indvarying <- FALSE
  m$pars$transform[m$pars$matrix == "PARS"] <- "param"
  m$pars$indvarying_study <- FALSE
  m$pars$indvarying_study[m$pars$param %in% varying] <- TRUE
  m
}

.run_one <- function(scenario, dims, opts) {
  design <- .simulate(scenario, dims$nstudy, dims$nper, dims$nobs,
    opts$seed + match(scenario, c(
      "univariate_drift", "bivariate_drift", "bivariate_drift_correlation"
    )))
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
  target_parameters <- .match_study_sd_parameters(
    names(design$truth), names(estimates)
  )
  estimated_sd <- rep(NA_real_, length(target_parameters))
  matched <- !is.na(target_parameters)
  estimated_sd[matched] <- estimates[target_parameters[matched]]
  identifiability <- fit$identifiability %||% NULL
  list(
    scenario = scenario, status = "completed",
    nstudy = dims$nstudy, nsubjects = dims$nstudy * dims$nper,
    nobs = nrow(design$data), converged = isTRUE(fit$estimate$converged),
    gradient_norm = .scalar(fit$estimate$gradient_norm),
    fit_seconds = unname(tm[["elapsed"]]), laplace_gap = .scalar(ck$gap),
    gap_per_subject = .scalar(ck$gap_per_subject), truth = design$truth,
    estimated_raw_study_sd = setNames(estimated_sd, names(design$truth)),
    identifiability = identifiability,
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
  scenarios <- c(
    "univariate_drift", "bivariate_drift", "bivariate_drift_correlation"
  )
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
      "study random effect. The original bivariate scenario permits study",
      "variation in two diagonal DRIFT parameters and fixes the T0 correlation",
      "at zero. The correlation scenario additionally varies a bounded T0",
      "correlation through a static PARS expression while retaining both",
      "marginal variances at exactly one; DIFFUSION is derived from the same",
      "three study-varying raw coordinates."), "",
    paste("An `NA` in the recovered-SD column means that the requested SD could",
      "not be matched to the fitted summary; it does not by itself diagnose",
      "non-identifiability. Identifiability warnings reported by ctsem are",
      "listed separately below."), "",
    "| Scenario | Status | Converged | Gradient norm | Fit s | Gap/person | Recovered study SDs |",
    "|---|---|---:|---:|---:|---:|---|"
  )
  for (x in out) {
    recovered <- if (is.null(x$estimated_raw_study_sd)) "" else paste(
      sprintf("%s=%.3f", names(x$estimated_raw_study_sd),
        x$estimated_raw_study_sd), collapse = "; "
    )
    lines <- c(lines, sprintf("| %s | %s | %s | %.6g | %.3f | %.6g | %s |",
      x$scenario, x$status, x$converged %||% NA,
      x$gradient_norm %||% NA, x$fit_seconds %||% NA,
      x$gap_per_subject %||% NA, recovered))
  }
  identification_notes <- unlist(lapply(out, function(x) {
    id <- x$identifiability %||% NULL
    bad <- id$parameters %||% id$not_estimable %||%
      id$notEstimable %||% character()
    if (!length(bad)) return(NULL)
    paste0("- `", x$scenario, "`: ", paste(bad, collapse = ", "))
  }))
  lines <- c(lines, "", "## Identifiability warnings", "",
    if (length(identification_notes)) identification_notes else
      "No explicit non-estimable parameter directions were returned.")
  writeLines(lines, file.path(opts$output, "REPORT.md"))
  print(out)
  message("Validation report: ", file.path(opts$output, "REPORT.md"))
  invisible(out)
}

if (sys.nframe() == 0L) .main(.options())
