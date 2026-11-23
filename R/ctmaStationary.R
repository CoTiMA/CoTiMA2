# The default path is deliberately a no-op. A stationary covariance is the
# covariance of the latent dynamic process, not the covariance of observed
# standardized data or the covariance of random intercepts.
.ctmaLogicalFlag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("'", name, "' must be one non-missing TRUE or FALSE.", call. = FALSE)
  }
  x
}

.ctmaStationaryOptions <- function(stationary, determinedDiffusions,
    stationary_missing = FALSE, backend_stationary = NULL,
    backend = "stan", optimize = TRUE, priors = FALSE, inits = NULL,
    stationaryVariance = NULL) {
  stationary <- .ctmaLogicalFlag(stationary, "stationary")
  determinedDiffusions <- .ctmaLogicalFlag(determinedDiffusions, "determinedDiffusions")
  if (!is.null(backend_stationary)) {
    backend_stationary <- .ctmaLogicalFlag(backend_stationary, "backendArgs$stationary")
    if (stationary_missing) stationary <- backend_stationary else if (
        !identical(stationary, backend_stationary)) {
      stop("Conflicting top-level and backend stationary settings.", call. = FALSE)
    }
  }
  if ((stationary || determinedDiffusions) && identical(backend, "julia")) {
    stop("These stationary options are not supported by the bundled backend='julia'. ",
      "Use backend='ctsemjulia' with compatible experimental ctsem, or backend='stan' ",
      "with a ctsem version that supports stationary initialization.", call. = FALSE)
  }
  if (!is.null(stationaryVariance)) {
    if (!is.numeric(stationaryVariance) || length(stationaryVariance) != 1L ||
        !is.finite(stationaryVariance) || stationaryVariance <= 0)
      stop("stationaryVariance must be NULL or one finite positive variance.", call. = FALSE)
    if (!stationary || !determinedDiffusions || backend != "ctsemjulia")
      stop("stationaryVariance requires stationary=TRUE, determinedDiffusions=TRUE and backend='ctsemjulia'.", call. = FALSE)
  }
  if (determinedDiffusions) {
    if (!stationary) stop("determinedDiffusions=TRUE currently requires stationary=TRUE. ",
      "The estimated stationary within-person covariance is not a fixed unit matrix.",
      call. = FALSE)
    if ((!isTRUE(optimize) || isTRUE(priors)) && backend != "ctsemjulia") stop(
      "determinedDiffusions Bayesian/MAP estimation requires backend='ctsemjulia'; other backends support maximum likelihood only.", call. = FALSE)
    if (!isTRUE(optimize) && !isTRUE(priors)) stop(
      "Bayesian determinedDiffusions requires priors=TRUE on the stationary raw coordinates.", call. = FALSE)
    if (isTRUE(priors)) message(
      "Determined-diffusion priors are native ctsem priors on the stationary raw coordinates (including log variance/shape and study contrasts). They are not transformed priors from an unrestricted diffusion model; no extra change-of-variables Jacobian is added.")
    if (!is.null(inits)) stop(
      "determinedDiffusions changes the raw parameter coordinates. ",
      "Old raw starting values cannot be reused; omit inits.", call. = FALSE)
  }
  list(stationary = stationary, determinedDiffusions = determinedDiffusions)
}

.ctmaDeterminedDiffusionModel <- function(model, determinedDiffusions = FALSE,
    study_effect_columns = NULL, study_design = NULL, stationaryVariance = NULL) {
  if (!isTRUE(determinedDiffusions)) return(model)
  if (!is.null(stationaryVariance)) return(.ctmaFixedStationaryVarianceModel(model, stationaryVariance))
  p <- model$pars
  n <- model$n.latent
  fail <- function(...) stop("determinedDiffusions: ", ..., call. = FALSE)
  if (!is.data.frame(p) || length(n) != 1L || !n %in% 1:2) {
    fail("only one or two dynamic latent processes are currently supported; ",
      "phantom random-intercept factors and higher dimensions are not supported.")
  }
  if (length(model$manifesttype) && any(model$manifesttype != 0)) {
    fail("this first implementation supports Gaussian manifests only; binary/ordinal measurements are not yet validated.")
  }
  # No untested expression graphs or existing auxiliary parameters are rewritten.
  if (any(p$matrix == "PARS" & (!is.na(p$param) | p$value != 0), na.rm = TRUE)) {
    fail("existing PARS/custom expression models require an explicit model-specific mapping.")
  }
  p <- p[p$matrix != "PARS", , drop = FALSE]
  dynamic <- p$matrix %in% c("DRIFT", "DIFFUSION")
  varying <- grep("^indvarying($|_)", names(p), value = TRUE)
  if (any(vapply(varying, function(v) any(p[[v]][dynamic] %in% TRUE), logical(1)))) {
    fail("individual/study random effects on DRIFT or DIFFUSION are not yet supported.")
  }
  if (any(vapply(varying, function(v) {
    any(p[[v]][p$matrix == "T0MEANS" & !is.na(p$param)] %in% TRUE)
  }, logical(1)))) {
    fail("random initial T0MEANS are not covered by this stationary parameterization. ",
      "For random manifest intercepts, use indVaryingT0=FALSE.")
  }
  if (any(vapply(varying, function(v) {
    any(p[[v]][p$matrix=="CINT" & !is.na(p$param)] %in% TRUE)
  },logical(1)))) {
    fail("random CINT effects are not yet validated; supported random intercepts use MANIFESTMEANS.")
  }
  ti <- grep("_effect$", names(p), value = TRUE)
  study_ti <- .ctmaDeterminedStudyEffects(p, ti, study_effect_columns,
    study_design, n, fail)
  free_label <- function(rows) {
    z <- p$param[rows]
    length(z) && !anyNA(z) && !anyDuplicated(z) &&
      all(grepl("^[[:alpha:]_][[:alnum:]_]*$", z))
  }
  dr <- which(p$matrix == "DRIFT")
  di <- which(p$matrix == "DIFFUSION" & p$row >= p$col)
  if (length(dr) != n*n || !free_label(dr) || !free_label(di)) {
    fail("all DRIFT cells and lower-triangular DIFFUSION cells must have distinct ",
      "free labels. Fixed, equal or expression-constrained process parameters are not rewritten.")
  }
  # Only the standard ctsem transforms have the tested scientific interpretation.
  allowed <- function(rows, kind) {
    for (i in rows) {
      s <- gsub("[[:space:]]", "", p$transform[i])
      s <- gsub("0.000001", "1e-06", s, fixed=TRUE)
      s <- gsub("0.0000000001", "1e-10", s, fixed=TRUE)
      ok <- if (kind == "drift") {
        if (p$row[i] == p$col[i]) s %in% c("-(1e-06+2*log1p_exp(-(2*param)))",
          "-(1e-6+2*log1p_exp(-(2*param)))") else s %in% c("param","2/(1+exp(-param))-1")
      } else {
        if (p$row[i] == p$col[i]) s %in% c("10*log1p_exp(2*param)+1e-10","1e-10+10*log1p_exp(2*param)") else s %in% c("param","2/(1+exp(-param))-1")
      }
      if (is.na(ok) || !ok) fail("custom process transform for ", p$matrix[i], "[",p$row[i],",",p$col[i],"] (",p$transform[i],") requires an explicit mapping.")
    }
  }
  allowed(dr, "drift"); allowed(di, "diffusion")
  # Switching covariance coordinates must not reinterpret fixed nonzero
  # off-diagonal measurement or population covariance constraints.
  if (any(p$matrix == "MANIFESTVAR" & p$row != p$col &
      (is.na(p$value) | p$value != 0))) {
    fail("nonzero/free measurement-error covariances are not supported in this first implementation.")
  }
  has_random <- any(vapply(varying, function(v) any(p[[v]] %in% TRUE), logical(1)))
  preserve_z <- has_random && identical(model$covmattransform, "z")
  if (has_random && !preserve_z && !identical(model$covmattransform, "cholesky")) {
    fail("random intercepts require the native z covariance construction (or an already explicit Cholesky model); covariance safeguards must not be silently changed.")
  }
  rp <- model$RAWPOPVAR
  if (!preserve_z && !is.null(rp) && length(rp) && is.matrix(rp)) {
    off <- rp[lower.tri(rp)]
    numeric_off <- suppressWarnings(as.numeric(off))
    if (any(!is.na(numeric_off) & numeric_off != 0)) {
      fail("fixed nonzero random-effect covariances need a model-specific covariance-coordinate conversion.")
    }
  }

  # Keep the original drift raw coordinates and their exact transforms.
  # For stable A and PD shape H (H11=1), solve AX+XA'=-H.
  # Sigma=v*X/trace(X), Q=v*H/trace(X) parameterizes precisely the admissible
  # stationary covariance cone. v is its total within-person variance.
  drift_names <- as.character(p$param[dr])
  if (any(drift_names %in% p$param[-dr]) ||
      any(p$param[di] %in% p$param[-which(p$matrix=="DIFFUSION")])) {
    fail("cross-matrix equality constraints need an explicit mapping.")
  }
  drift_transforms <- as.character(p$transform[dr])
  cell <- function(r,c) {
    at <- which(p$matrix=="DRIFT" & p$row==r & p$col==c)
    paste0("(",p$param[at],")")
  }
  A <- matrix("0",n,n)
  for (i in dr) A[p$row[i],p$col[i]] <- paste0("1*",p$param[i])
  if (n == 1L) {
    extra <- "ctma_eq_logvariance"
    aux <- c(drift_names,extra)
    aux_transforms <- c(drift_transforms,"param")
    rate <- cell(1,1)
    L <- matrix(paste0("sqrt(-2*",rate,"*exp(",extra,"))"))
    S <- matrix(paste0("sqrt(exp(",extra,"))"))
    starts <- c(rep(0,length(dr)),log(.5))
  } else {
    extra <- c("ctma_eq_logtotalvariance","ctma_eq_shape21","ctma_eq_logshape22")
    aux <- c(drift_names,extra)
    aux_transforms <- c(drift_transforms,rep("param",3))
    a<-cell(1,1);b<-cell(1,2);c<-cell(2,1);d<-cell(2,2)
    h12<-extra[2]
    h22<-paste0("(",h12,"*",h12,"+exp(2*",extra[3],"))")
    den<-paste0("((",a,"+",d,")*(",a,"*",d,"-",b,"*",c,"))")
    y<-paste0("(",c,"*",d,"/2+",a,"*",b,"*",h22,"/2-",
      a,"*",d,"*",h12,")/",den)
    x<-paste0("(-0.5-",b,"*(",y,"))/",a)
    z<-paste0("(-",h22,"/2-",c,"*(",y,"))/",d)
    tr<-paste0("((",x,")+(",z,"))")
    # For valid stationary A, trace(X)>0. abs only protects Q evaluation
    # during invalid proposals; stationary initial covariance still rejects
    # unstable A rather than turning it into a different admissible model.
    scale<-paste0("(exp(",extra[1],")/abs(",tr,"))")
    s11<-paste0("(",scale,"*(",x,"))")
    s12<-paste0("(",scale,"*(",y,"))")
    s22<-paste0("(",scale,"*(",z,"))")
    # Q=scale*H is algebraically -A Sigma-Sigma A', without cancellation.
    l11<-paste0("sqrt(",scale,")")
    L<-matrix(c(l11,paste0(h12,"*",l11),"0",
      paste0("exp(",extra[3],")*",l11)),2)
    t11<-paste0("sqrt(",s11,")")
    t21<-paste0("(",s12,")/",t11)
    S<-matrix(c(t11,t21,"0",paste0("sqrt(",s22,"-(",t21,")*(",t21,"))")),2)
    starts<-c(rep(0,length(dr)),log(1),0,0)
  }
  if (any(extra %in% p$param)) fail("reserved auxiliary parameter names already exist.")
  original <- p[p$matrix %in% c("DRIFT", "DIFFUSION", "T0VAR"), , drop = FALSE]
  write_cells <- function(mat, cells) {
    ix <- which(p$matrix == mat)
    for (i in ix) {
      x <- cells[p$row[i], p$col[i]]
      isnum <- !is.na(suppressWarnings(as.numeric(x)))
      p$param[i] <<- if (isnum) NA_character_ else x
      p$value[i] <<- if (isnum) as.numeric(x) else NA_real_
      p$transform[i] <<- NA_character_
      for (v in varying) p[[v]][i] <<- FALSE
      for (v in ti) p[[v]][i] <<- FALSE
    }
  }
  # Preserve the model-wide z construction for random effects. Only process
  # covariance cells are re-expressed: z uses SDs and Fisher-z in dimension2.
  # This keeps native RI correlation-bound penalties and diagnostics active.
  if (preserve_z && n == 2L) {
    # Work directly with covariance entries. Re-expanding Cholesky factors
    # into covariances duplicates large expressions (and can exceed R's
    # 10000-byte environment-name limit in ctsem's matrix setup).
    fisher <- function(rho) paste0("(log(1+",rho,")-log(1-",rho,"))/2")
    rho_q <- paste0("(",h12,"/sqrt(",h22,"))")
    rho_s <- paste0("((",y,")/sqrt((",x,")*(",z,")))")
    # Positive scale cancels in both correlations; marginal SDs still use it.
    L <- matrix(c(paste0("sqrt(",scale,")"),fisher(rho_q),"0",
      paste0("sqrt(",scale,"*",h22,")")),2)
    S <- matrix(c(paste0("sqrt(",s11,")"),fisher(rho_s),"0",
      paste0("sqrt(",s22,")")),2)
  }
  write_cells("DRIFT", A); write_cells("DIFFUSION", L)
  # These are explicit Sigma cells as well as the stationary distribution's
  # parameterization; ctsem stationary=TRUE derives the same Sigma from A,Q.
  write_cells("T0VAR", S)
  for (i in seq_along(aux)) {
    row <- p[di[1], , drop = FALSE]
    row$matrix <- "PARS"; row$row <- i; row$col <- 1L
    row$param <- aux[i]; row$value <- NA_real_
    row$transform <- aux_transforms[i]; row$indvarying <- FALSE; row$sdscale <- 1
    for (v in varying) row[[v]] <- FALSE
    for (v in ti) row[[v]] <- FALSE
    if (aux[i] %in% extra) for (v in study_ti) row[[v]] <- TRUE
    p <- rbind(p, row)
  }
  model$pars <- p
  model$covmattransform <- if (preserve_z) "z" else "cholesky"
  model$determinedDiffusions <- list(
    parameterization = "drift-preserving admissible stationary within-person covariance",
    dimension = n, original = original, drift_parameters = drift_names,
    prior_coordinates = "native ctsem priors on drift raw coordinates, log stationary variance/shape and study contrasts; not prior-equivalent to unrestricted DIFFUSION",
    study_effect_columns = study_ti,
    study_covariances = if (length(study_ti)) "fully study-specific" else "common",
    drift_constraints_preserved = TRUE,
    random_covariance_construction_preserved = has_random,
    auxiliary = aux, suggested_raw = stats::setNames(starts, aux),
    note = paste("Diffusion covariance is derived from drift and an estimated stationary",
      "within-person covariance, not a unit observed covariance. Original drift",
      "parameters and transforms, including negative diagonal constraints, are preserved.",
      "The covariance coordinates do not reduce the number of free parameters."))
  message("Note: Determined diffusion preserves the original drift parameters and ",
    "constraints. It estimates a stationary within-person covariance; ",
    "random-intercept covariance remains separate.")
  if (length(study_ti)) .ctma_print_message("Study contrasts now act on the stationary within-person covariance coordinates. Drift is common unless study-level random DRIFT effects are explicitly selected; derived diffusion is study-specific. Raw TI coefficients use the new coordinates.")
  model
}

.ctmaStationaryResult <- function(result, stationary, determinedDiffusions, stationaryVariance = NULL) {
  if (!stationary && !determinedDiffusions) return(result)
  note <- if (determinedDiffusions) paste(
    "Diffusion is derived from drift and an estimated stationary within-person",
    "covariance; observed standardization does not fix that covariance to one.",
    "Original drift parameters and transforms are preserved. Random-intercept",
    "covariance is separate; identification is not guaranteed.") else paste(
    "Stationary initial distribution requested: ctsem derives initial moments",
    "from drift and diffusion; this does not fix within-person variance to one.")
  if (!is.null(stationaryVariance)) note <- paste("Stationary univariate process variance fixed to", stationaryVariance,
    "; diffusion variance rate is -2 * DRIFT * stationaryVariance. Original DRIFT transform is preserved.")
  spec <- list(stationary = stationary, determinedDiffusions = determinedDiffusions,
    stationaryVariance = stationaryVariance, note = note)
  chart <- result$ctModel$determinedDiffusions
  if (is.list(chart) && !is.null(chart$study_covariances)) {
    spec$study_covariances <- chart$study_covariances
    spec$study_effect_columns <- chart$study_effect_columns
  }
  if (determinedDiffusions) spec$prior_coordinates <- paste(
    "When priors are enabled: native ctsem priors on the stationary raw coordinates,",
    "including drift, log variance/shape and study contrasts. Not equivalent to",
    "priors on unrestricted diffusion coordinates; no additional Jacobian.")
  result$stationary_specification <- spec
  for (nm in c("argumentList", "argumentListOfCtmaInit")) {
    if (is.list(result[[nm]])) {
      result[[nm]]$stationary <- stationary
      result[[nm]]$determinedDiffusions <- determinedDiffusions
      result[[nm]]$stationaryVariance <- stationaryVariance
    }
  }
  if (is.list(result$summary)) result$summary$stationary_specification <- spec
  result
}

# Published ctsem 3.11.1 explicitly disables its stationary flag. For our
# determined model, Sigma is already an explicit T0VAR expression. With zero
# CINT the stationary mean is zero too, so using these initial cells with the
# native flag off is exact, not a covariance approximation or fixed-Sigma fit.
.ctmaStationaryFitArgs <- function(args) {
  if (!isTRUE(args$stationary)) return(args)
  native_formals <- names(formals(ctsem::ctFit))
  backend <- if (is.null(args$backend)) "stan" else args$backend
  old_block <- grepl("Stationary option temporarily unavailable",
    paste(deparse(body(ctsem::ctFit)),collapse="\n"), fixed=TRUE)
  derived <- is.list(args$model$determinedDiffusions)
  if (derived && (identical(backend,"stan") || old_block)) {
    m <- args$model
    p <- m$pars
    ci <- p[p$matrix=="CINT",,drop=FALSE]
    if (nrow(ci) != m$n.latent || anyNA(ci$value) || any(ci$value != 0)) {
      stop("The Stan compatibility path for determinedDiffusions currently requires ",
        "CINT fixed to zero. Use compatible native ctsem Julia for nonzero/free CINT.",
        call.=FALSE)
    }
    ix <- which(p$matrix=="T0MEANS")
    p$param[ix] <- NA_character_; p$value[ix] <- 0
    p$transform[ix] <- NA_character_
    for (v in grep("^indvarying($|_)|_effect$",names(p),value=TRUE))p[[v]][ix]<-FALSE
    m$pars <- p
    m$determinedDiffusions$stationary_implementation <- "explicit initial covariance and zero mean"
    args$model <- m
    args$stationary <- FALSE
    if ("backend" %in% native_formals) args$backend <- "stan"
    message("Note: Stationary initialization is specified explicitly through the ",
      "estimated covariance and zero mean; ctsem's unavailable/native-only ",
      "stationary switch is not needed on this Stan path.")
  } else if (old_block) {
    stop("This installed ctsem version disables stationary=TRUE. ",
      "Use compatible experimental ctsem with CoTiMA2 backend='ctsemjulia', ",
      "or the supported zero-CINT determinedDiffusions=TRUE parameterization. ",
      "The stationary request was not silently ignored.",call.=FALSE)
  } else if (identical(backend,"stan") && "backend" %in% native_formals) {
    args$backend <- "stan"
  }
  args
}

# An explicit FALSE disables only individual initial-mean random effects.
# NULL/TRUE retain the existing automatic choice; do not alter T0VAR or RI means.
.ctmaApplyT0MeanVariation <- function(model, indVaryingT0=NULL) {
  if (is.null(indVaryingT0)) return(model)
  indVaryingT0 <- .ctmaLogicalFlag(indVaryingT0, "indVaryingT0")
  if (!indVaryingT0) {
    rows <- model$pars$matrix %in% "T0MEANS"
    model$pars$indvarying[rows] <- FALSE
  }
  model
}


.ctmaDeterminedStudyEffects <- function(p, ti, study_columns, design, n, fail) {
  dynamic <- p$matrix %in% c("DRIFT", "DIFFUSION")
  active <- ti[vapply(ti, function(v) any(p[[v]][dynamic] %in% TRUE), logical(1))]
  if (!length(active)) return(character())
  if (any(vapply(ti, function(v) any(p[[v]][p$matrix == "DRIFT"] %in% TRUE), logical(1))))
    fail("study/predictor effects on DRIFT cannot be moved without changing invariance constraints; the study-specific covariance path requires common drift.")
  if (is.null(study_columns) || !length(study_columns) ||
      anyDuplicated(study_columns) || !all(study_columns %in% ti) ||
      !setequal(active, study_columns))
    fail("study/predictor effects on DRIFT or DIFFUSION require an explicit saturated study design to preserve invariance constraints. General moderators and partially varying covariance structures are not rewritten.")
  di <- which(p$matrix == "DIFFUSION" & p$row >= p$col)
  if (length(di) != n * (n + 1L) / 2L ||
      !all(vapply(study_columns, function(v) all(p[[v]][di] %in% TRUE), logical(1))))
    fail("all free diffusion covariance coordinates must vary across all study contrasts. Partial covariance invariance constraints require an explicit mapping.")
  if (!is.matrix(design) || ncol(design) != length(study_columns) ||
      !is.numeric(design) || any(!is.finite(design)))
    fail("a finite numeric study-contrast design is required.")
  patterns <- unique(cbind(Intercept = 1, design))
  if (nrow(patterns) != ncol(patterns) || qr(patterns)$rank != ncol(patterns))
    fail("study contrasts must span one unrestricted covariance per study; the design is not saturated and full rank.")
  # Full-rank saturated dummy/weighted contrasts specify arbitrary coordinate
  # values for each study. Thus a bijective covariance chart preserves the
  # model set and the number of parameters (not raw TI coefficient meanings).
  study_columns
}
