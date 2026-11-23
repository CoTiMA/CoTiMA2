# Internal research adapter for the Lohmann, Luedtke & Robitzsch supplement.
# No @export: intentionally accessed through CoTiMA2:::ctmaLohmannFit().
# Original Stan strings are read as literals, never by sourcing the compiler.
ctmaLohmannFit <- function(primaryStudies = NULL, reStructure = "URE", n.latent = 2L,
    supplementDirectory = NULL, correlationMatrices = NULL, sampleSizes = NULL,
    timeLag = NULL, samplingCovariances = NULL, stanData = NULL,
    fit = TRUE, chains = 4L, iter = 2000L, warmup = floor(iter / 2),
    coresToUse = 1L, seed = 20261106L, control = list(), init = 0,
    stanModel = NULL, activeDirectory = NULL, ...) {
  if (length(list(...))) stop("Unknown arguments; ctmaFit model settings are not accepted.")
  reStructure <- toupper(reStructure)
  if (length(reStructure) != 1L || is.na(reStructure) ||
      !reStructure %in% c("FE", "URE", "ERE", "CRE")) stop("reStructure must be FE, URE, ERE or CRE.")
  if (!is.logical(fit) || length(fit) != 1L || is.na(fit)) stop("fit must be TRUE or FALSE.")
  if (!is.null(activeDirectory) && (!is.character(activeDirectory) ||
      length(activeDirectory) != 1L || is.na(activeDirectory) || !dir.exists(activeDirectory)))
    stop("activeDirectory must be an existing directory.")
  if (!is.numeric(n.latent) || length(n.latent) != 1L || is.na(n.latent) ||
      !n.latent %in% c(1,2)) stop("n.latent must be 1 or 2.")
  source <- if (n.latent == 2L) .ctmaLohmannSource(supplementDirectory, reStructure) else
    list(code = .ctmaLohmannUnivariateCode(reStructure), path = NA_character_,
      md5 = NA_character_)
  if (!is.null(stanData)) {
    if (any(!vapply(list(primaryStudies, correlationMatrices, sampleSizes,
        timeLag, samplingCovariances), is.null, logical(1))))
      stop("Supply stanData alone, or supply study inputs, not both.")
    dat <- stanData
  } else {
    if (!is.null(primaryStudies)) {
      if (!is.null(correlationMatrices) || !is.null(sampleSizes) || !is.null(timeLag))
        stop("Do not combine primaryStudies with matrices, sampleSizes or timeLag.")
      p <- primaryStudies
      if (is.null(p$empcovs)) p <- ctmaFitToPrep(p)
      correlationMatrices <- p$empcovs
      sampleSizes <- unlist(p$sampleSizes, use.names = FALSE)
      timeLag <- vapply(p$deltas, function(x) {
        if (!is.numeric(x) || !length(x) || any(!is.finite(x)) ||
            any(x <= 0) || max(abs(x - x[1])) > 1e-8 * max(1, abs(x[1])))
          stop("The original model requires equally spaced waves within each study.")
        x[1]
      }, numeric(1))
      if (length(p$deltas) != length(correlationMatrices) ||
          any(vapply(seq_along(p$deltas), function(i)
            length(p$deltas[[i]]) + 1L != nrow(correlationMatrices[[i]]) / n.latent, logical(1))))
        stop("Study intervals and correlation-matrix dimensions do not match.")
    }
    dat <- .ctmaLohmannData(correlationMatrices, sampleSizes, timeLag, samplingCovariances, n.latent)
  }
  # rstan accepts lists of matrices as arrays; normalize the authors' format.
  if (is.list(dat$V_i) && !is.data.frame(dat$V_i)) {
    vv <- dat$V_i
    if (!length(vv) || !all(vapply(vv, is.matrix, logical(1))) ||
        length(unique(vapply(vv, function(x) paste(dim(x), collapse=","), character(1)))) != 1L)
      stop("V_i must contain equally sized padded matrices.")
    dat$V_i <- array(NA_real_, c(length(vv), nrow(vv[[1]]), ncol(vv[[1]])))
    for (s in seq_along(vv)) dat$V_i[s,,] <- vv[[s]]
  }
  .ctmaLohmannValidateData(dat, n.latent)
  note <- paste("Literal supplement replication, not a validated stationary CT model:",
    "the original residual-covariance formula is retained;",
    "negative drift diagonals do not guarantee stability or valid diffusion.",
    "FE additionally uses unbounded correlation coordinates and N(0,100) priors.",
    "Inspect inst/research/ctmaLohmannFit.md before interpreting results.")
  if (n.latent == 1L) note <- paste(
    "Univariate extension of the Lohmann approach, not an original-supplement replication.",
    "Uses a stationary unit-total-variance process with random drift and trait proportion.",
    "The URE/ERE/CRE priors follow the corresponding bivariate conventions.")
  warning(note, call. = FALSE)
  out <- list(call = match.call(), reStructure = reStructure, n.latent = n.latent, stanData = dat,
    modelCode = source$code, source = source[c("path", "md5")],
    notes = note, fit = NULL, summary = NULL,
    parameter_order = c("A11_raw", "A21", "A12", "A22_raw",
      "trait_X_logit", "trait_Y_logit", "trait_correlation_raw", "within_correlation_raw"))
  if (n.latent == 1L) out$parameter_order <- c("A11_raw", "trait_X_logit")
  if (fit) {
    scalarInt <- function(x, name, lower) {
      if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
          x < lower || x != floor(x) || x > .Machine$integer.max)
        stop(name, " must be an integer >= ", lower, ".")
      as.integer(x)
    }
    chains <- scalarInt(chains, "chains", 1)
    iter <- scalarInt(iter, "iter", 2)
    warmup <- scalarInt(warmup, "warmup", 0)
    coresToUse <- scalarInt(coresToUse, "coresToUse", 1)
    seed <- scalarInt(seed, "seed", 1)
    if (warmup >= iter) stop("warmup must be smaller than iter.")
    if (!is.list(control)) stop("control must be a list of Stan NUTS controls.")
    if (!requireNamespace("rstan", quietly = TRUE)) stop("Install rstan to fit this model.")
    if (is.null(stanModel)) {
      stanModel <- rstan::stan_model(model_code = source$code,
        model_name = paste0("CoTiMA2_Lohmann_", n.latent, "D_", reStructure), auto_write = FALSE)
    } else if (!methods::is(stanModel, "stanmodel") ||
        !identical(stanModel@model_code, source$code))
      stop("stanModel must be compiled from this exact supplement modelCode.")
    message("Lohmann ", reStructure, ": sampling the correlation likelihood with Stan/NUTS.")
    out$fit <- rstan::sampling(stanModel, data = dat, chains = chains,
      iter = iter, warmup = warmup, cores = min(coresToUse, chains),
      seed = seed, control = control, init = init)
    out$summary <- rstan::summary(out$fit)$summary
  }
  class(out) <- "ctmaLohmannFit"
  if (!is.null(activeDirectory)) {
    if (length(activeDirectory) != 1L || !dir.exists(activeDirectory))
      stop("activeDirectory must be an existing directory.")
    saveRDS(out, file.path(activeDirectory, paste0("ctmaLohmannFit_", n.latent, "D_", reStructure, ".rds")))
  }
  out
}

.ctmaLohmannSource <- function(directory, structure) {
  if (!is.character(directory) || length(directory) != 1L || is.na(directory))
    stop("Supply supplementDirectory: the downloaded original supplement folder or modelcompiler.R file.")
  files <- if (file.exists(directory) && !dir.exists(directory)) directory else
    list.files(directory, pattern = "^modelcompiler[.]R$", recursive = TRUE, full.names = TRUE)
  if (!length(files)) stop("No modelcompiler.R found in supplementDirectory.")
  hashes <- unname(tools::md5sum(files))
  if (length(unique(hashes)) != 1L) stop("Different compiler files found; supply one explicit modelcompiler.R path.")
  expressions <- parse(files[1], keep.source = FALSE)
  key <- paste0("OS_", structure, "_char")
  values <- lapply(expressions, function(x) {
    if (is.call(x) && as.character(x[[1]]) %in% c("=", "<-") &&
        identical(x[[2]], as.name(key)) && is.character(x[[3]]) && length(x[[3]]) == 1L) x[[3]] else NULL
  })
  values <- Filter(Negate(is.null), values)
  if (length(values) != 1L) stop("Expected exactly one literal Stan string named ", key, ".")
  list(code = values[[1]], path = normalizePath(files[1]), md5 = hashes[1])
}

.ctmaLohmannData <- function(matrices, n, lag, V, n.latent = 2L) {
  if (!is.list(matrices) || !length(matrices)) stop("Supply a nonempty list of correlationMatrices.")
  k <- length(matrices)
  if (!is.numeric(n) || length(n) != k || any(!is.finite(n)) || any(n <= 1))
    stop("sampleSizes must contain one finite N > 1 per study.")
  if (!is.numeric(lag) || length(lag) != k || any(!is.finite(lag)) || any(lag <= 0))
    stop("timeLag must contain one positive finite interval per study.")
  for (R in matrices) {
    if (!is.matrix(R) || !is.numeric(R) || nrow(R) != ncol(R) ||
        nrow(R) < 2L * n.latent || nrow(R) %% n.latent || any(!is.finite(R)) ||
        !isTRUE(all.equal(R, t(R), check.attributes = FALSE)) ||
        any(abs(diag(R)-1) > 1e-8) || any(abs(R) > 1 + 1e-8))
      stop("Use complete symmetric correlation matrices, diagonal 1, in wave-major order (X1,X2,... for one variable; X1,Y1,X2,Y2,... for two).")
  }
  if (!is.list(V) || length(V) != k)
    stop("Supply samplingCovariances, one known asymptotic correlation covariance matrix per study. No automatic approximation is chosen.")
  dims <- vapply(matrices, nrow, integer(1))
  nc <- dims * (dims-1L)/2L
  m <- max(nc)
  ri <- matrix(-99, k, m)
  vi <- array(-99, c(k,m,m))
  for (s in seq_len(k)) {
    if (!is.matrix(V[[s]]) || !is.numeric(V[[s]]) ||
        !identical(dim(V[[s]]), rep(as.integer(nc[s]),2L)))
      stop("samplingCovariances dimensions do not match the observed lower triangle.")
    ri[s,seq_len(nc[s])] <- matrices[[s]][lower.tri(matrices[[s]])]
    vi[s,seq_len(nc[s]),seq_len(nc[s])] <- V[[s]]
  }
  list(Nstudies=k, maxwaves=max(dims)/n.latent, Nlatent=max(dims)+n.latent,
    Nmanifest=max(dims), max_cors=m, n_cors=as.integer(nc),
    Nmeasurements=as.integer(dims), TimeLag=lag, r_i=ri, V_i=vi)
}

.ctmaLohmannValidateData <- function(d, n.latent = 2L) {
  required <- c("Nstudies","maxwaves","Nlatent","Nmanifest","max_cors",
    "n_cors","Nmeasurements","TimeLag","r_i","V_i")
  if (!is.list(d) || !all(required %in% names(d))) stop("Incomplete stanData.")
  for (key in required[1:5]) {
    x <- d[[key]]
    if (!is.numeric(x) || length(x)!=1L || !is.finite(x) || x<1 || x!=floor(x))
      stop("Invalid dimension: ", key)
  }
  k <- d$Nstudies; m <- d$max_cors
  if (d$maxwaves < 2 || d$Nmanifest != n.latent*d$maxwaves ||
      d$Nlatent != d$Nmanifest+n.latent || m != choose(d$Nmanifest,2)) stop("Inconsistent Stan dimensions.")
  for (key in c("n_cors","Nmeasurements","TimeLag"))
    if (!is.numeric(d[[key]]) || length(d[[key]])!=k || any(!is.finite(d[[key]]))) stop("Invalid ",key)
  if (any(d$Nmeasurements < 2*n.latent | d$Nmeasurements > d$Nmanifest |
      d$Nmeasurements %% n.latent != 0) || any(d$n_cors != choose(d$Nmeasurements,2)) ||
      any(d$TimeLag <= 0) || !is.matrix(d$r_i) ||
      !identical(dim(d$r_i), as.integer(c(k,m))) ||
      !identical(dim(d$V_i), as.integer(c(k,m,m)))) stop("Inconsistent study dimensions or intervals.")
  for (s in seq_len(k)) {
    idx <- seq_len(d$n_cors[s]); r <- d$r_i[s,idx]; v <- matrix(d$V_i[s,idx,idx], length(idx), length(idx))
    if (any(!is.finite(r)) || any(abs(r)>1) || any(!is.finite(v)) ||
        !isTRUE(all.equal(v,t(v),check.attributes=FALSE)) ||
        inherits(try(chol(v),silent=TRUE),"try-error")) stop("Invalid correlations or non-positive-definite sampling covariance in study ",s)
  }
  invisible(TRUE)
}

# Authored univariate extension; the original supplement is bivariate only.
# One stable OU within-person process plus an independent time-constant trait.
.ctmaLohmannUnivariateCode <- function(structure) {
  stopifnot(structure %in% c("FE", "URE", "ERE", "CRE"))
  parameters <- switch(structure,
    FE = "",
    URE = "vector<lower=0>[2] Theta_SD; array[Nstudies] vector[2] Theta_s;",
    ERE = "vector<lower=0>[2] Theta_SD; real cor_raw; matrix[2,Nstudies] z_Theta;",
    CRE = "vector<lower=0>[2] Theta_SD; cholesky_factor_corr[2] V_Theta; matrix[2,Nstudies] z_Theta;")
  transformed <- switch(structure,
    FE = "matrix[2,Nstudies] deviations = rep_matrix(0,2,Nstudies);",
    URE = "matrix[2,Nstudies] deviations; for (s in 1:Nstudies) deviations[,s] = Theta_s[s];",
    ERE = paste("matrix[2,2] Omega; matrix[2,Nstudies] deviations;",
      "real corRE = 2*inv_logit(cor_raw)-1; Omega = rep_matrix(corRE,2,2);",
      "Omega[1,1]=1; Omega[2,2]=1;",
      "deviations=diag_pre_multiply(Theta_SD,cholesky_decompose(Omega))*z_Theta;"),
    CRE = "matrix[2,Nstudies] deviations = diag_pre_multiply(Theta_SD,V_Theta)*z_Theta;")
  priors <- switch(structure,
    FE = "Theta_AR ~ normal(0,100); Theta_RI ~ normal(0,100);",
    URE = paste("Theta_AR ~ normal(0,1); Theta_RI ~ normal(0,1); Theta_SD ~ cauchy(0,2.5);",
      "for (s in 1:Nstudies) Theta_s[s] ~ normal(rep_vector(0,2),Theta_SD);"),
    ERE = paste("Theta_AR ~ normal(0,1); Theta_RI ~ normal(0,1); Theta_SD ~ cauchy(0,2.5);",
      "to_vector(z_Theta) ~ normal(0,1); cor_raw ~ normal(0,1);"),
    CRE = paste("Theta_AR ~ normal(0,1); Theta_RI ~ normal(0,1); Theta_SD ~ cauchy(0,2.5);",
      "to_vector(z_Theta) ~ normal(0,1); V_Theta ~ lkj_corr_cholesky(2);"))
  paste0("// CoTiMA2 univariate extension, not original Lohmann code.\n",
    "data { int<lower=1> Nstudies; int<lower=2> maxwaves; int<lower=1> max_cors;\n",
    "array[Nstudies] int<lower=1> n_cors; array[Nstudies] int<lower=2> Nmeasurements;\n",
    "vector<lower=0>[Nstudies] TimeLag; matrix[Nstudies,max_cors] r_i;\n",
    "array[Nstudies] matrix[max_cors,max_cors] V_i; }\n",
    "parameters { real Theta_AR; real Theta_RI; ",parameters," }\n",
    "transformed parameters { ",transformed,"\n",
    "vector[Nstudies] drift; vector[Nstudies] trait; vector[Nstudies] diffusion;\n",
    "for (s in 1:Nstudies) { drift[s]=-exp(Theta_AR+deviations[1,s]);\n",
    "trait[s]=inv_logit(Theta_RI+deviations[2,s]);\n",
    "diffusion[s]=-2*drift[s]*(1-trait[s]); } }\n",
    "model { ",priors,"\n",
    "for (s in 1:Nstudies) { int nc=n_cors[s]; vector[nc] rho; int idx=1;\n",
    "for (j in 1:(Nmeasurements[s]-1)) for (i in (j+1):Nmeasurements[s]) {\n",
    "rho[idx]=trait[s]+(1-trait[s])*exp(drift[s]*TimeLag[s]*(i-j)); idx+=1; }\n",
    "r_i[s,1:nc]' ~ multi_normal(rho,V_i[s][1:nc,1:nc]); } }\n")
}
