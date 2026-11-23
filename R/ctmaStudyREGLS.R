# Internal second-stage meta-analysis of one transformed parameter.
.ctmaStudyREGLS <- function(estimates, covariance, level = 0.95, scaleTime = 1) {
  if (!is.numeric(scaleTime) || length(scaleTime) != 1L ||
      !is.finite(scaleTime) || scaleTime <= 0)
    stop("scaleTime must be one finite positive number.", call. = FALSE)
  k <- length(estimates)
  if (!is.numeric(estimates) || k < 2L || any(!is.finite(estimates)) ||
      !is.matrix(covariance) || !is.numeric(covariance) ||
      !identical(dim(covariance), c(k, k)) || any(!is.finite(covariance)) ||
      !isSymmetric(unname(covariance)))
    stop("Supply at least two finite study estimates and a matching symmetric covariance matrix.")
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) stop("level must be between zero and one.")
  if (!is.null(names(estimates)) && !is.null(dimnames(covariance)) &&
      (!identical(names(estimates), rownames(covariance)) ||
       !identical(names(estimates), colnames(covariance))))
    stop("Study estimate names and covariance labels must match in order.")
  ones <- rep(1, k)
  gls <- function(sigma, tau2) {
    r <- chol(sigma + diag(tau2, k))
    a <- forwardsolve(t(r), ones)
    b <- forwardsolve(t(r), estimates)
    precision <- sum(a^2)
    mu <- sum(a*b)/precision
    q <- sum((b - mu*a)^2)
    list(mean = mu, se = 1/sqrt(precision), q = q,
         reml = (2*sum(log(diag(r))) + log(precision) + q)/2)
  }
  full <- unname(covariance)
  r <- tryCatch(chol(full), error = function(e)
    stop("The study covariance must be positive definite.", call. = FALSE))
  diagonal <- diag(diag(full), k)
  df <- k - 1L
  upper <- sum((estimates - mean(estimates))^2)/df
  moments <- function(sigma, r) {
    w <- chol2inv(r)
    u <- rowSums(w)
    q <- gls(sigma, 0)$q
    dl <- max(0, (q - df)/(sum(diag(w)) - sum(u^2)/sum(u)))
    pm <- if (q <= df) 0 else stats::uniroot(
      function(tau2) gls(sigma, tau2)$q - df,
      c(0, upper), tol = max(.Machine$double.xmin, upper*1e-12)
    )$root
    c(dl = dl, pm = pm)
  }
  f <- moments(full, r)
  d <- moments(diagonal, chol(diagonal))
  # Search log-spaced intervals to allow multiple restricted-likelihood minima.
  objective <- function(tau2) gls(full, tau2)$reml
  # In the residual contrast space, the likelihood is increasing beyond
  # the squared length of the centered estimates.
  bound <- upper*df
  reml <- 0
  if (bound > 0) {
    grid <- c(0, bound*10^seq(-12, 0, length.out = 97L))
    values <- vapply(grid, objective, numeric(1))
    minima <- which(values[-c(1L, length(grid))] < values[-c(length(grid)-1L, length(grid))] &
                    values[-c(1L, length(grid))] <= values[-c(1L, 2L)]) + 1L
    intervals <- unique(c(2L, minima, length(grid)-1L))
    candidates <- vapply(intervals, function(i)
      stats::optimize(objective, grid[c(i-1L, i+1L)],
                      tol = max(.Machine$double.xmin, grid[i+1L]*1e-10))$minimum,
      numeric(1))
    candidates <- c(0, candidates, bound)
    reml <- candidates[which.min(vapply(candidates, objective, numeric(1)))]
  }
  tau2 <- c(dl_full = f[["dl"]], dl_diag = d[["dl"]],
            pm_full = f[["pm"]], pm_diag = d[["pm"]], reml_full = reml)
  critical <- stats::qnorm((1 + level) / 2)
  pooled_by_method <- do.call(rbind, lapply(names(tau2), function(method) {
    sigma <- if (endsWith(method, "_diag")) diagonal else full
    fit <- gls(sigma, tau2[[method]])
    half_width <- critical * fit$se
    c(mean = fit$mean, se = fit$se,
      lower = fit$mean - half_width, upper = fit$mean + half_width)
  }))
  rownames(pooled_by_method) <- names(tau2)
  pooled <- pooled_by_method["pm_full", ]
  result <- list(tau2 = tau2, fixed_mean = gls(full, 0)$mean,
       pooled = pooled,
       pooled_by_method = data.frame(tau2 = tau2, pooled_by_method),
       q_zero = c(full = gls(full, 0)$q, diagonal = gls(diagonal, 0)$q),
       q_pooled = gls(full, tau2[["pm_full"]])$q)
  if (scaleTime != 1) {
    original <- result[c("tau2", "fixed_mean", "pooled", "pooled_by_method")]
    original$tau2 <- original$tau2 * scaleTime^2
    original$fixed_mean <- original$fixed_mean * scaleTime
    original$pooled <- original$pooled * scaleTime
    original$pooled_by_method$tau2 <- original$tau2
    columns <- c("mean", "se", "lower", "upper")
    original$pooled_by_method[columns] <-
      original$pooled_by_method[columns] * scaleTime
    original$scaleTime <- scaleTime
    result$originalTime <- original
  }
  structure(result, class = "ctmaStudyREGLS")
}

summary.ctmaStudyREGLS <- function(object, ...) {
  structure(object[intersect(c("pooled_by_method", "pooled", "q_pooled", "originalTime"),
                             names(object))],
            class = "ctmaStudyREGLS")
}

print.ctmaStudyREGLS <- function(x, ...) {
  view <- unclass(x)
  view$originalTime <- NULL
  print(view, ...)
  if (!is.null(x$originalTime)) cat(
    "Original-time GLS results are available under $originalTime.\n")
  cat("\nDL = DerSimonian\u2013Laird; PM = Paule\u2013Mandel; ",
      "REML = Restricted Maximum Likelihood; ",
      "diag = diagonal covariance; full = full covariance.\n", sep = "")
  invisible(x)
}
