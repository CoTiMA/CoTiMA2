# Optimize the sum of independent person-level objectives, then let ctsem
# refine the solution and compute uncertainty on the complete dataset.
.ctmaNativeLikelihoodFit <- function(args, compute) {
  if (!isTRUE(args$optimize) || isTRUE(args$priors))
    stop("Distributed native likelihood requires ML without priors.")
  if (any(vapply(compute$remote, function(h) h$authentication != "auto", logical(1))))
    stop("Distributed native likelihood requires noninteractive SSH authentication.")
  prepare <- args
  prepare$fit <- FALSE
  full <- .ctmaNativeCtsemJuliaFit(prepare)
  levels <- full$laplace$levels
  if (length(levels) && (length(levels) != 1L || levels[[1L]]$name != "id"))
    stop("Distributed likelihood currently supports independent person-level random effects.")
  if (length(full$ti_missing$parameter) > 0L)
    stop("Distributed likelihood requires complete time-independent predictors.")

  if (!is.null(full$laplace) && length(levels) && !identical(full$intoverpop, "laplace"))
    stop("Distributed native likelihood requires Laplace integration for person-level random effects.")

  profiles <- c(list(list(host = "localhost", workers = compute$local,
                         threads = compute$threads)), compute$remote)
  ids <- unique(args$datalong$id)
  slots <- rep(seq_along(profiles), vapply(profiles, function(h) h$workers, integer(1)))
  slots <- slots[order(stats::ave(seq_along(slots), slots, FUN = seq_along), slots)]
  slots <- sort(utils::head(slots, length(ids)))
  if (length(slots) < 2L) stop("Distributed likelihood requires at least two workers and persons.")
  groups <- split(ids, rep(seq_along(slots), length.out = length(ids)))
  layout <- function(spec) {
    list(spec$parameter_table, spec$ti_effects,
         lapply(spec$laplace$levels, function(level)
           level[setdiff(names(level), c("group", "ngroups", "labels"))]))
  }
  shards <- lapply(groups, function(group) {
    prepare$datalong <- args$datalong[args$datalong$id %in% group, , drop = FALSE]
    shard <- .ctmaNativeCtsemJuliaFit(prepare)
    if (!identical(layout(shard), layout(full)))
      stop("A worker partition changed the native parameter layout.")
    shard
  })
  npar <- max(c(0L, full$parameter_table$parnumber, full$laplace$npar,
                full$ti_effects$coefficient), na.rm = TRUE)
  if (!is.finite(npar) || npar < 1L) stop("No free native parameters to optimize.")
  start <- args$inits
  derived <- FALSE
  if (is.null(start)) {
    start <- numeric(npar)
    data_start <- get0(".ctDataStart", asNamespace("ctsem"), inherits = FALSE)
    if (is.function(data_start) && !identical(args$optimcontrol$datastart, FALSE)) {
      candidate <- tryCatch(data_start(args$datalong, args$model, full, npar),
                            error = function(e) NULL)
      if (is.numeric(candidate) && length(candidate) == npar) {
        usable <- is.finite(candidate)
        start[usable] <- candidate[usable]
        derived <- any(usable)
      }
    }
  }
  if (!is.numeric(start) || length(start) != npar || any(!is.finite(start)))
    stop("Supply a finite native raw starting vector matching the full model.")
  evaluate <- function(spec, raw, cores) {
    out <- getExportedValue("ctsem", "ctJuliaEvaluate")(spec, raw, cores = cores)
    value <- c(out$value, out$gradient)
    if (length(value) != length(raw) + 1L)
      stop("Incompatible native objective evaluation.")
    value
  }
  environment(evaluate) <- baseenv()
  cluster <- structure(list(), class = c("SOCKcluster", "cluster"))
  on.exit(if (length(cluster)) parallel::stopCluster(cluster), add = TRUE)
  for (i in unique(slots)) {
    h <- profiles[[i]]
    opts <- list(workers = rep(h$host, sum(slots == i)), useXDR = TRUE,
                 rscript_envs = c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
                                  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"))
    if (i == 1L) opts$rscript_libs <- .libPaths() else {
      opts$user <- h$user
      opts$rscript <- h$rscript
      opts$revtunnel <- TRUE
      opts$rshopts <- c("-p", h$port, "-o", "BatchMode=yes", "-o", "ExitOnForwardFailure=yes")
    }
    nodes <- do.call(parallelly::makeClusterPSOCK, opts)
    cluster <- structure(c(unclass(cluster), unclass(nodes)), class = class(cluster))
  }
  cores <- lapply(profiles[slots], function(h) if (is.null(h$threads)) compute$threads else h$threads)
  setup <- function(spec, cores) {
    assign(".cotima_likelihood", list(spec = spec, cores = cores), .GlobalEnv)
    d <- utils::packageDescription("ctsem")
    c(d$Version, d$RemoteSha)
  }
  environment(setup) <- baseenv()
  versions <- parallel::clusterMap(cluster, setup, shards, cores, RECYCLE = FALSE)
  d <- utils::packageDescription("ctsem")
  if (is.null(d$RemoteSha) || !all(vapply(versions, identical, logical(1), c(d$Version, d$RemoteSha))))
    stop("Distributed likelihood requires matching, Git-pinned ctsem installations.")
  worker <- function(raw, evaluate) {
    state <- get(".cotima_likelihood", .GlobalEnv)
    evaluate(state$spec, raw, state$cores)
  }
  environment(worker) <- baseenv()
  previous <- value <- NULL
  combined <- function(raw) {
    if (!identical(raw, previous)) {
      value <<- Reduce(`+`, parallel::clusterCall(cluster, worker, raw, evaluate))
      previous <<- raw
    }
    value
  }
  check <- function(raw) {
    reference <- evaluate(full, raw, args$cores)
    delta <- abs(reference - combined(raw))
    if (any(!is.finite(delta)) || any(delta > 1e-7 + 1e-9 * abs(reference)))
      stop("Distributed likelihood or gradient disagrees with the full dataset.")
    reference
  }
  if (derived && any(!is.finite(evaluate(full, start, args$cores)))) {
    start <- numeric(npar)
    derived <- FALSE
  }
  baseline <- check(start)
  check(start + seq(-0.1, 0.1, length.out = length(start)))
  message("Optimizing native likelihood across ", length(slots), " workers.")
  timer <- proc.time()[["elapsed"]]
  opt <- stats::optim(start, function(raw) {
    value <- combined(raw)
    if (all(is.finite(value))) -value[1L] else 1e100
  }, function(raw) {
    value <- combined(raw)
    if (all(is.finite(value))) -value[-1L] else numeric(length(raw))
  }, method = "BFGS",
                     control = list(maxit = 200L, reltol = 1e-10))
  elapsed <- proc.time()[["elapsed"]] - timer
  checked <- check(opt$par)
  parallel::stopCluster(cluster)
  cluster <- list()
  rejection <- .ctmaNativeLikelihoodRejection(opt$convergence, baseline, checked)
  if (is.null(rejection)) {
    message("Refining the distributed solution on the full dataset.")
    args$inits <- opt$par
  } else {
    warning("Distributed candidate rejected (", rejection,
            "). Continuing with the native ctsem starting strategy.", call. = FALSE)
  }
  result <- .ctmaNativeCtsemJuliaFit(args)
  result$distributed_likelihood <- list(workers = length(slots), seconds = elapsed,
    accepted = is.null(rejection), rejection = rejection,
    start = if (derived) "data" else if (is.null(prepare$inits)) "zero" else "supplied",
    start_objective = baseline[1L], optimizer_code = opt$convergence, objective = checked[1L],
    max_gradient = max(abs(checked[-1L])))
  result
}


.ctmaNativeLikelihoodRejection <- function(code, baseline, candidate) {
  if (code != 0L) return("optimizer did not converge")
  if (candidate[1L] < baseline[1L] - 1e-7) return("objective worsened")
  if (max(abs(candidate[-1L])) > max(0.01, 0.1 * max(abs(baseline[-1L]))))
    return("gradient remains large relative to the start")
  NULL
}
