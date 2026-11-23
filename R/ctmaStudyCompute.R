.ctmaStudyWorkerFit <- function(args, cores) {
  args$cores <- cores
  notes <- character()
  fit <- tryCatch(withCallingHandlers(
    do.call(CoTiMA2:::.ctmaInitCtsemJuliaStudyFit, args),
    warning = function(w) {
      notes <<- c(notes, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ), error = function(e) e)
  list(fit = fit, warnings = notes, host = Sys.info()[["nodename"]])
}

.ctmaRunStudyCompute <- function(requests, compute, studies) {
  ids <- which(!vapply(requests, is.null, logical(1)))
  result <- vector("list", length(requests))
  if (!length(ids)) return(result)
  if (!requireNamespace("parallelly", quietly = TRUE))
    stop("Study scheduling requires the parallelly package.")
  profiles <- c(list(list(host = "localhost", workers = compute$local,
                          threads = compute$threads)), compute$remote)
  # Round robin allocation so short queues can use more than one machine.
  hosts <- rep(seq_along(profiles), vapply(profiles, function(h) h$workers, integer(1)))
  slots <- order(stats::ave(seq_along(hosts), hosts, FUN = seq_along), hosts)
  counts <- tabulate(hosts[utils::head(slots, length(ids))], nbins = length(profiles))
  cluster <- structure(list(), class = c("SOCKcluster", "cluster"))
  on.exit(if (length(cluster)) parallel::stopCluster(cluster), add = TRUE)
  for (i in which(counts > 0L)) {
    h <- profiles[[i]]
    cores <- if (is.null(h$threads)) compute$threads else h$threads
    options <- list(workers = rep(h$host, counts[i]), useXDR = TRUE,
      rscript_envs = c(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
                       MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"))
    if (i == 1L) {
      options$rscript_libs <- .libPaths()
    } else {
      options$user <- h$user
      options$rscript <- h$rscript
      options$revtunnel <- TRUE
      options$rshopts <- c("-p", as.character(h$port),
        "-o", "ExitOnForwardFailure=yes",
        "-o", if (h$authentication == "auto") "BatchMode=yes" else "BatchMode=no")
    }
    nodes <- do.call(parallelly::makeClusterPSOCK, options)
    cluster <- structure(c(unclass(cluster), unclass(nodes)), class = class(cluster))
    initialize <- function(cores) {
      library(CoTiMA2)
      assign(".ctmaStudyCores", cores, envir = .GlobalEnv)
      list(cotima = as.character(utils::packageVersion("CoTiMA2")),
           ctsem = as.character(utils::packageVersion("ctsem")),
           sha = utils::packageDescription("ctsem")$RemoteSha,
           available = CoTiMA2:::.ctmaCtsemJuliaLaplaceStatus()$available)
    }
    environment(initialize) <- baseenv()
    versions <- parallel::clusterCall(nodes, initialize, cores)
    expected_sha <- utils::packageDescription("ctsem")$RemoteSha
    if (!all(vapply(versions, function(v)
      isTRUE(v$available) &&
      identical(v$cotima, as.character(utils::packageVersion("CoTiMA2"))) &&
      identical(v$ctsem, as.character(utils::packageVersion("ctsem"))) &&
      (is.null(expected_sha) || identical(v$sha, expected_sha)), logical(1))))
      stop("Study worker package versions or native ctsem API differ on ", h$host)
  }

  run <- function(args) {
    CoTiMA2:::.ctmaStudyWorkerFit(args, get(".ctmaStudyCores", envir = .GlobalEnv))
  }
  environment(run) <- baseenv()
  fitted <- parallel::clusterApplyLB(cluster, requests[ids], run)
  for (j in seq_along(ids)) {
    i <- ids[j]
    label <- .ctmaInitStudyWarningLabel(i, length(requests), studies[[i]])
    for (note in fitted[[j]]$warnings) warning(paste0(label, note), call. = FALSE)
    value <- fitted[[j]]$fit
    if (inherits(value, "error"))
      stop(label, fitted[[j]]$host, ": ", conditionMessage(value), call. = FALSE)
    result[[i]] <- value
  }
  result
}
