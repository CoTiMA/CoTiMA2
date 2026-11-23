# One joint posterior draw per row and one fitted study profile per column.
.ctmaExtractNativeDrift <- function(fe, n_studies, cores = 2L) {
  fit <- fe[["studyFitList"]]
  if (!inherits(fit, "ctJuliaFit"))
    stop("Expected a pooled native ctsem fit in fe$studyFitList.")
  raw <- fit$estimate$rawposterior
  if (is.null(fit$sample) || !is.matrix(raw) || nrow(raw) < 2L)
    stop("Expected joint posterior draws from a sampled fit.")

  ns <- asNamespace("ctsem")
  model <- get(".ctFitModelObject", ns)(fit)
  ids <- get(".ctFitIdMap", ns)(fit)$original
  ti <- get(".ctFitTIpredData", ns)(fit)
  k <- as.integer(n_studies)
  if (length(k) != 1L || is.na(k) || k < 2L ||
      ncol(ti) != k - 1L || nrow(ti) != length(ids) ||
      any(!is.finite(ti)) || anyNA(ids) || anyDuplicated(ids) ||
      !identical(colnames(ti), model$TIpredNames))
    stop("Expected exactly K - 1 study contrasts and no extra TI predictors.")

  first <- which(!duplicated(as.data.frame(ti)))
  if (length(first) != k)
    stop("The number of distinct fitted contrast profiles differs from K.")

  # ctExtract evaluates selected subjects in fitted internal ID order.
  selected <- ctsem::ctExtract(fit, subjectMatrices = TRUE,
    subjects = first, realid = FALSE, nsamples = "all", cores = cores)
  drift <- selected[["subj_DRIFT"]]
  if (!identical(dim(drift), c(nrow(raw), k, 1L, 1L)) ||
      !isTRUE(all.equal(unname(selected$rawpars), unname(raw),
                        tolerance = 0, check.attributes = FALSE)) ||
      any(!is.finite(drift)))
    stop("Expected all joint draws and one univariate DRIFT per study profile.")

  profiles <- paste0("profile_", seq_len(k))
  ans <- matrix(drift, nrow(raw), k, dimnames = list(NULL, profiles))
  attr(ans, "study_profiles") <- data.frame(
    profile = profiles, representative_id = ids[first],
    ti[first, , drop = FALSE], check.names = FALSE)
  ans
}
