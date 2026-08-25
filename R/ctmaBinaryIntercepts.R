# Validate and normalize user-supplied location parameter specifications.
.ctma_normalize_location_specs <- function(cint, manifestMeans,
    nlatent, nmanifest) {
  nlatent <- as.integer(nlatent)[1L]
  nmanifest <- as.integer(nmanifest)[1L]
  if (is.na(nlatent) || nlatent < 1L || is.na(nmanifest) || nmanifest < 1L) {
    stop("Latent and manifest dimensions must be positive integers.",
      call. = FALSE)
  }
  # Keep scalar zero scalar because ctmaLabels handles that public shorthand;
  # the model constructors expand it to the required matrix dimensions.
  if (is.null(cint)) cint <- 0
  if (is.null(manifestMeans)) manifestMeans <- 0
  if (length(cint) > 1L && length(cint) != nlatent) stop(
    "You specified CINTs, but their length does not correspond to the ",
    "number of latent variables.", call. = FALSE)
  if (length(manifestMeans) > 1L && length(manifestMeans) != nmanifest) stop(
    "You specified manifestMeans, but their length does not correspond to ",
    "the number of manifest variables.", call. = FALSE)
  list(cint = cint, manifestMeans = manifestMeans)
}

# Validate location parameters for binary manifest variables using the final
# ctsem parameter table. A CINT counts only when it is connected to the
# manifest variable through an active loading.
.ctma_check_binary_intercepts <- function(model, binaries,
    process_intercepts = FALSE, inform = TRUE) {
  if (is.null(binaries)) return(invisible(model))
  if (is.null(model$pars) || !is.data.frame(model$pars)) stop(
    "The ctsem model has no parameter table for binary-intercept validation.",
    call. = FALSE)

  nmanifest <- as.integer(model$n.manifest)[1L]
  nlatent <- as.integer(model$n.latent)[1L]
  if (!is.finite(nmanifest) || nmanifest < 1L ||
      !is.finite(nlatent) || nlatent < 1L) stop(
    "The ctsem model does not report valid manifest and latent dimensions.",
    call. = FALSE)
  if (length(binaries) != nmanifest) stop(
    "binaries must contain exactly one entry per manifest variable.",
    call. = FALSE)
  binary_manifest <- which(as.integer(binaries) == 1L)
  if (!length(binary_manifest) || isTRUE(process_intercepts))
    return(invisible(model))

  pars <- model$pars
  free_by_row <- function(matrix_name, nr) {
    answer <- rep(FALSE, nr)
    selected <- pars[!is.na(pars$matrix) & pars$matrix == matrix_name,
      , drop = FALSE]
    if (!nrow(selected)) return(answer)
    for (i in seq_len(nrow(selected))) {
      r <- as.integer(selected$row[i])
      if (!is.na(r) && r >= 1L && r <= nr)
        answer[r] <- answer[r] || !is.na(selected$param[i])
    }
    answer
  }

  manifest_mean_free <- free_by_row("MANIFESTMEANS", nmanifest)
  cint_free <- free_by_row("CINT", nlatent)
  lambda_active <- matrix(FALSE, nmanifest, nlatent)
  lambda <- pars[!is.na(pars$matrix) & pars$matrix == "LAMBDA",
    , drop = FALSE]
  if (nrow(lambda)) for (i in seq_len(nrow(lambda))) {
    r <- as.integer(lambda$row[i]); cc <- as.integer(lambda$col[i])
    if (is.na(r) || is.na(cc) || r < 1L || r > nmanifest ||
        cc < 1L || cc > nlatent) next
    is_free <- !is.na(lambda$param[i])
    is_fixed_nonzero <- !is.na(lambda$value[i]) && lambda$value[i] != 0
    lambda_active[r, cc] <- is_free || is_fixed_nonzero
  }
  cint_support <- vapply(seq_len(nmanifest), function(j)
    any(lambda_active[j, ] & cint_free), logical(1L))

  has_cint <- cint_support[binary_manifest]
  has_manifest_mean <- manifest_mean_free[binary_manifest]
  labels <- model$manifestNames
  if (is.null(labels) || length(labels) < nmanifest)
    labels <- paste0("manifest ", seq_len(nmanifest))
  describe <- function(index) paste(labels[index], collapse = ", ")

  missing <- !has_cint & !has_manifest_mean
  if (any(missing)) stop(
    "Binary variable(s) ", describe(binary_manifest[missing]),
    " have neither a freely estimated CINT nor a freely estimated ",
    "MANIFESTMEANS parameter. Set the corresponding CINT or manifestMeans ",
    "parameter free (for example with 'auto').",
    call. = FALSE)

  manifest_only <- !has_cint & has_manifest_mean
  if (isTRUE(inform) && any(manifest_only)) message(
    "Binary variable(s) ", describe(binary_manifest[manifest_only]),
    " use freely estimated MANIFESTMEANS while their connected latent CINTs ",
    "are fixed.")

  both <- has_cint & has_manifest_mean
  if (isTRUE(inform) && any(both)) warning(
    "Binary variable(s) ", describe(binary_manifest[both]),
    " have both freely estimated MANIFESTMEANS and connected CINT ",
    "parameters. Check location identification; normally one should be fixed.",
    call. = FALSE)
  invisible(model)
}
