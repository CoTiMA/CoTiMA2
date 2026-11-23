# Add only explicitly selected population correlations.
.ctmaTablesRIRequested <- function(parameters) {
  any(grepl("^ri(t0)?cor($|\\[)", parameters, ignore.case = TRUE))
}

.ctmaTablesAddRI <- function(d, fit, stored, parameters) {
  if (!.ctmaTablesRIRequested(parameters)) return(d)
  native <- .ctmaTablesNative(fit)
  pars <- fit$ctModel$pars
  if (is.null(pars)) pars <- native$ctstanmodelbase$pars
  if (is.null(pars)) pars <- native$model$pars
  if (!all(c("matrix", "param", "indvarying") %in% names(pars)))
    .ctmaTablesStop("Random-intercept correlations require the saved model parameter specification.")

  ri <- stored$randomIntercepts
  saved <- if (is.list(ri)) ri$popcor_mean else NULL
  engine <- stored
  if (is.null(engine$rawpopcorr) && !is.null(native)) {
    engine <- native$resultsSummary
    if (is.null(engine$rawpopcorr)) engine <- .ctmaTablesSummary(native)
  }
  labels <- if (is.matrix(saved)) rownames(ri$popsd) else rownames(engine$popsd)
  if (is.null(labels) || anyDuplicated(labels) || length(labels) < 2L)
    .ctmaTablesStop("No named population random-effect correlation summaries are available for this fit.")

  varying <- pars[pars$indvarying %in% TRUE, , drop = FALSE]
  kind <- toupper(varying$matrix[match(labels, varying$param)])
  intercept <- kind %in% c("MANIFESTMEANS", "CINT")
  t0 <- kind %in% "T0MEANS"
  pair <- which(lower.tri(matrix(FALSE, length(labels), length(labels))), arr.ind = TRUE)
  bothRI <- intercept[pair[,1L]] & intercept[pair[,2L]]
  withT0 <- (intercept | t0)[pair[,1L]] & (intercept | t0)[pair[,2L]] & !bothRI
  pair <- pair[bothRI | withT0, , drop = FALSE]
  if (!nrow(pair)) .ctmaTablesStop("No eligible CINT, MANIFESTMEANS or T0MEANS correlation pairs are available.")
  group <- ifelse(intercept[pair[,1L]] & intercept[pair[,2L]], "RICOR", "RIT0COR")
  left <- labels[pair[,1L]]
  right <- labels[pair[,2L]]
  rows <- data.frame(matrix = group, row = pair[,1L], col = pair[,2L],
                     parameter = paste(left, right, sep = " / "))
  if (is.matrix(saved)) {
    if (!identical(dim(saved), rep(length(labels), 2L)))
      .ctmaTablesStop("Saved random-intercept correlations do not match the parameter labels.")
    fields <- c(mean = "popcor_mean", sd = "model_popcor_sd",
                "2.5%" = "model_popcor_025", "50%" = "model_popcor_50",
                "97.5%" = "model_popcor_975")
    for (name in names(fields)) {
      value <- ri[[fields[[name]]]]
      if (is.matrix(value) && identical(dim(value), dim(saved))) rows[[name]] <- value[pair]
    }
  } else {
    # A raw correlation is unchanged by positive affine transforms on both margins.
    transforms <- varying$transform[match(unique(c(left, right)), varying$param)]
    transform <- gsub("[[:space:]()]", "", as.character(transforms))
    affine <- !is.na(transform) & grepl("^([0-9.]+\\*)?param$", transform)
    scale <- suppressWarnings(as.numeric(sub("\\*param$", "", transform)))
    affine <- affine & (transform %in% "param" | (is.finite(scale) & scale > 0))
    if (!all(affine))
      .ctmaTablesStop("Reporting native random-intercept correlations requires verified positive linear transforms of the selected random effects.")
    cor <- engine$rawpopcorr
    idx <- match(paste0(left, "__", right), rownames(cor))
    swap <- is.na(idx)
    idx[swap] <- match(paste0(right[swap], "__", left[swap]), rownames(cor))
    if (is.null(dim(cor)) || anyNA(idx) || !"mean" %in% colnames(cor))
      .ctmaTablesStop("The saved ctsem correlation table does not contain the selected parameter pairs.")
    stats <- intersect(c("mean", "sd", "2.5%", "50%", "97.5%", "n_eff", "Rhat"), colnames(cor))
    for (name in stats) rows[[name]] <- cor[idx, name]
  }

  # Reuse the existing statistic names so estimates and correlations share columns.
  for (name in setdiff(names(rows), c("matrix", "row", "col", "parameter"))) {
    matchName <- match(tolower(name), tolower(names(d)))
    if (!is.na(matchName)) names(rows)[names(rows) == name] <- names(d)[matchName]
  }
  for (name in setdiff(names(rows), names(d))) d[[name]] <- NA_real_
  for (name in setdiff(names(d), names(rows))) rows[[name]] <- NA
  rbind(d, rows[, names(d), drop = FALSE])
}
