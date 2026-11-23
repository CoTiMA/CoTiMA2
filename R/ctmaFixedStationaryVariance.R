.ctmaFixedStationaryVarianceModel <- function(model, variance) {
  fail <- function(...) stop("stationaryVariance: ", ..., call. = FALSE)
  if (!is.numeric(variance) || length(variance) != 1L ||
      !is.finite(variance) || variance <= 0) fail("supply a finite positive variance.")
  p <- model$pars
  if (!is.data.frame(p) || model$n.latent != 1L ||
      any(model$manifesttype != 0)) fail("currently requires one Gaussian latent process.")
  if (any(p$matrix == "PARS")) fail("existing PARS models require an explicit mapping.")
  varying <- grep("^indvarying($|_)", names(p), value = TRUE)
  if (any(vapply(varying, function(v) any(p[[v]] %in% TRUE), logical(1))))
    fail("person-level/random matrix effects are not supported; select study DRIFT effects through backendArgs after model construction.")
  for (m in c("T0MEANS", "CINT", "MANIFESTMEANS")) {
    v <- p$value[p$matrix == m]
    if (!length(v) || anyNA(v) || any(v != 0)) fail("requires fixed zero ", m, ".")
  }
  dr <- which(p$matrix == "DRIFT")
  di <- which(p$matrix == "DIFFUSION")
  t0 <- which(p$matrix == "T0VAR")
  if (length(dr) != 1L || length(di) != 1L || length(t0) != 1L)
    fail("expected scalar DRIFT, DIFFUSION and T0VAR.")
  label <- p$param[dr]
  if (is.na(label) || !grepl("^[[:alpha:]_][[:alnum:]_]*$", label) ||
      any(p$param[-dr] == label, na.rm = TRUE)) fail("DRIFT must have its own free parameter label.")
  transform <- gsub("[[:space:]]", "", p$transform[dr])
  transform <- gsub("0.000001", "1e-06", transform, fixed = TRUE)
  if (is.na(transform) || !transform %in% c("-(1e-06+2*log1p_exp(-(2*param)))",
      "-(1e-6+2*log1p_exp(-(2*param)))", "-exp(param)"))
    fail("DRIFT must use a supported strictly negative transform.")
  if (!is.na(p$value[di])) fail("remove fixed DIFFUSION; it is derived from DRIFT and stationaryVariance.")
  if (!is.na(p$value[t0]) && abs(p$value[t0]^2 - variance) > 1e-10 * max(1, variance))
    fail("fixed T0VAR is inconsistent with stationaryVariance (T0VAR is a square-root parameter).")
  ti <- grep("_effect$", names(p), value = TRUE)
  if (any(vapply(ti, function(v) any(p[[v]][c(di, t0)] %in% TRUE), logical(1))))
    fail("DIFFUSION and T0VAR cannot have separate study effects; use studyVarying='DRIFT' or character().")
  original <- p
  aux <- p[dr, , drop = FALSE]
  aux$matrix <- "PARS"
  aux$row <- aux$col <- 1L
  write_cell <- function(i, expression = NA_character_, value = NA_real_) {
    p$param[i] <<- expression
    p$value[i] <<- value
    p$transform[i] <<- NA_character_
    for (v in c(varying, ti)) p[[v]][i] <<- FALSE
  }
  write_cell(dr, paste0("1*", label))
  write_cell(di, paste0("sqrt(-2*", label, "*", format(variance, digits = 17L), ")"))
  write_cell(t0, value = sqrt(variance))
  model$pars <- rbind(p, aux)
  model$determinedDiffusions <- list(
    parameterization = "fixed univariate stationary process variance",
    stationaryVariance = variance, drift_parameters = label,
    original = original,
    study_effect_columns = ti[vapply(ti, function(v) isTRUE(aux[[v]]), logical(1))],
    study_covariances = "fixed stationary variance; diffusion follows each study drift",
    auxiliary = label, suggested_raw = stats::setNames(0, label),
    note = "The DRIFT coordinate, transformation and study contrasts are preserved. No free variance coordinate is added.")
  model
}
