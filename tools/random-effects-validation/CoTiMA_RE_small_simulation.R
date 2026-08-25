# Small recovery study for CoTiMA study-level random effects
#
# Purpose:
#   Test whether ctmaRePrep() and ctmaReFit() propagate a study-level random
#   effect into DRIFT and recover clearly separated study parameters.
#
# Design choices:
#   * one latent/manifest variable;
#   * four studies with moderate continuous-time drift values;
#   * a common diffusion parameter (no study-specific drift--diffusion ridge);
#   * six equally spaced observations per person;
#   * no individual-level random parameters;
#   * independent, centered study effects.

set.seed(20260810)

library(ctsem)
library(CoTiMA2)

trueDrift <- c(-0.15, -0.30, -0.55, -0.85)
numberStudies <- length(trueDrift)
personsPerStudy <- 60L
numberWaves <- 6L
commonDiffusion <- 0.35
initialVariance <- 1.00

# Generate each study from its own fixed-parameter continuous-time model.
studyData <- lapply(seq_len(numberStudies), function(study) {
  generatingModel <- ctModel(
    type = "ct",
    n.latent = 1,
    n.manifest = 1,
    Tpoints = numberWaves,
    manifestNames = "V1",
    latentNames = "eta1",
    LAMBDA = matrix(1, 1, 1),
    T0MEANS = matrix(0, 1, 1),
    T0VAR = matrix(initialVariance, 1, 1),
    DRIFT = matrix(trueDrift[study], 1, 1),
    DIFFUSION = matrix(commonDiffusion, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1),
    MANIFESTVAR = matrix(0, 1, 1),
    CINT = matrix(0, 1, 1)
  )

  generated <- ctGenerate(
    generatingModel,
    n.subjects = personsPerStudy,
    Tpoints = numberWaves,
    dtmean = 1,
    logdtsd = 0,
    wide = FALSE
  )

  generated <- as.data.frame(generated)
  generated$id <- generated$id + (study - 1L) * personsPerStudy
  generated$study <- study

  # k - 1 ordinary dummy variables. Study 4 is the reference study.
  for (j in seq_len(numberStudies - 1L))
    generated[[paste0("TI", j)]] <- as.numeric(study == j)

  generated
})

datalong <- do.call(rbind, studyData)
row.names(datalong) <- NULL

# Create the freely estimated ctsem model used as the starting FI model.
analysisModel <- ctModel(
  type = "ct",
  n.latent = 1,
  n.manifest = 1,
  Tpoints = numberWaves,
  manifestNames = "V1",
  latentNames = "eta1",
  LAMBDA = matrix(1, 1, 1),
  T0MEANS = matrix(0, 1, 1),
  T0VAR = matrix("T0var_eta1", 1, 1),
  DRIFT = matrix("V1toV1", 1, 1),
  DIFFUSION = matrix("diff_eta1", 1, 1),
  MANIFESTMEANS = matrix(0, 1, 1),
  MANIFESTVAR = matrix(0, 1, 1),
  CINT = matrix(0, 1, 1),
  n.TIpred = numberStudies - 1L,
  TIpredNames = paste0("TI", seq_len(numberStudies - 1L))
)

# No legacy study-dummy effects are needed in this clean hierarchical test.
effectColumns <- paste0("TI", seq_len(numberStudies - 1L), "_effect")
analysisModel$pars[analysisModel$pars$matrix == "DRIFT", effectColumns] <- FALSE
analysisModel$pars[analysisModel$pars$matrix == "DIFFUSION", effectColumns] <- FALSE
analysisModel$pars[analysisModel$pars$matrix == "T0VAR", effectColumns] <- FALSE

# Build, but do not estimate, ctsem's Stan representation.
baseFit <- ctFit(
  datalong = datalong,
  model = analysisModel,
  fit = FALSE,
  optimize = FALSE,
  priors = TRUE,
  verbose = 0
)

# Add the argument structure normally supplied by CoTiMA.
baseFit$argumentList <- list(
  randomEffect = list(
    reDRIFT = matrix(TRUE, 1, 1),
    reDIFFUSION = matrix(FALSE, 1, 1),
    reT0var = matrix(FALSE, 1, 1)
  )
)

preparedSimulation <- ctmaRePrep(
  model = baseFit,
  study_id = rep(seq_len(numberStudies), each = personsPerStudy),
  reDRIFT = matrix(TRUE, 1, 1),
  reDIFFUSION = matrix(FALSE, 1, 1),
  reT0VAR = matrix(FALSE, 1, 1),
  covariance = "independent",
  parameterization = "centered"
)

# Plausible raw initialization for the ctsem DRIFT transformation.
inverseDriftTransform <- function(drift) {
  target <- (-drift - 1e-6) / 2
  -0.5 * log(expm1(target))
}

populationDrift <- mean(trueDrift)
rawDriftInitial <- inverseDriftTransform(populationDrift)

initialValues <- function(chain_id = 1L) {
  list(
    # There are three free population parameters in this model: DRIFT,
    # DIFFUSION, and T0VAR. The latter two values are moderate raw-scale
    # starting points; exact initial values are not treated as known.
    rawpopmeans = array(
      c(rawDriftInitial, -1.5, -0.75),
      dim = 3L
    ),
    studyre_sd = array(0.25, dim = 1L),
    studyre_effect = matrix(0, numberStudies, 1L)
  )
}

simulationFit <- ctmaReFit(
  preparedSimulation,
  estimation = "bayes",
  chains = 2,
  iter = 200,
  cores = 2,
  refresh = 10,
  fit_args = list(
    warmup = 100,
    init = initialValues,
    control = list(
      adapt_delta = 0.90,
      max_treedepth = 10
    )
  )
)

simulationSummary <- summary(simulationFit)
estimatedDrift <- subset(simulationSummary$study_parameters, matrix == "DRIFT")

recovery <- data.frame(
  study = seq_len(numberStudies),
  true = trueDrift,
  posterior_mean = estimatedDrift$mean,
  posterior_median = estimatedDrift$median,
  lower = estimatedDrift$lower,
  upper = estimatedDrift$upper
)

samplerParameters <- rstan::get_sampler_params(
  simulationFit$engine_fit,
  inc_warmup = FALSE
)

diagnostics <- do.call(rbind, lapply(samplerParameters, function(x) {
  c(
    divergences = sum(x[, "divergent__"]),
    proportion_at_maximum_treedepth = mean(x[, "treedepth__"] >= 10),
    mean_acceptance = mean(x[, "accept_stat__"]),
    mean_leapfrog_steps = mean(x[, "n_leapfrog__"])
  )
}))

print(recovery)
print(diagnostics)
print(rstan::summary(
  simulationFit$engine_fit,
  pars = c("rawpopmeans", "studyre_sd", "studyre_effect")
)$summary)

saveRDS(
  list(
    design = list(
      trueDrift = trueDrift,
      commonDiffusion = commonDiffusion,
      personsPerStudy = personsPerStudy,
      numberWaves = numberWaves
    ),
    recovery = recovery,
    diagnostics = diagnostics,
    fit = simulationFit
  ),
  "CoTiMA_RE_small_simulation_results.rds"
)
