test_that("ctmaBiG common- and random-effects calculations match formulas", {
  effect <- matrix(c(-0.10, -0.22, -0.18, -0.35, -0.12, -0.27),
    ncol = 1, dimnames = list(NULL, "V1toV1"))
  standard_error <- matrix(c(0.04, 0.08, 0.05, 0.10, 0.06, 0.07),
    ncol = 1)

  result <- CoTiMA2:::.ctmaBiGMetaAnalysis(effect, standard_error)
  weight <- 1 / standard_error[, 1]^2
  common_mean <- sum(weight * effect[, 1]) / sum(weight)
  common_variance <- 1 / sum(weight)
  q <- sum(weight * (effect[, 1] - common_mean)^2)
  df <- nrow(effect) - 1
  c_value <- sum(weight) - sum(weight^2) / sum(weight)
  tau2 <- max(0, (q - df) / c_value)
  random_weight <- 1 / (standard_error[, 1]^2 + tau2)
  random_mean <- sum(random_weight * effect[, 1]) / sum(random_weight)
  random_variance <- 1 / sum(random_weight)
  prediction_half_width <- stats::qt(0.975, df = nrow(effect) - 2) *
    sqrt(tau2 + random_variance)

  expect_equal(result$common["FixedEffect_Drift", 1], common_mean)
  expect_equal(result$common["FixedEffect_DriftVariance", 1],
    common_variance)
  expect_equal(result$common["FixedEffect_DriftProb", 1],
    2 * stats::pnorm(-abs(common_mean / sqrt(common_variance))))
  expect_equal(result$heterogeneity["Q_Drift", 1], q)
  expect_equal(result$heterogeneity["Heterogeneity_df", 1], df)
  expect_equal(result$heterogeneity["tau2Drift", 1], tau2)
  expect_equal(result$random["RandomEffecttot_Drift", 1], random_mean)
  expect_equal(result$random["RandomEffecttot_DriftVariance", 1],
    random_variance)
  expect_equal(result$random["RandomEffecttot_DriftLowerLimitPI", 1],
    random_mean - prediction_half_width)
  expect_equal(result$random["RandomEffecttot_DriftUpperLimitPI", 1],
    random_mean + prediction_half_width)
})

test_that("ctmaBiG truncates impossible negative heterogeneity", {
  effect <- matrix(c(0.100, 0.101, 0.099, 0.102), ncol = 1)
  standard_error <- matrix(rep(0.1, 4), ncol = 1)
  result <- CoTiMA2:::.ctmaBiGMetaAnalysis(effect, standard_error)

  expect_equal(result$heterogeneity["tau2Drift", 1], 0)
  expect_equal(result$heterogeneity["I2_Drift", 1], 0)
  expect_equal(result$heterogeneity["H2_Drift", 1], 1)
  expect_equal(result$random["RandomEffecttot_Drift", 1],
    result$common["FixedEffect_Drift", 1])
})

test_that("ctmaBiG uses parameter-specific study counts", {
  effect <- cbind(
    complete = c(-0.1, -0.2, -0.3, -0.4),
    incomplete = c(-0.1, NA, -0.3, -0.4),
    insufficient = c(NA, NA, -0.3, NA)
  )
  standard_error <- matrix(0.1, nrow = 4, ncol = 3,
    dimnames = list(NULL, colnames(effect)))
  result <- CoTiMA2:::.ctmaBiGMetaAnalysis(effect, standard_error)

  expect_equal(unname(result$k), c(4L, 3L, 1L))
  expect_equal(unname(result$heterogeneity["Heterogeneity_df", 1:2]),
    c(3, 2))
  expect_equal(unname(result$heterogeneity["AnalysedStudies", ]),
    c(4, 3, 1))
  expect_true(all(is.na(result$common[, "insufficient"])))
})

test_that("ctmaBiG meta-analysis validates its inputs", {
  expect_error(
    CoTiMA2:::.ctmaBiGMetaAnalysis(matrix(1, 2, 1), matrix(1, 1, 2)),
    "identical dimensions"
  )
  expect_error(
    CoTiMA2:::.ctmaBiGMetaAnalysis(matrix(1, 2, 1), matrix(1, 2, 1),
      confidence = 1),
    "between zero and one"
  )
})

test_that("ctmaBiG extracts identical transformed values from Stan and Julia", {
  make_draws <- function(offset) {
    array(seq_len(20) / 10 + offset, dim = c(5, 2, 2))
  }
  transformedpars <- list(
    pop_DRIFT = make_draws(0),
    pop_DIFFUSIONcov = make_draws(10),
    pop_T0cov = make_draws(20)
  )
  stan_fit <- structure(
    list(stanfit = list(transformedpars = transformedpars)),
    class = "ctStanFit"
  )
  julia_fit <- structure(list(ctsem_fit = stan_fit), class = "ctmaJuliaFit")
  parameter_names <- paste0("p", seq_len(12))

  stan <- CoTiMA2:::.ctmaBiGExtractStudy(
    stan_fit, n_latent = 2, parameter_names = parameter_names)
  julia <- CoTiMA2:::.ctmaBiGExtractStudy(
    julia_fit, n_latent = 2, parameter_names = parameter_names)

  expect_identical(stan$backend, "stan")
  expect_identical(julia$backend, "julia")
  expect_equal(julia$estimate, stan$estimate)
  expect_equal(julia$standard_error, stan$standard_error)
  expect_identical(names(julia$estimate), parameter_names)
})

test_that("ctmaBiG study extraction handles mixed backends and exclusions", {
  draws <- array(seq_len(20) / 10, dim = c(5, 2, 2))
  transformedpars <- list(
    pop_DRIFT = draws,
    pop_DIFFUSIONcov = draws + 10,
    pop_T0VAR = draws + 20
  )
  stan_fit <- structure(
    list(stanfit = list(transformedpars = transformedpars)),
    class = "ctStanFit"
  )
  julia_fit <- structure(list(ctsem_fit = stan_fit), class = "ctmaJuliaFit")
  parameter_names <- paste0("p", seq_len(4))
  init_fit <- list(
    n.latent = 2L,
    studyFitList = list(stan_fit, julia_fit, stan_fit),
    modelResults = list(
      DRIFT = list(stats::setNames(rep(0, 4), parameter_names)),
      DIFFUSION = list(stats::setNames(rep(0, 4), parameter_names)),
      T0VAR = list(stats::setNames(rep(0, 4), parameter_names))
    )
  )

  result <- CoTiMA2:::.ctmaBiGExtractStudies(init_fit, exclude_studies = 2)

  expect_identical(result$included_studies, c(1L, 3L))
  expect_identical(result$backends, c("stan", NA_character_, "stan"))
  expect_true(all(is.na(result$estimates[2, ])))
  expect_equal(result$estimates[1, ], result$estimates[3, ])
  expect_error(
    CoTiMA2:::.ctmaBiGExtractStudies(init_fit, exclude_studies = 4),
    "valid study indices"
  )
})

test_that("ctmaBiG explains missing Julia transformed results", {
  julia_fit <- structure(list(ctsem_fit = NULL), class = "ctmaJuliaFit")
  expect_error(
    CoTiMA2:::.ctmaBiGExtractStudy(julia_fit, 1),
    "uncertainty=TRUE.*transformed=TRUE"
  )
})
