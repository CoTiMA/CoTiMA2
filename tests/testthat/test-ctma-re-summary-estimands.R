test_that("hierarchical summaries distinguish raw centres and study averages", {
  spec <- data.frame(
    type = 1L, index = 1L, matrix = "DIFFUSION", row = 1L, col = 1L,
    parameter = "diffusion", block = "DIFFUSION", block_id = 1L,
    transform = "exp(param)", t0_row = 0L, t0_col = 0L,
    stringsAsFactors = FALSE
  )
  effects <- array(c(0, 0, log(3), log(2)), dim = c(2, 2, 1))
  fit <- structure(list(
    optimizer_uncertainty = "posterior",
    prepared = list(
      specification = spec,
      study_labels = c("study1", "study2"),
      covariance_blocks = data.frame(structure = "independent")
    ),
    draws = list(
      studyre_effect = effects,
      studyre_sd = matrix(c(0.1, 0.2), ncol = 1),
      rawpopmeans = matrix(c(0, log(2)), ncol = 1),
      studyre_T0VAR = NULL,
      studyre_T0cov = NULL
    )
  ), class = c("CoTiMAFit", "CoTiMAStudyREFit", "ctmaStudyREFit"))

  expect_output(
    result <- summary(fit),
    "Hierarchical summary estimands"
  )
  expect_equal(result$population_parameters_raw_centre$mean, 1.5)
  expect_equal(result$population_parameters_study_average$mean, 2.5)
  expect_identical(
    result$population_parameters,
    result$population_parameters_raw_centre
  )
  expect_match(result$estimand_interpretation, "meta-analytic population")
})

test_that("T0 covariance summaries average jointly transformed study matrices", {
  spec <- data.frame(
    type = 2L, index = 1L, matrix = "T0VAR", row = 1L, col = 1L,
    parameter = "T0 variance", block = "T0VAR", block_id = 1L,
    transform = "param", t0_row = 1L, t0_col = 1L,
    stringsAsFactors = FALSE
  )
  t0var <- array(0, dim = c(2, 2, 1, 1))
  t0var[, 1, 1, 1] <- c(1, 3)
  t0var[, 2, 1, 1] <- c(5, 7)
  t0cov <- t0var * 2
  fit <- structure(list(
    optimizer_uncertainty = "posterior",
    prepared = list(
      specification = spec,
      study_labels = c("study1", "study2"),
      covariance_blocks = data.frame(structure = "independent")
    ),
    draws = list(
      studyre_effect = array(0, dim = c(2, 2, 1)),
      studyre_sd = matrix(c(0.1, 0.2), ncol = 1),
      rawpopmeans = matrix(c(0, 0), ncol = 1),
      studyre_T0VAR = t0var,
      studyre_T0cov = t0cov
    )
  ), class = c("CoTiMAFit", "CoTiMAStudyREFit", "ctmaStudyREFit"))

  invisible(capture.output(result <- summary(fit)))
  expect_equal(result$population_T0VAR_study_average$mean, 4)
  expect_equal(result$population_T0cov_study_average$mean, 8)
})
