test_that("ctmaRedHet compares residual rather than total moderator variation", {
  null <- list(
    effects = matrix(c(-2, -1, 1, 2), ncol = 1,
      dimnames = list(NULL, "V1toV1")),
    standard_errors = matrix(rep(0.2, 4), ncol = 1),
    parameter_names = "V1toV1"
  )
  # A moderator may explain a large shift in the fitted study means.  Only the
  # remaining WEC deviation is supplied here, so that explained shift must not
  # be counted as residual heterogeneity.
  moderated <- list(
    effects = matrix(c(-0.2, -0.1, 0.1, 0.2), ncol = 1,
      dimnames = list(NULL, "V1toV1")),
    standard_errors = matrix(rep(0.2, 4), ncol = 1),
    parameter_names = "V1toV1"
  )
  null$analysis <- .ctmaBiGMetaAnalysis(null$effects, null$standard_errors)
  moderated$analysis <- .ctmaBiGMetaAnalysis(
    moderated$effects, moderated$standard_errors)
  result <- .ctmaRedHetReduction(null, moderated)

  expect_gt(result$proportional_reduction, 0)
  expect_equal(result$descriptive_proportional_reduction, 0.99,
    tolerance = 1e-12)
  expect_match(result$scale, "residual heterogeneity")
})

test_that("ctmaRedHet retains increases and handles zero denominators", {
  make_input <- function(values, standard_error) {
    x <- list(
      effects = matrix(values, ncol = 1,
        dimnames = list(NULL, "V1toV1")),
      standard_errors = matrix(rep(standard_error, length(values)), ncol = 1),
      parameter_names = "V1toV1"
    )
    x$analysis <- .ctmaBiGMetaAnalysis(x$effects, x$standard_errors)
    x
  }
  smaller <- make_input(c(-0.2, -0.1, 0.1, 0.2), 0.01)
  larger <- make_input(c(-2, -1, 1, 2), 0.01)
  increased <- .ctmaRedHetReduction(smaller, larger)
  expect_lt(increased$proportional_reduction, 0)

  homogeneous <- make_input(rep(0, 4), 1)
  undefined <- .ctmaRedHetReduction(homogeneous, homogeneous)
  expect_true(is.na(undefined$proportional_reduction))
})

test_that("ctmaRedHet WEC extraction is backend neutral", {
  fit <- list(
    n.latent = 1L,
    studyList = vector("list", 3L),
    parameterNames = list(DRIFT = "V1toV1"),
    modelResults = list(
      DRIFTCoeffViaWECMean = as.list(c(-0.2, 0, 0.2)),
      DRIFTCoeffViaWECSD = as.list(rep(0.1, 3))
    )
  )
  out <- .ctmaRedHetWEC(fit, "fit")
  expect_equal(dim(out$effects), c(3L, 1L))
  expect_equal(out$effects[, 1], c(-0.2, 0, 0.2))
  expect_equal(out$standard_errors[, 1], rep(0.1, 3))
})

test_that("ctmaRedHet continuous output retains the legacy matrix orientation", {
  make_fit <- function(scale = 1, moderated = FALSE) list(
    n.latent = 1L,
    studyList = vector("list", 4L),
    parameterNames = list(DRIFT = "V1toV1"),
    modelResults = list(
      DRIFTCoeffViaWECMean = as.list(scale * c(-2, -1, 1, 2)),
      DRIFTCoeffViaWECSD = as.list(rep(0.2, 4))
    ),
    argumentList = list(
      WEC = TRUE,
      ind.mod.number = integer(),
      mod.number = if (moderated) 1L else integer(),
      scaleTime = 1
    )
  )
  result <- .ctmaRedHetContinuous(
    make_fit(), make_fit(0.1, moderated = TRUE), NULL, 4L, quote(test()))
  fixed <- result$summary$continuousTime[[
    "fixedEffects (w/o Mod | with Mod)"
  ]]
  expect_equal(dim(fixed), c(8L, 2L))
  expect_s3_class(result, "ctmaRedHet")
  expect_gt(result$heterogeneity_reduction$proportional_reduction, 0)
})

test_that("discrete residuals are computed after the matrix exponential", {
  fake_draws <- list(
    total = list(matrix(c(-0.5, -0.4), ncol = 1)),
    predicted = list(matrix(c(-0.2, -0.2), ncol = 1)),
    n_latent = 1L
  )
  testthat::local_mocked_bindings(
    .ctmaRedHetDriftDraws = function(...) fake_draws
  )
  fit <- list(n.latent = 1L, parameterNames = list(DRIFT = "V1toV1"))
  result <- .ctmaRedHetDiscreteEffects(fit, dt = 2)
  expected <- exp(2 * c(-0.5, -0.4)) - exp(2 * c(-0.2, -0.2))
  expect_equal(unname(result$effects[1, 1]), mean(expected))
  expect_equal(unname(result$standard_errors[1, 1]), stats::sd(expected))
  expect_false(isTRUE(all.equal(result$effects[1, 1],
    exp(2 * mean(c(-0.3, -0.2))))))
})
