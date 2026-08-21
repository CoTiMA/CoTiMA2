.make_re_red_het_fit <- function(sd_draws, order = c("drift", "diff"),
    studies = c("A", "B", "C"), priors = list(sd_prior = c(.25, .25)),
    moderator = FALSE) {
  specs <- list(
    drift = data.frame(
      re = 1L, matrix = "DRIFT", row = 1L, col = 1L,
      parameter = "V1toV1", type = 1L, index = 1L,
      transform = "param", block = "dynamic", stringsAsFactors = FALSE
    ),
    diff = data.frame(
      re = 2L, matrix = "DIFFUSION", row = 1L, col = 1L,
      parameter = "diff_eta1", type = 1L, index = 2L,
      transform = "exp(param)", block = "dynamic", stringsAsFactors = FALSE
    )
  )
  specification <- do.call(rbind, specs[order])
  rownames(specification) <- NULL
  setup <- matrix(0L, nrow = 2, ncol = 1,
    dimnames = list(NULL, "climate"))
  if (moderator) setup[, 1] <- seq_len(2)
  structure(list(
    estimation = "bayes",
    backend = "julia",
    draws = list(studyre_sd = as.matrix(sd_draws)[, match(order,
      c("drift", "diff")), drop = FALSE]),
    prepared = list(
      specification = specification,
      study_labels = studies,
      covariance_blocks = data.frame(
        block = "dynamic", structure = "independent"),
      priors = priors,
      standata = list(TIPREDEFFECTsetup = setup),
      original = list(ctstanmodelbase = list(TIpredNames = "climate")),
      study_dummy_names = character()
    )
  ), class = c("CoTiMAFit", "CoTiMAStudyREFit", "ctmaStudyREFit"))
}

test_that("ctmaReRedHet computes exact parameter and overall reductions", {
  null <- .make_re_red_het_fit(cbind(rep(2, 100), rep(4, 100)))
  moderated <- .make_re_red_het_fit(
    cbind(rep(1, 100), rep(2, 100)), moderator = TRUE)

  result <- suppressMessages(ctmaReRedHet(
    null, moderated, comparison_draws = 500, seed = 12))

  expect_s3_class(result, "ctmaReRedHet")
  expect_equal(result$parameter_reduction$unconditional_variance, c(4, 16))
  expect_equal(result$parameter_reduction$moderated_variance, c(1, 4))
  expect_equal(result$parameter_reduction$proportional_reduction, c(.75, .75))
  expect_equal(result$parameter_reduction$lower, c(.75, .75))
  expect_equal(result$overall_reduction$proportional_reduction, .75)
  expect_equal(result$parameter_reduction$probability_above_threshold,
    c(1, 1))
  expect_identical(result$moderated_moderators, "climate")
  expect_length(result$unconditional_moderators, 0L)
  expect_identical(summary(result), result)
  expect_output(print(result), "Overall")
})

test_that("ctmaReRedHet semantically matches reordered random effects", {
  null <- .make_re_red_het_fit(cbind(rep(2, 10), rep(4, 10)))
  moderated <- .make_re_red_het_fit(
    cbind(rep(1, 10), rep(2, 10)), order = c("diff", "drift"))

  result <- suppressMessages(ctmaReRedHet(null, moderated,
    comparison_draws = 50))
  expect_equal(result$parameter_reduction$proportional_reduction, c(.75, .75))
})

test_that("ctmaReRedHet reports negative reductions", {
  null <- .make_re_red_het_fit(cbind(rep(1, 20), rep(2, 20)))
  moderated <- .make_re_red_het_fit(cbind(rep(2, 20), rep(3, 20)))
  result <- suppressMessages(ctmaReRedHet(null, moderated,
    comparison_draws = 100))
  expect_true(all(result$parameter_reduction$proportional_reduction < 0))
  expect_equal(result$parameter_reduction$probability_above_threshold,
    c(0, 0))
})

test_that("ctmaReRedHet enforces comparable models", {
  null <- .make_re_red_het_fit(cbind(2, 4))
  different_studies <- .make_re_red_het_fit(cbind(1, 2), studies = c("A", "B"))
  expect_error(
    suppressMessages(ctmaReRedHet(null, different_studies)),
    "same studies"
  )

  different_spec <- .make_re_red_het_fit(cbind(1, 2))
  different_spec$prepared$specification$parameter[1] <- "other"
  expect_error(
    suppressMessages(ctmaReRedHet(null, different_spec)),
    "same random-effect specification"
  )

  zero <- .make_re_red_het_fit(cbind(0, 2))
  expect_error(
    suppressMessages(ctmaReRedHet(zero, null)),
    "zero variance"
  )

  posterior <- .make_re_red_het_fit(cbind(rep(2, 10), rep(4, 10)))
  expect_error(
    suppressMessages(ctmaReRedHet(null, posterior)),
    "two posterior fits or two point-estimate fits"
  )
})

test_that("ctmaReRedHet warns about different priors", {
  null <- .make_re_red_het_fit(cbind(rep(2, 10), rep(4, 10)))
  moderated <- .make_re_red_het_fit(
    cbind(rep(1, 10), rep(2, 10)), priors = list(sd_prior = c(.5, .5)))
  expect_warning(
    suppressMessages(ctmaReRedHet(null, moderated, comparison_draws = 20)),
    "different random-effect priors"
  )
})
