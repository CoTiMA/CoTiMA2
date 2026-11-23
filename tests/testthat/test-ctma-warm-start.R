make_warm_start_fixture <- function(
    nsubjects = 1000L, nparams = 10L, pooled = TRUE) {
  rows <- rep(seq_len(nsubjects), each = 3L)
  study <- rep(rep(0:1, each = nsubjects / 2L), each = 3L)
  binary_person <- rep(c(rep(0, nsubjects - 30L), rep(1, 30L)), each = 3L)
  fit <- list(standata = list(
    nsubjects = nsubjects,
    nparams = nparams,
    subject = rows,
    tipredsdata = matrix(rep(study[seq(1L, length(study), by = 3L)]),
      ncol = 1L),
    manifesttype = c(1L, 0L),
    Y = cbind(binary_person, seq_along(rows) / length(rows))
  ))
  source <- list(
    studyList = list(list(), list()),
    argumentList = list(
      randomInterceptStructure = if (pooled) "pooled" else "none"
    )
  )
  list(fit = fit, source = source)
}

test_that("warm-start modes are validated", {
  expect_identical(CoTiMA2:::.ctmaWarmStartMode("auto"), "auto")
  expect_identical(CoTiMA2:::.ctmaWarmStartMode(TRUE), "always")
  expect_identical(CoTiMA2:::.ctmaWarmStartMode(FALSE), "never")
  expect_error(CoTiMA2:::.ctmaWarmStartMode("sometimes"), "warm_start")
})

test_that("warm-start planner enlarges and stratifies complete subjects", {
  fixture <- make_warm_start_fixture()
  plan <- CoTiMA2:::.ctmaWarmStartPlan(
    fixture$fit, fixture$source, fraction = 0.05,
    subjects_per_parameter = 20L, minimum_per_study = 25L,
    maximum_subjects = 500L, binary_minimum = 20L, seed = 42L
  )
  repeated <- CoTiMA2:::.ctmaWarmStartPlan(
    fixture$fit, fixture$source, fraction = 0.05,
    subjects_per_parameter = 20L, minimum_per_study = 25L,
    maximum_subjects = 500L, binary_minimum = 20L, seed = 42L
  )
  expect_true(plan$use)
  expect_gte(plan$selected, 200L)
  expect_lt(plan$selected, 1000L)
  expect_true(all(plan$study_counts >= 25L))
  expect_gte(plan$binary_requirements[["1:0"]][["selected"]], 20L)
  expect_gte(plan$binary_requirements[["1:1"]][["selected"]], 20L)
  expect_identical(plan$subjects, repeated$subjects)
  expect_true(all(tabulate(
    fixture$fit$standata$subject, nbins = 1000L
  )[plan$subjects] >= 3L))
})

test_that("automatic warm start skips pilots that approach full data", {
  fixture <- make_warm_start_fixture(nsubjects = 100L, nparams = 5L)
  plan <- CoTiMA2:::.ctmaWarmStartPlan(
    fixture$fit, fixture$source, mode = "auto",
    subjects_per_parameter = 20L, minimum_per_study = 25L
  )
  expect_false(plan$use)
  expect_match(plan$reason, "80%|full dataset")
})

test_that("Julia ML backend exposes guarded warm-start defaults", {
  defaults <- ctmaBackendDefaults("julia", "ML", "ctmaFit")
  expect_identical(defaults$warm_start, "auto")
  expect_identical(defaults$warm_start_fraction, 0.10)
  expect_identical(defaults$warm_start_subjects_per_parameter, 20L)
  expect_identical(defaults$warm_start_minimum_per_study, 25L)
  expect_identical(defaults$warm_start_maximum_subjects, 5000L)
  expect_identical(defaults$warm_start_binary_minimum, 20L)
  expect_identical(defaults$warm_start_maxiters, 250L)
})
