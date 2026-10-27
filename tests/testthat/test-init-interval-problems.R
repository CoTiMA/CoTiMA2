test_that("interval problems can start with diffusion or T0VAR", {
  result <- .ctmaInitRecordIntervalProblem(
    list(), 1, "DIFFUSION", 2L, c("variance_X", "covariance_XY")
  )
  expect_length(result, 1L)
  expect_match(result[[1]], "Study 1.*DIFFUSION.*covariance_XY")
  result <- .ctmaInitRecordIntervalProblem(
    result, 2, "T0VAR", 1L, "initial_variance"
  )
  expect_length(result, 2L)
  expect_match(result[[2]], "Study 2.*T0VAR.*initial_variance")
  expect_match(result[[1]], "Study 1")
  expect_identical(
    .ctmaInitRecordIntervalProblem(result, 3, "DRIFT", integer(), NULL),
    result
  )
  expect_match(.ctmaInitRecordIntervalProblem(
    list(), 1, "T0VAR", 2L, character()
  )[[1]], "T0VAR\\[2\\]")
})
