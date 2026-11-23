test_that("TI summary rows match exact predictor prefixes", {
  labels <- c("tip_TI10_V1toV1", "tip_TI1_V1toV1", "tip_TI2_TI1",
              "tip_TI1_V2toV2", "tip_TI100_V2toV2")
  expect_identical(.ctmaTipredSummaryRows(labels, "TI1"), c(2L, 4L))
  expect_identical(.ctmaTipredSummaryRows(labels, c("TI1", "TI10")), c(2L, 4L, 1L))
  expect_identical(.ctmaTipredSummaryRows(labels, "TI2"), 3L)
  expect_identical(.ctmaTipredSummaryRows(labels, "TI3"), integer())
})
