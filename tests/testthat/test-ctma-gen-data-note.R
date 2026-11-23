test_that("ctmaGenData extraction note uses the directory and result name", {
  expect_message(
    CoTiMA2:::.ctmaGenDataExtractionNote("/tmp/generated/", "studies"),
    paste0(
      "ctmaExtract\\(activeDirectory = \"/tmp/generated/\", ",
      "ctmaGenDataList = studies, envir = globalenv\\(\\)\\)"
    )
  )
})

test_that("ctmaGenData extraction note explains its fallback name", {
  expect_message(
    CoTiMA2:::.ctmaGenDataExtractionNote("/tmp/generated/"),
    "replace 'generatedData'.*specify resultName"
  )
})

test_that("ctmaGenData validates resultName before generating data", {
  expect_error(
    ctmaGenData(resultName = character()),
    "'resultName' must be NULL or one non-empty character string"
  )
})

test_that("ctmaGenData preserves the legacy positional arguments", {
  expect_identical(
    names(formals(ctmaGenData))[23:29],
    c("T0means", "T0var", "TIpreds", "tpoints", "tpointTargets", "useRawData", "resultName")
  )
})
