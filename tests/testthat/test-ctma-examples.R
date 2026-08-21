test_that("example objects require explicit loading and unload safely", {
  directory <- withr::local_tempdir()
  exampleA <- list(value = 1)
  save(exampleA, file = file.path(directory, "exampleA.rda"))
  testthat::local_mocked_bindings(
    .ctmaExampleDirectory = function() directory,
    .ctmaVersionNote = function() invisible(NULL)
  )
  target <- new.env(parent = emptyenv())
  expect_identical(ctmaLoadExamples("exampleA", target), "exampleA")
  expect_true(exists("exampleA", target, inherits = FALSE))
  removed <- ctmaUnloadExamples(envir = target)
  expect_identical(removed$removed, "exampleA")
  expect_false(exists("exampleA", target, inherits = FALSE))

  ctmaLoadExamples("exampleA", target)
  target$exampleA$value <- 2
  kept <- ctmaUnloadExamples(envir = target)
  expect_identical(kept$modified, "exampleA")
  expect_identical(target$exampleA$value, 2)
})

test_that("example loading never overwrites silently", {
  directory <- withr::local_tempdir()
  exampleA <- 1
  save(exampleA, file = file.path(directory, "exampleA.rda"))
  testthat::local_mocked_bindings(
    .ctmaExampleDirectory = function() directory,
    .ctmaVersionNote = function() invisible(NULL)
  )
  target <- new.env(parent = emptyenv())
  target$exampleA <- 99
  expect_error(
    ctmaLoadExamples("exampleA", target),
    "Remove or rename.*overwrite = TRUE.*ctmaSecuredObjects"
  )
  expect_identical(target$exampleA, 99)
  expect_message(
    ctmaLoadExamples("exampleA", target, overwrite = TRUE),
    "Secured 1 existing object"
  )
  expect_identical(target$exampleA, 1)
  expect_identical(target$ctmaSecuredObjects$exampleA, 99)
  ctmaUnloadExamples("exampleA", target)
  expect_identical(target$exampleA, 99)
  expect_identical(target$ctmaSecuredObjects$exampleA, 99)
})

test_that("repeated overwrites preserve every secured value", {
  directory <- withr::local_tempdir()
  exampleA <- 1
  save(exampleA, file = file.path(directory, "exampleA.rda"))
  testthat::local_mocked_bindings(
    .ctmaExampleDirectory = function() directory,
    .ctmaVersionNote = function() invisible(NULL)
  )
  target <- new.env(parent = emptyenv())
  target$exampleA <- 98
  suppressMessages(ctmaLoadExamples("exampleA", target, overwrite = TRUE))
  target$exampleA <- 99
  suppressMessages(ctmaLoadExamples("exampleA", target, overwrite = TRUE))
  expect_identical(target$ctmaSecuredObjects$exampleA, 98)
  expect_identical(target$ctmaSecuredObjects$exampleA.1, 99)
})

test_that("all packaged complex examples unload except a changed object", {
  target <- new.env(parent = emptyenv())
  loaded <- suppressMessages(ctmaLoadExamples(envir = target))
  expect_gt(length(loaded), 100L)
  target$targetVariables1 <- c(target$targetVariables1, "changed")
  result <- suppressMessages(ctmaUnloadExamples(envir = target))
  expect_setequal(result$modified, "targetVariables1")
  expect_setequal(result$removed, setdiff(loaded, "targetVariables1"))
  expect_identical(ls(target, all.names = TRUE), "targetVariables1")
})
