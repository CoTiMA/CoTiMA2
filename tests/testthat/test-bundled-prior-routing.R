test_that("bundled MAP requests fail before setup while ML and NUTS remain allowed", {
  expect_error(.ctmaCheckBundledPriorRequest("JULIA", TRUE, TRUE), "does not support MAP")
  expect_error(ctmaFit(backend="julia", optimize=TRUE, priors=TRUE), "does not support MAP")
  expect_error(ctmaInit(backend="julia", optimize=TRUE, priors=TRUE), "does not support MAP")
  expect_no_error(.ctmaCheckBundledPriorRequest("julia", TRUE, FALSE))
  expect_no_error(.ctmaCheckBundledPriorRequest("julia", FALSE, TRUE))
  expect_no_error(.ctmaCheckBundledPriorRequest("ctsemjulia", TRUE, TRUE))
  expect_no_error(.ctmaCheckBundledPriorRequest("stan", TRUE, TRUE))
})
