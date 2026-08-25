test_that("location specifications are normalized and dimension checked", {
  normalized <- CoTiMA2:::.ctma_normalize_location_specs(NULL, NULL, 2, 3)
  expect_identical(normalized$cint, 0)
  expect_identical(normalized$manifestMeans, 0)
  expect_identical(
    CoTiMA2:::.ctma_normalize_location_specs("auto", "auto", 2, 3),
    list(cint = "auto", manifestMeans = "auto")
  )
  expect_error(
    CoTiMA2:::.ctma_normalize_location_specs(c(0, 0, 0), 0, 2, 3),
    "number of latent"
  )
  expect_error(
    CoTiMA2:::.ctma_normalize_location_specs(0, c(0, 0), 2, 3),
    "number of manifest"
  )
})

test_that("binary intercept validation uses final free parameters", {
  make_model <- function(cint = 0, manifest_mean = 0) {
    ctsem::ctModel(
      type = "ct", n.latent = 1, n.manifest = 1,
      latentNames = "eta", manifestNames = "y",
      LAMBDA = matrix(1, 1, 1), CINT = matrix(cint, 1, 1),
      MANIFESTMEANS = matrix(manifest_mean, 1, 1)
    )
  }
  expect_error(
    CoTiMA2:::.ctma_check_binary_intercepts(make_model(), 1),
    "neither a freely estimated CINT"
  )
  expect_no_error(CoTiMA2:::.ctma_check_binary_intercepts(
    make_model(cint = "cint"), 1
  ))
  expect_message(CoTiMA2:::.ctma_check_binary_intercepts(
    make_model(manifest_mean = "mean"), 1
  ), "MANIFESTMEANS")
  expect_warning(CoTiMA2:::.ctma_check_binary_intercepts(
    make_model(cint = "cint", manifest_mean = "mean"), 1
  ), "both freely estimated")
  expect_error(CoTiMA2:::.ctma_check_binary_intercepts(make_model(), c(1, 0)),
    "exactly one entry")
})

test_that("a free binary CINT does not require indVarying", {
  model <- ctsem::ctModel(
    type = "ct", n.latent = 1, n.manifest = 1,
    latentNames = "eta", manifestNames = "y",
    LAMBDA = matrix(1, 1, 1), CINT = matrix("cint", 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1)
  )
  model$pars$indvarying <- FALSE
  expect_no_error(CoTiMA2:::.ctma_check_binary_intercepts(model, 1))
})

test_that("random-intercept and indVarying location models are exempt", {
  model <- ctsem::ctModel(
    type = "ct", n.latent = 1, n.manifest = 1,
    latentNames = "eta", manifestNames = "y",
    LAMBDA = matrix(1, 1, 1), CINT = matrix(0, 1, 1),
    MANIFESTMEANS = matrix(0, 1, 1)
  )
  expect_no_error(CoTiMA2:::.ctma_check_binary_intercepts(
    model, 1, process_intercepts = TRUE
  ))
})
