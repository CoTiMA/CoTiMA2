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

test_that("Julia indVarying uses the pooled phantom representation", {
  translated <- CoTiMA2:::.ctma_julia_pooled_intercepts(
    "CINT", FALSE, "julia", announce = FALSE
  )
  expect_identical(translated$indVarying, FALSE)
  expect_identical(translated$randomIntercepts, "CINT")
  expect_true(translated$pooled)

  stan <- CoTiMA2:::.ctma_julia_pooled_intercepts(
    "CINT", FALSE, "stan", announce = FALSE
  )
  expect_identical(stan$indVarying, "CINT")
  expect_identical(stan$randomIntercepts, FALSE)
  expect_false(stan$pooled)

  unrestricted <- CoTiMA2:::.ctma_julia_pooled_intercepts(
    FALSE, "CINT", "julia", announce = FALSE
  )
  expect_false(unrestricted$pooled)
})

test_that("pooling fixes only TI effects on phantom covariances", {
  t0var <- matrix(0, 4, 4)
  t0var[lower.tri(t0var, diag = TRUE)] <- paste0(
    "v", seq_len(sum(lower.tri(t0var, diag = TRUE)))
  )
  model <- ctsem::ctModel(
    type = "ct", n.latent = 4, n.manifest = 2, n.TIpred = 2,
    latentNames = c("eta1", "eta2", "eta1_cint", "eta2_cint"),
    manifestNames = c("y1", "y2"),
    LAMBDA = cbind(diag(2), matrix(0, 2, 2)),
    T0MEANS = matrix(paste0("m", 1:4), 4, 1),
    T0VAR = t0var
  )
  effect_columns <- grep("_effect$", names(model$pars), value = TRUE)
  model$pars[, effect_columns] <- TRUE
  original_parameters <- model$pars$param
  pooled <- CoTiMA2:::.ctma_pool_phantom_covariance(model, 2L)
  phantom_covariance <- pooled$pars$matrix == "T0VAR" &
    (pooled$pars$row > 2L | pooled$pars$col > 2L)
  ordinary_covariance <- pooled$pars$matrix == "T0VAR" &
    pooled$pars$row <= 2L & pooled$pars$col <= 2L
  phantom_means <- pooled$pars$matrix == "T0MEANS" & pooled$pars$row > 2L

  expect_true(all(!unlist(pooled$pars[phantom_covariance, effect_columns])))
  expect_true(all(unlist(pooled$pars[ordinary_covariance, effect_columns])))
  expect_true(all(unlist(pooled$pars[phantom_means, effect_columns])))
  expect_identical(pooled$pars$param, original_parameters)
  expect_true(any(!is.na(pooled$pars$param[phantom_covariance])))
})
