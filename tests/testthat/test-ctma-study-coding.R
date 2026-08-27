test_that("study coding defaults preserve study structure", {
  resolved <- CoTiMA2:::.ctma_resolve_study_coding(
    studyCoding = "dummy",
    studyVarying = c("CINT", "DIFFUSION", "T0VAR"),
    WEC = FALSE
  )
  expect_identical(resolved$coding, "dummy")
  expect_identical(resolved$varying, c("CINT", "DIFFUSION", "T0VAR"))
})

test_that("legacy WEC maps to parameter-selective weighted coding", {
  expect_warning(
    resolved <- CoTiMA2:::.ctma_resolve_study_coding(
      studyCoding = "dummy", studyVarying = NULL, WEC = TRUE
    ),
    "WEC=TRUE.*deprecated"
  )
  expect_identical(resolved$coding, "weighted")
  expect_identical(resolved$varying, c("CINT", "DIFFUSION", "T0VAR"))
  expect_false("DRIFT" %in% resolved$varying)
})

test_that("invalid study coding specifications fail early", {
  expect_error(
    CoTiMA2:::.ctma_resolve_study_coding("reference", NULL, FALSE),
    "studyCoding"
  )
  expect_error(
    CoTiMA2:::.ctma_resolve_study_coding("dummy", "UNKNOWN", FALSE),
    "Unsupported 'studyVarying'"
  )
})

test_that("parameter-selective study effects affect only requested families", {
  model <- list(pars = data.frame(
    matrix = c("DRIFT", "DIFFUSION", "T0VAR", "CINT", "T0MEANS"),
    param = c("V1toV1", "diff_eta1", "T0cov_eta1", "cintV1", "MeanV1"),
    TI1_effect = TRUE,
    TI2_effect = TRUE,
    stringsAsFactors = FALSE
  ))
  selected <- CoTiMA2:::.ctma_apply_study_varying(
    model,
    study_effect_columns = c("TI1_effect", "TI2_effect"),
    studyVarying = c("CINT", "DIFFUSION", "T0VAR")
  )
  active <- selected$pars$matrix[
    selected$pars$TI1_effect | selected$pars$TI2_effect
  ]
  expect_setequal(active, c("CINT", "DIFFUSION", "T0VAR"))
  expect_false(any(selected$pars[selected$pars$matrix == "DRIFT", 3:4]))
})

test_that("phantom CINT coding targets only phantom intercept means", {
  model <- list(pars = data.frame(
    matrix = c("T0MEANS", "T0MEANS", "T0VAR", "DRIFT"),
    param = c("MeanV1", "MeanV1_cint", "Cov_cint1", "V1toV1"),
    TI1_effect = TRUE,
    stringsAsFactors = FALSE
  ))
  selected <- CoTiMA2:::.ctma_apply_study_varying(
    model,
    study_effect_columns = "TI1_effect",
    studyVarying = "CINT",
    phantom_intercepts = TRUE
  )
  expect_false(selected$pars$TI1_effect[1])
  expect_true(selected$pars$TI1_effect[2])
  expect_false(selected$pars$TI1_effect[3])
  expect_false(selected$pars$TI1_effect[4])
})
