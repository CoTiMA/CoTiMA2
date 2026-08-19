test_that("CoTiMA-RE covariance aliases are normalized", {
  expect_identical(CoTiMA:::.ctma_covariance_name("ure"), "independent")
  expect_identical(CoTiMA:::.ctma_covariance_name("ere"), "equicorrelated")
  expect_identical(CoTiMA:::.ctma_covariance_name("equal"), "equicorrelated")
  expect_identical(CoTiMA:::.ctma_covariance_name("cre"), "full")
  expect_error(CoTiMA:::.ctma_covariance_name("unknown"),
    "Unknown covariance structure")
})

test_that("unfitted CoTiMA wrappers are distinguished from ctsem setups", {
  wrapper <- structure(list(
    data = data.frame(id = 1),
    ctModel = list(pars = data.frame())
  ), class = "CoTiMAFit")
  expect_true(CoTiMA:::.ctma_is_unfitted_cotima(wrapper))

  ctsem_setup <- list(
    data = wrapper$data,
    ctModel = wrapper$ctModel,
    standata = list(),
    stanmodeltext = "parameters {} model {}"
  )
  expect_false(CoTiMA:::.ctma_is_unfitted_cotima(ctsem_setup))
  expect_identical(
    CoTiMA:::.ctma_generate_ctsem_setup(ctsem_setup, verbose = FALSE),
    ctsem_setup
  )
})

test_that("ctmaReFit exposes the optional Julia backend controls", {
  arguments <- names(formals(ctmaReFit))
  expect_true(all(c("backend", "backendArgs") %in% arguments))
  expect_identical(eval(formals(ctmaReFit)$backend), c("stan", "julia"))
  expect_true(file.exists(system.file(
    "julia", "CtsemJulia.jl", package = "CoTiMA"
  )))
})

test_that("automatic matrix blocks count correlation parameters correctly", {
  spec <- data.frame(
    re = 1:8,
    matrix = c(rep("DRIFT", 4), rep("DIFFUSION", 3), "T0VAR"),
    row = c(1, 2, 1, 2, 1, 2, 2, 2),
    col = c(1, 1, 2, 2, 1, 1, 2, 1),
    parameter = paste0("p", 1:8),
    stringsAsFactors = FALSE
  )

  result <- CoTiMA:::.ctma_re_blocks(
    spec = spec,
    covariance = "block",
    covarianceBlocks = list(
      DRIFT = "equicorrelated",
      DIFFUSION = "independent",
      T0VAR = "full"
    ),
    nlatent = 2,
    verbose = FALSE
  )

  expect_identical(result$blocks$block,
    c("DRIFT", "DIFFUSION", "T0VAR"))
  expect_identical(result$blocks$size, c(4L, 3L, 1L))
  expect_identical(result$blocks$structure,
    c("equicorrelated", "independent", "independent"))
  code <- CoTiMA:::.ctma_stan_block_code(result$blocks)
  expect_true(code$has_corr)
  expect_identical(code$ncorr, 1L)
  expect_match(code$par, "studyre_rho_b1", fixed = TRUE)
})

test_that("custom blocks require exactly the selected cells", {
  spec <- data.frame(
    re = 1L, matrix = "DRIFT", row = 1L, col = 1L,
    parameter = "drift_1_1", stringsAsFactors = FALSE
  )
  assignment <- matrix(NA_character_, 2, 2)
  assignment[1, 1] <- "dynamic"
  assignment[2, 2] <- "not_selected"

  expect_error(
    CoTiMA:::.ctma_re_blocks(
      spec = spec,
      covariance = "block",
      reBlocks = list(DRIFT = assignment),
      blockStructures = c(dynamic = "independent"),
      nlatent = 2,
      verbose = FALSE
    ),
    "not selected"
  )
})

test_that("MAP and conditional ML set opposite prior switches", {
  standata <- list(priors = 1L, studyre_use_hyperpriors = 1L)
  map_data <- CoTiMA:::.ctma_optimizer_data(standata, "map")
  ml_data <- CoTiMA:::.ctma_optimizer_data(standata, "ml")

  expect_identical(map_data$priors, 1L)
  expect_identical(map_data$studyre_use_hyperpriors, 1L)
  expect_identical(ml_data$priors, 0L)
  expect_identical(ml_data$studyre_use_hyperpriors, 0L)
  expect_error(CoTiMA:::.ctma_optimizer_data(list(), "map"),
    "expected 'priors' switch", fixed = TRUE)
})

test_that("optimizer output receives a leading draw dimension", {
  result <- CoTiMA:::.ctma_add_draw_dimension(list(
    scalar = 2,
    vector = c(1, 2),
    matrix = matrix(1:4, 2, 2)
  ))

  expect_identical(dim(result$scalar), c(1L, 1L))
  expect_identical(dim(result$vector), c(1L, 2L))
  expect_identical(dim(result$matrix), c(1L, 2L, 2L))
})

test_that("transformed random-effect SD respects the parameter transform", {
  alpha <- c(-0.5, 0.2)
  tau <- c(0.3, 0.7)

  expect_equal(
    CoTiMA:::.ctma_transformed_re_sd(alpha, tau, "param"),
    tau,
    tolerance = 1e-10
  )

  # If X is normal, exp(X) is lognormal with this analytic SD.
  expected <- sqrt(
    (exp(tau^2) - 1) * exp(2 * alpha + tau^2)
  )
  expect_equal(
    CoTiMA:::.ctma_transformed_re_sd(alpha, tau, "exp(param)"),
    expected,
    tolerance = 1e-8
  )
})

test_that("pseudo-process T0 parameters are recognized as random intercept variation", {
  pars <- data.frame(
    matrix = c("T0MEANS", "T0MEANS", "T0VAR", "T0VAR", "T0VAR"),
    row = c(1L, 2L, 1L, 2L, 2L),
    col = c(1L, 1L, 1L, 1L, 2L),
    param = c("MeanV1", "MeanV1_cint", "T0var_V1",
              "Cov_T0eta1_cint1", "Cov_cint1"),
    stringsAsFactors = FALSE
  )

  ri <- CoTiMA:::.ctma_ri_pseudo_indices(pars)
  expect_identical(ri, 2L)
  expect_null(CoTiMA:::.ctma_ri_t0_interpretation(1L, 1L, ri))
  expect_match(
    CoTiMA:::.ctma_ri_t0_interpretation(2L, 2L, ri),
    "variance of a person-level random-intercept"
  )
  expect_match(
    CoTiMA:::.ctma_ri_t0_interpretation(2L, 1L, ri),
    "substantive initial state"
  )
})
