test_that("Stan pilot diagnostics summarize NUTS geometry", {
  chain1 <- cbind(
    divergent__ = c(0, 0),
    treedepth__ = c(9, 10),
    n_leapfrog__ = c(511, 1023),
    accept_stat__ = c(0.8, 0.9)
  )
  chain2 <- cbind(
    divergent__ = c(0, 1),
    treedepth__ = c(10, 8),
    n_leapfrog__ = c(1023, 255),
    accept_stat__ = c(0.85, 0.95)
  )
  summary_matrix <- matrix(
    c(
      0, 1, 100, 1.01,
      0, 1, 50, 1.05,
      0, 1, 10, 1.20
    ),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(
      c("theta", "phi", "llrow[1]"),
      c("mean", "sd", "n_eff", "Rhat")
    )
  )

  out <- CoTiMA2:::.ctmaSummarizeStanPilot(
    list(chain1, chain2), summary_matrix,
    max_treedepth = 10L, elapsed_seconds = 20
  )

  expect_identical(out$chains, 2L)
  expect_identical(out$postwarmup_per_chain, 2L)
  expect_equal(out$seconds_per_postwarmup_iteration, 10)
  expect_equal(out$divergences, 1)
  expect_equal(out$treedepth_saturated, 2)
  expect_equal(out$treedepth_saturation_fraction, 0.5)
  expect_equal(out$largest_substantive_rhat, 1.05)
  expect_equal(out$smallest_substantive_neff, 50)
  expect_equal(out$substantive_parameters, 2)
})

test_that("study-coding benchmark validates controlled arguments early", {
  skip_if_not_installed("rstan")

  expect_error(
    ctmaStudyCodingBenchmark(list(), iter = 100L, warmup = 100L),
    "smaller than"
  )
  expect_error(
    ctmaStudyCodingBenchmark(
      list(), fitArgs = list(studyCoding = "dummy")
    ),
    "must not be supplied"
  )
  expect_error(
    ctmaStudyCodingBenchmark(
      list(), backendArgs = list(control = list(max_treedepth = 8L))
    ),
    "Do not repeat"
  )
})
