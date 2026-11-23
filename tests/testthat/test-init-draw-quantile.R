test_that("undefined random-effect correlation draws retain unavailable quantiles", {
  for(p in c(.025,.5,.975)) {
    expect_identical(.ctmaInitDrawQuantile(1:100,p),stats::quantile(1:100,p))
    expect_identical(.ctmaInitDrawQuantile(c(1,NA,2),p),NA_real_)
    expect_identical(.ctmaInitDrawQuantile(c(1,NaN,2),p),NA_real_)
    expect_identical(.ctmaInitDrawQuantile(c(1,Inf,2),p),NA_real_)
    expect_identical(.ctmaInitDrawQuantile(numeric(),p),NA_real_)
  }
})
