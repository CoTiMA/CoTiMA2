test_that('retired public switches are absent throughout the package', {
  retired <- c('activateRPB', 'experimental')
  affected <- c('ctmaAllInvFit', 'ctmaBiG', 'ctmaBiGOMX', 'ctmaEqual',
    'ctmaFit', 'ctmaInit', 'ctmaLCS', 'ctmaOptimizeFit', 'ctmaOptimizeInit',
    'ctmaPlot', 'ctmaPower', 'ctmaPRaw', 'ctmaRedHet', 'ctmaSaveFile')
  ns <- asNamespace('CoTiMA2')
  for (name in affected) {
    fun <- get(name, ns)
    expect_false(any(names(formals(fun)) %in% retired), info = name)
    expect_false(any(all.names(body(fun)) %in% retired), info = name)
    expect_false('RPushbullet' %in% all.names(body(fun)), info = name)
  }
  expect_error(ctmaInit(activateRPB = FALSE), 'no longer supports')
  expect_error(ctmaFit(experimental = FALSE), 'no longer supports')
})

test_that('default pseudo-data generation retains the standard covariance target', {
  sigma <- matrix(c(1, .3, .3, 1), 2)
  set.seed(516)
  x <- ctmaPRaw(empCovMat = sigma, empN = 80)
  set.seed(516)
  y <- ctmaPRaw(empCovMat = sigma, empN = 80)
  expect_identical(x, y)
  expect_equal(sum(complete.cases(x$data)), 80)
  expect_equal(unname(cov(x$data, use = "complete.obs")), sigma, ignore_attr = TRUE, tolerance = 1e-10)
})
