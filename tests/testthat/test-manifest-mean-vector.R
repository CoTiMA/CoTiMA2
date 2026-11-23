test_that("manifest mean vectors preserve mixed labels and fixed values", {
  x <- ctmaLabels(n.latent=2,n.manifest=2,manifestMeans=c("mmV1",0))
  expect_identical(x$manifestMeansParams,c("mmV1","0"))
  x <- ctmaLabels(n.latent=2,n.manifest=2,manifestMeans=c(0,.5))
  expect_identical(x$manifestMeansParams,c(0,.5))
  expect_identical(ctmaLabels(n.latent=2,n.manifest=2,manifestMeans="auto")$manifestMeansParams,"auto")
  expect_equal(ctmaLabels(n.latent=2,n.manifest=2,manifestMeans=0)$manifestMeansParams,0)
})
