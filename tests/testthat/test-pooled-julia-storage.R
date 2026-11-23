test_that('pooled storage preserves engine object and metadata without nesting primary fits', {
  native <- structure(list(coefficients=list(raw=c(a=.1,b=.2)),
    information=list(covariance=diag(2)), optimization=list(loglik=-12)),class='ctmaJuliaFit')
  setup <- list(studyList=list(list(originalStudyNo=1)),emprawList=list(data.frame(V1_T0=1,V1_T1=2,dT1=3)))
  out <- CoTiMA2:::.ctmaWrapPooledJuliaFit(native,setup)
  expect_s3_class(out,'CoTiMAFit')
  expect_false(inherits(out,'ctmaJuliaFit'))
  expect_identical(out$studyFitList,native)
  expect_identical(out$emprawList,setup$emprawList)
  expect_identical(coef(out),coef(native))
  expect_identical(vcov(out),vcov(native))
  expect_identical(unserialize(serialize(out,NULL)),out)
})
