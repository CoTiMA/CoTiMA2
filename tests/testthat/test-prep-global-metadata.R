test_that("global moderator dictionaries are not indexed as studies", {
 x<-list(n.studies=3L,studyNumbers=as.list(c(1,4,313)),empcovs=rep(list(diag(2)),3),
  moderatorLabels=c("a","b","c"),moderatorValues=list("continuous",1:3,1:2))
 y<-CoTiMA2:::.ctmaPrepExcludeStudies(x,c(1,4))
 expect_identical(y$moderatorLabels,x$moderatorLabels)
 expect_identical(y$moderatorValues,x$moderatorValues)
})
