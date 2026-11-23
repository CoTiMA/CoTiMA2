order_init <- function() {
  ids <- c(10,20,30)
  studies <- lapply(ids,function(i)list(originalStudyNo=i,studyNumber=i,
    sampleSize=i,delta_t=c(i,i+1),timePoints=3L))
  list(studyList=studies,studyFitList=as.list(ids+1),emprawList=as.list(ids+2),
    ind.mod.List=as.list(ids+3),modelResults=list(DRIFT=as.list(ids+4),CINT=list()),
    primaryStudyList=list(studyNumbers=as.list(ids),sampleSizes=as.list(ids),
      moderators=as.list(ids+5),emprawList=as.list(ids+2),n.studies=3L,
      moderatorLabels=list('a','b','c'),moderatorValues=list(1,2,3)),
    statisticsList=list(originalStudyNumbers=ids,allSampleSizes=ids,
      overallSampleSize=sum(ids),allTpoints=rep(3,3),allDeltas=1:6,lambda=diag(2)),
    summary=list(estimates=matrix(ids,ncol=1,dimnames=list(paste('Study No',ids),'b')),
      confidenceIntervals=matrix(ids,ncol=1,dimnames=list(paste('Study No',ids),'lo'))),
    studyDesign=data.frame(study=ids,nWaves=3L),n.studies=3L)
}
test_that('study order and all aligned caches follow original identities', {
  x<-order_init();before<-x;set.seed(11);seed<-.Random.seed
  y<-.ctmaFitStudyOrder(x,c(30,10,20))
  expect_identical(x,before);expect_identical(.Random.seed,seed)
  expect_equal(vapply(y$studyList,`[[`,numeric(1),'originalStudyNo'),c(30,10,20))
  expect_equal(unlist(y$studyFitList),c(31,11,21))
  expect_equal(unlist(y$emprawList),c(32,12,22))
  expect_equal(unlist(y$ind.mod.List),c(33,13,23))
  expect_equal(unlist(y$modelResults$DRIFT),c(34,14,24))
  expect_equal(unlist(y$primaryStudyList$moderators),c(35,15,25))
  expect_identical(y$primaryStudyList$moderatorLabels,x$primaryStudyList$moderatorLabels)
  expect_equal(y$summary$estimates[,1],c('Study No 30'=30,'Study No 10'=10,'Study No 20'=20))
  expect_equal(y$studyDesign$study,c(30,10,20))
  expect_equal(y$statisticsList$allSampleSizes,c(30,10,20))
  expect_identical(y$statisticsList$lambda,x$statisticsList$lambda)
})
test_that('selection happens once, preserves dimensions and recomputes totals', {
  x<-order_init();y<-.ctmaFitStudyOrder(x,20)
  expect_length(y$studyList,1);expect_length(y$emprawList,1)
  expect_length(y$primaryStudyList$sampleSizes,1)
  expect_equal(dim(y$summary$estimates),c(1,1))
  expect_equal(y$statisticsList$overallSampleSize,20)
  expect_equal(y$statisticsList$overallTpoints,3)
  expect_equal(y$statisticsList$allDeltas,c(20,21))
  expect_equal(y$n.studies,1);expect_equal(y$primaryStudyList$n.studies,1)
  # A retained full-size preparation is matched by its own IDs, not truncated positions.
  x$studyList<-x$studyList[c(1,3)];x$studyFitList<-x$studyFitList[c(1,3)]
  x$emprawList<-x$emprawList[c(1,3)];x$ind.mod.List<-x$ind.mod.List[c(1,3)]
  x$modelResults$DRIFT<-x$modelResults$DRIFT[c(1,3)]
  y<-.ctmaFitStudyOrder(x,c(30,10))
  expect_equal(unlist(y$primaryStudyList$studyNumbers),c(30,10))
  expect_equal(unlist(y$primaryStudyList$moderators),c(35,15))
})
test_that('ambiguous or unknown identities and misaligned caches are rejected', {
  x<-order_init()
  for(ids in list(c(10,10),NA_real_,99,numeric()))
    expect_error(.ctmaFitStudyOrder(x,ids),'unique, nonmissing')
  x$emprawList<-list(1)
  expect_error(.ctmaFitStudyOrder(x,10),'emprawList.*not aligned')
})

test_that('unnamed historical confidence-interval rows follow initialization order', {
  x<-order_init();rownames(x$summary$confidenceIntervals)<-NULL
  y<-.ctmaFitStudyOrder(x,c(30,10))
  expect_equal(unname(y$summary$confidenceIntervals[,1]),c(30,10))
  expect_equal(dim(y$summary$confidenceIntervals),c(2,1))
})
