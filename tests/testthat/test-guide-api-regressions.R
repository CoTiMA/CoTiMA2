test_that('invariant and power argument paths apply the retired-control guard',{
 original<-.ctmaRetiredOptimControls
 checker<-function(x,backend)stop('Unrecognised optimcontrol name(s): ',names(x))
 testthat::local_mocked_bindings(.ctmaRetiredOptimControls=function(args)original(args,checker))
 for(fun in list(.ctmaResolveAllInvStanArgs,.ctmaResolvePowerStanArgs)){
  z<-fun();expect_null(z$optimcontrol$parstepsAutoModel);expect_null(z$optimcontrol$groupFreeThreshold)
  expect_error(fun(overrides=list(optimcontrol=list(parstepsAutoModel=TRUE))),'not silently discarded')
 }
})
test_that('model comparison reads native numeric metadata without print labels',{
 f<-list(summary=list(minus2ll=-10,n.parameters=5))
 g<-list(summary=list(minus2ll=20,n.parameters=3))
 z<-ctmaCompFit(f,g)[[1]]
 expect_true(any(grepl('Diff_Minus2LL:.*30$',z)))
 expect_true(any(grepl('Model 1.*Model 2',z)))
 expect_error(ctmaCompFit(list(),g),'finite -2 log likelihoods')
 expect_error(.ctmaComparisonStats(list(summary=list(model=c('a','b'),minus2ll=1,n.parameters=2)),'x'),'labels')
})
test_that('legacy vector model comparisons retain pair counts',{
 f<-list(summary=list(model=c('a','b'),minus2ll=c(10,20),n.parameters=c(3,4)))
 g<-list(summary=list(model='c',minus2ll=30,n.parameters=2))
 expect_length(ctmaCompFit(f,g)[[1]],14L)
})
test_that('shuffle matches original identities and repairs subset prep caches',{
 ids<-c(10,20,30)
 studies<-lapply(ids,function(i)list(originalStudyNo=i,studyNumber=i,startValues=i/10,delta_t=i,sampleSize=i*10))
 init<-list(studyList=studies,primaryStudyList=list(studyNumbers=as.list(ids),startValues=as.list(1:6),deltas=as.list(ids),emprawList=as.list(1:6)),emprawList=as.list(ids),studyFitList=as.list(ids+1))
 set.seed(11);seed<-.Random.seed
 z<-.ctmaOptimizeShuffleInit(init,c(30,10,20))
 expect_identical(.Random.seed,seed)
 expect_equal(vapply(z$studyList,function(x)x$originalStudyNo,numeric(1)),c(30,10,20))
 expect_equal(unlist(z$primaryStudyList$startValues),c(3,1,2))
 expect_equal(unlist(z$emprawList),c(30,10,20))
 expect_equal(unlist(z$primaryStudyList$emprawList),c(30,10,20))
 expect_equal(unlist(z$studyFitList),c(31,11,21))
 expect_equal(unlist(init$primaryStudyList$startValues),1:6)
 expect_error(.ctmaOptimizeShuffleInit(init,c(10,10,30)),'permutation')
 init$studyFitList<-list(1)
 expect_error(.ctmaOptimizeShuffleInit(init,c(30,10,20)),'align studyFitList')
})

test_that("equality refits preserve the recorded model time scale",{
 expect_equal(.ctmaEqualTimeScale(list(argumentList=list(scaleTime=1/48),summary=list())),1/48)
 expect_equal(.ctmaEqualTimeScale(list(summary=list(scaledTime=.25))),.25)
 expect_equal(.ctmaEqualTimeScale(list(summary=list(scaleTime=.5))),.5)
 expect_equal(.ctmaEqualTimeScale(list()),1)
 expect_error(.ctmaEqualTimeScale(list(argumentList=list(scaleTime=0))),"invalid time scale")
})
