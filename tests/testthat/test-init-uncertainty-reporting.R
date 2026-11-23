test_that("explicit no-width flags override numerical matrix widths, not estimates", {
 p <- data.frame(matrix=c("DRIFT","DRIFT","DIFFUSION"),row=c(1,2,1),col=c(1,1,1),param=c("a","b","q"))
 tab <- data.frame(matrix=c("DRIFT","DRIFT","DIFFUSIONcov","T0cov"),row=c(1,2,1,1),col=1,Mean=1:4,sd=.1,`2.5%`=0,`50%`=1,`97.5%`=2,check.names=FALSE)
 f<-list(model=list(pars=p,stationary=TRUE),uncertainty=list(intervalcheck=list(unidentified="a")),resultsSummary=list(parmatrices=tab))
 s<-.ctmaInitMaskNativeSummary(f)
 expect_identical(s$parmatrices$Mean,tab$Mean)
 expect_true(is.na(s$parmatrices$sd[1]))
 expect_false(is.na(s$parmatrices$sd[2]))
 expect_false(is.na(s$parmatrices$sd[3]))
 expect_true(is.na(s$parmatrices$sd[4]))
 f$uncertainty$intervalcheck$unidentified<-"q"
 expect_true(is.na(.ctmaInitMaskNativeSummary(f)$parmatrices$sd[3]))
 f$uncertainty$intervalcheck<-NULL
 f$identifiability<-list(nweak=1,parameters="a")
 expect_identical(.ctmaInitMaskNativeSummary(f),f$resultsSummary)
})

test_that("wave screen counts observed unique occasions, not rows or global times", {
 d<-data.frame(id=rep(1:3,each=2),time=c(0,1,0,2,0,3),V1=1)
 m<-list(manifestNames="V1",pars=data.frame(matrix="MANIFESTMEANS",indvarying=TRUE))
 expect_match(.ctmaInitWaveDiagnostic(d,m),"two observed waves")
 expect_null(.ctmaInitWaveDiagnostic(d,list(manifestNames="V1")))
 expect_match(.ctmaInitWaveDiagnostic(d,list(manifestNames="V1"),"CINT"),"constraints")
 d<-rbind(d,d[1,]);expect_match(.ctmaInitWaveDiagnostic(d,m),"two observed waves")
 d<-rbind(d,data.frame(id=1,time=2,V1=NA_real_))
 expect_match(.ctmaInitWaveDiagnostic(d,m),"two observed waves")
 d$V1[nrow(d)]<-1
 expect_null(.ctmaInitWaveDiagnostic(d,m))
})

test_that("saved-summary refresh preserves points and likelihood with mixed uncertainty", {
 tab<-data.frame(matrix="DRIFT",row=c(1,2,1,2),col=c(1,1,2,2),Mean=1:4,sd=.1)
 f<-list(model=list(pars=data.frame(matrix="DRIFT",row=1,col=2,param="a")),
  uncertainty=list(intervalcheck=list(unidentified="a")),resultsSummary=list(parmatrices=tab))
 a<-matrix(as.character(1:9),1,dimnames=list("Study No 17",NULL))
 x<-list(n.studies=1,n.latent=2,studyFitList=list(f),summary=list(estimates=a,
  drift_estimates_original_time_scale=a,confidenceIntervals=matrix(1,1,24),
  drift_CI_original_time_scale=matrix(1,1,8),minus2ll=99))
 y<-.ctmaRefreshInitUncertainty(x)
 expect_identical(y$summary$minus2ll,99)
 expect_identical(y$summary$estimates[,seq(2,8,2)],x$summary$estimates[,seq(2,8,2)])
 expect_true(is.na(y$summary$estimates[1,7]))
 expect_false(is.na(y$summary$estimates[1,3]))
 expect_true(all(is.na(y$summary$confidenceIntervals[1,5:6])))
 expect_match(y$summary$identification_notes,"Study 17")
 expect_identical(.ctmaRefreshInitUncertainty(y),y)
})
