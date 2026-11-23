optimize_fixture <- function() {
  ids<-c(11,22,33);group<-rep(1:3,c(2,3,1));long<-data.frame(id=rep(seq_along(group),each=2),time=rep(c(0,1),length(group)))
  long$Y<-seq_len(nrow(long))/10
  long$TI1<-as.numeric(rep(group,each=2)==1);long$TI2<-as.numeric(rep(group,each=2)==2)
  p<-data.frame(matrix='DRIFT',row=1L,col=1L,param='a',value=NA_real_,
    indvarying=FALSE,TI1_effect=TRUE,TI2_effect=TRUE)
  model<-list(pars=p,TIpredNames=c('TI1','TI2'),n.TIpred=2L)
  spec<-list(parameter_table=transform(p,parnumber=1L),laplace=NULL,
    ti_effects=data.frame(parameter=c(1L,1L),predictor=1:2,coefficient=2:3),TIpredNames=c('TI1','TI2'))
  native<-structure(list(model_spec=spec,estimate=list(raw=c(-.3,.1,.2),loglik=-10),
    args=list(resolved=list(optimize=TRUE,priors=FALSE,intoverpop='auto')),
    optim=list(converged=FALSE,gradient_norm=1)),class=c('ctJuliaFit','ctFit'))
  structure(list(studyFitList=native,ctModel=model,data=long,
    studyList=lapply(ids,function(i)list(originalStudyNo=i)),
    argumentList=list(backend='ctsemjulia',optimize=TRUE,priors=FALSE,
      studyCoding='dummy',scaleTime=1/48,scaleTI=FALSE),
    summary=list(loglik=-10,minus2ll=20)),class='CoTiMAFit')
}
test_that('default gambling keeps legacy dispatch, defaults and RNG',{
  expect_true(eval(formals(ctmaOptimizeFit)$gambling))
  # Modern dispatch must never consume RNG on a legacy call, even if invalid.
  set.seed(13);before<-.Random.seed
  expect_error(ctmaOptimizeFit(reFits=1L),'ctmaFitFit|ctmaIniFit')
  expect_identical(.Random.seed,before)
  expect_error(ctmaOptimizeFit(reFits=1,candidateStudies=1),'gambling=FALSE')
})
test_that('systematic contexts preserve data and use original IDs',{
  f<-optimize_fixture();c<-.ctmaOptimizeContext(f)
  expect_identical(c$ids,c('11','22','33'))
  expect_identical(c$counts,c(2L,3L,1L));expect_identical(c$scaleTime,1/48)
  expect_error(.ctmaOptimizeContext(f,backend='stan'),'retain')
  f$argumentList$priors<-TRUE;expect_error(.ctmaOptimizeContext(f),'ML only')
})
test_that('exclusion and reference plans are explicit, immutable and deduplicated',{
  c<-.ctmaOptimizeContext(optimize_fixture())
  p<-.ctmaOptimizePlan(c,c(11,22),c('single','blocks'),list(c(11)),reference='largestN')
  expect_identical(p[[1]]$reference,'33');expect_identical(p[[2]]$reference,'22')
  expect_length(p,4L) # current, largest, omit11, omit22; duplicate block omitted
  expect_error(.ctmaOptimizePlan(c,999,'single'),'original study IDs')
  expect_error(.ctmaOptimizePlan(c,NULL,'single'),'Explicit')
  expect_error(.ctmaOptimizePlan(c,c(11,22,33),'all'),'fewer than two')
  expect_error(.ctmaOptimizePlan(c,11,'single',maxCases=1),'maxCases')
  expect_identical(c$ids,c('11','22','33'))
})
test_that('rebasing preserves raw predictors and does not mutate observations',{
  c<-.ctmaOptimizeContext(optimize_fixture())
  p<-.ctmaOptimizePlan(c,reference='largestN')[[2]]
  x<-.ctmaOptimizeCase(c,p)
  expect_identical(vapply(x$studies,function(s)as.character(s$originalStudyNo),character(1)),c('11','33','22'))
  expect_equal(x$targets[['base:a']],-.1)
  expect_equal(x$targets[['ti:a|TI1']],-.1)
  expect_equal(x$targets[['ti:a|TI2']],-.2)
  expect_identical(x$data[,c('id','time','Y')],c$data[,c('id','time','Y')])
  p<-.ctmaOptimizePlan(c,11,'single')[[2]];x<-.ctmaOptimizeCase(c,p)
  expect_equal(sort(unique(x$data$id)),3:6)
  expect_length(x$model$TIpredNames,1);expect_equal(x$targets[['base:a']],-.3)
  expect_equal(x$targets[['ti:a|TI1']],.2)
  expect_false('TI2'%in%names(x$data))
})
test_that('partial constraints are not silently freed by reference changes',{
  f<-optimize_fixture();f$ctModel$pars$TI2_effect<-FALSE
  # Inactive TI slot must not exist in the original fitted layout.
  f$studyFitList$model_spec$ti_effects<-f$studyFitList$model_spec$ti_effects[1,,drop=FALSE]
  f$studyFitList$estimate$raw<-f$studyFitList$estimate$raw[1:2]
  c<-.ctmaOptimizeContext(f)
  expect_error(.ctmaOptimizeCase(c,.ctmaOptimizePlan(c,reference='largestN')[[2]]),'partially constrained')
})
test_that('authoritative Stan TI indices override misleading covariance names',{
  p<-data.frame(matrix=c('DRIFT','DIFFUSION'),row=1L,col=1L,param=c('a','q'),value=NA_real_)
  d<-list(nparams=2L,nindvarying=0L,nindvaryingoffdiagonals=0L,
    matsetup=cbind(1L,1L,1:2),TIPREDEFFECTsetup=matrix(c(1,3,2,4),2),indvaryingindex=integer())
  n<-structure(list(ctstanmodel=list(pars=p),ctstanmodelbase=list(TIpredNames=c('TI1','TI2')),
    standata=d,stanfit=list(rawest=1:6)),class='ctStanFit')
  map<-.ctmaOptimizeLayout(n)
  expect_identical(map$keys,c('base:a','base:q','ti:a|TI1','ti:a|TI2','ti:q|TI1','ti:q|TI2'))
})
test_that('pilots never certify uncertainty and failures remain explicit',{
  good<-list(assessment=list(loglik=-10,bad=FALSE,converged=TRUE,uncertaintyValidated=TRUE))
  expect_identical(.ctmaOptimizeStatus(good,TRUE),'screened_inconclusive')
  expect_identical(.ctmaOptimizeStatus(good,FALSE),'numerically_checked')
  expect_identical(.ctmaOptimizeStatus(list(error='cap',timedOut=TRUE),TRUE),'timeout')
  expect_error(.ctmaOptimizeControls(list(bad=TRUE)),'Unknown')
  expect_error(.ctmaOptimizeControls(list(maxiter=0)),'maxiter')
})
test_that('systematic orchestration keeps full-data best and checkpoints all warnings',{
  f<-optimize_fixture();dir<-tempfile();on.exit(unlink(dir,recursive=TRUE))
  calls<-0L
  local_mocked_bindings(.ctmaOptimizeRun=function(payload,timeout,logfile){
    calls<<-calls+1L
    ll<-if(payload$case$case$fullData)-11 else -1
    fit<-f;fit$summary<-list(loglik=ll,minus2ll=-2*ll)
    list(fit=fit,point=c('base:a'=-.3,'ti:a|TI1'=.1,'ti:a|TI2'=.2),
      assessment=list(loglik=ll,converged=FALSE,maxGradient=.001,
        stationary=FALSE,weakDirections=2L,repaired=TRUE,bad=FALSE,uncertaintyValidated=FALSE),
      startLoglik=-12,warnings=c('first','second'))
  })
  x<-ctmaOptimizeFit(ctmaFitFit=f,gambling=FALSE,candidateStudies=11,
    exclusion='single',quickCheck=FALSE,reFits=1L,iter=10L,verbose=0,
    checkpointDirectory=dir,timeBudget=30)
  expect_identical(x$bestFit,f) # original -10 beats every full-data restart -11
  expect_null(x$bestConvergedFit)
  expect_equal(calls,2L);expect_length(x$warnings[[1]]$messages,2L)
  expect_equal(x$caseDiagnostics$loglik,c(-11,-1))
  expect_equal(x$all_minus2ll,22) # excluded-study criterion NEVER enters winner vector
  oldCalls<-calls
  y<-ctmaOptimizeFit(ctmaFitFit=f,gambling=FALSE,candidateStudies=11,
    exclusion='single',quickCheck=FALSE,reFits=1L,iter=10L,verbose=0,
    checkpointDirectory=dir,timeBudget=30,resume=TRUE)
  expect_equal(calls,oldCalls);expect_identical(y$caseDiagnostics$status,x$caseDiagnostics$status)
  expect_error(ctmaOptimizeFit(ctmaFitFit=f,gambling=FALSE,candidateStudies=22,
    exclusion='single',quickCheck=FALSE,reFits=1L,iter=10L,
    checkpointDirectory=dir,timeBudget=30,resume=TRUE),'manifest differs')
})
test_that('systematic failures and quick checks do not replace the full-data fit',{
  f<-optimize_fixture();dir<-tempfile();on.exit(unlink(dir,recursive=TRUE))
  local_mocked_bindings(.ctmaOptimizeRun=function(...)list(error='bad pilot',timedOut=TRUE,warnings=character()))
  x<-ctmaOptimizeFit(ctmaFitFit=f,gambling=FALSE,checkpointDirectory=dir,timeBudget=10,verbose=0)
  expect_identical(x$bestFit,f);expect_identical(x$caseDiagnostics$status,'timeout')
  expect_length(x$errors,1L)
})
test_that('scaled dummy rebasing preserves each retained study predictor',{
 f<-optimize_fixture();subject<-!duplicated(f$data$id)
 old<-as.matrix(f$data[,c('TI1','TI2')]);means<-colMeans(old[subject,,drop=FALSE]);sds<-apply(old[subject,,drop=FALSE],2,sd)
 f$data[,c('TI1','TI2')]<-sweep(sweep(old,2,means,'-'),2,sds,'/')
 f$argumentList$scaleTI<-TRUE
 context<-.ctmaOptimizeContext(f)
 for(p in .ctmaOptimizePlan(context,11,'single',reference='largestN')) {
   z<-.ctmaOptimizeCase(context,p)
   order<-vapply(z$studies,function(x)as.character(x$originalStudyNo),character(1))
   got<-vapply(order,function(id){
     at<-which(context$ids[context$group]%in%id)[1]
     if(!p$fullData)at<-which(z$data$id==context$data$id[at])[1]
     z$targets[['base:a']]+sum(vapply(z$model$TIpredNames,function(n)z$data[[n]][at]*z$targets[[paste0('ti:a|',n)]],numeric(1)))
   },numeric(1))
   want<-as.numeric(context$design[match(order,context$ids),,drop=FALSE]%*%context$layout$raw)
   expect_equal(unname(got),want,tolerance=1e-12)
 }
})
test_that('two-study contexts and excluded references are explicit',{
 f<-optimize_fixture();c<-.ctmaOptimizeContext(f)
 p<-.ctmaOptimizePlan(c,33,'single')[[2]]
 expect_identical(p$reference,'22')
 z<-.ctmaOptimizeCase(c,p)
 f$data<-z$data;f$ctModel<-z$model;f$studyList<-z$studies
 f$studyFitList$model_spec$TIpredNames<-'TI1'
 f$studyFitList$model_spec$ti_effects<-data.frame(parameter=1L,predictor=1L,coefficient=2L)
 f$studyFitList$estimate$raw<-unname(unlist(z$targets))
 expect_identical(.ctmaOptimizeContext(f)$ids,c('11','22'))
})
test_that('converged baseline needs no summary criterion for comparison',{
 f<-optimize_fixture();f$summary<-NULL;f$studyFitList$optim$converged<-TRUE
 dir<-tempfile();on.exit(unlink(dir,recursive=TRUE))
 local_mocked_bindings(.ctmaOptimizeRun=function(...)list(fit=f,
  assessment=list(loglik=-11,converged=TRUE,maxGradient=.001,stationary=NA,
   weakDirections=0L,repaired=NA,bad=FALSE,uncertaintyValidated=FALSE),warnings=character()))
 x<-ctmaOptimizeFit(ctmaFitFit=f,gambling=FALSE,quickCheck=FALSE,checkpointDirectory=dir,verbose=0)
 expect_identical(x$bestConvergedFit,f)
 expect_identical(x$bestFit,f)
})
test_that('Stan fit=FALSE layouts are adapted before scientific key extraction',{
 p<-data.frame(matrix='DRIFT',row=1L,col=1L,param='a',value=NA_real_)
 d<-list(nparams=1L,nindvarying=0L,nindvaryingoffdiagonals=0L,
  matsetup=matrix(c(1L,1L,1L),nrow=1),TIPREDEFFECTsetup=matrix(1L),indvaryingindex=integer())
 prepared<-list(ctstanmodel=list(pars=p),ctstanmodelbase=list(TIpredNames='TI1'),standata=d)
 expect_identical(.ctmaOptimizePreparedLayout(prepared,'stan')$keys,c('base:a','ti:a|TI1'))
})
test_that('resume fingerprints include the actual starting point',{
 f<-optimize_fixture();a<-.ctmaOptimizeContext(f)$fingerprint
 f$studyFitList$estimate$raw[1]<--.4
 expect_false(identical(a,.ctmaOptimizeContext(f)$fingerprint))
})
test_that('systematic controls preserve caller RNG and honor full draw counts',{
 f<-optimize_fixture();dir<-tempfile();on.exit(unlink(dir,recursive=TRUE))
 set.seed(29);before<-.Random.seed;seen<-NULL
 local_mocked_bindings(.ctmaOptimizeRun=function(payload,...){seen<<-payload$finishsamples
  list(error='intentional test failure',warnings=character())})
 ctmaOptimizeFit(ctmaFitFit=f,gambling=FALSE,quickCheck=FALSE,
  finishsamples=777,checkpointDirectory=dir,verbose=0)
 expect_identical(.Random.seed,before);expect_equal(seen,777)
})
test_that('published Stan iteration limits fall back explicitly to wall-time caps',{
 expect_warning(x<-.ctmaOptimizeSetIterations(list(maxiter=9),'stan',10,function(estonly=FALSE)NULL),'wall-time cap')
 expect_null(x$maxiter)
 expect_equal(.ctmaOptimizeSetIterations(list(),'stan',10,function(maxiter=1)NULL)$maxiter,10)
 expect_equal(.ctmaOptimizeSetIterations(list(),'ctsemjulia',10)$maxiter,10)
})
test_that('legacy mode does not silently ignore new budget/checkpoint controls',{
 expect_error(ctmaOptimizeFit(timeBudget=5),'gambling=FALSE')
 expect_error(ctmaOptimizeFit(checkpointDirectory=tempfile()),'gambling=FALSE')
})
test_that('auto references rank actual sample size and skip current reference',{
 c<-.ctmaOptimizeContext(optimize_fixture())
 p<-.ctmaOptimizePlan(c,references='AUTO')
 a<-attr(p,'automaticSelection')
 expect_identical(a$referenceCandidates,c('22','11'))
 expect_length(p,3L)
 expect_false(any(a$studyDiagnostics$nominated_exclusion))
})
test_that('auto exclusions require localized stored flags, not small N',{
 c<-.ctmaOptimizeContext(optimize_fixture())
 p<-.ctmaOptimizePlan(c,candidates='auto',exclusion=c('single','all'))
 expect_length(p,1L)
 expect_length(attr(p,'automaticSelection')$candidateStudies,0L)
 c$native$identifiability<-list(directions=list(
  list(parameters=c('a','rawtipredeffect_a_TI2'),loadings=c(.8,.2)),
  list(parameters='rawtipredeffect_a_TI1',loadings=1)))
 p<-.ctmaOptimizePlan(c,candidates='auto',references='auto',exclusion=c('single','all'))
 a<-attr(p,'automaticSelection')
 expect_identical(a$candidateStudies,'11')
 expect_identical(a$referenceCandidates,'22')
 expect_equal(a$studyDiagnostics$localized_weak_directions,c(1,0,0))
 expect_equal(sum(vapply(p,function(x)!x$fullData,logical(1))),1)
})
test_that('auto primary diagnostics are aligned by original ID',{
 f<-optimize_fixture()
 initial<-list(studyList=rev(f$studyList),studyFitList=list(
  list(optim=list(converged=TRUE)),list(optim=list(converged=FALSE)),list(optim=list(converged=TRUE))))
 c<-.ctmaOptimizeContext(f,initial)
 a<-.ctmaOptimizeAutoCandidates(c,'auto','auto')
 expect_identical(a$candidateStudies,'22')
 expect_identical(a$referenceCandidates,'11')
 expect_equal(which(a$studyDiagnostics$primary_nonconverged),2)
})
test_that('auto planning handles partial constraints and rejects spurious labels',{
 c<-.ctmaOptimizeContext(optimize_fixture());c$model$pars$TI2_effect<-FALSE
 c$native$identifiability<-list(directions=list(
  list(parameters='unrelated_TI1',loadings=1),
  list(parameters='rawtipredeffect_a_TI1',loadings=NA_real_)))
 a<-.ctmaOptimizeAutoCandidates(c,'auto','auto')
 expect_length(a$referenceCandidates,0L);expect_length(a$candidateStudies,0L)
 expect_match(paste(a$notes,collapse=' '),'partial')
})
test_that('automatic nomination reports do not silently request exclusions',{
 f<-optimize_fixture();f$studyFitList$identifiability<-list(directions=list(
   list(parameters='rawtipredeffect_a_TI1',loadings=1)))
 dir<-tempfile();on.exit(unlink(dir,recursive=TRUE))
 called<-list()
 local_mocked_bindings(.ctmaOptimizeRun=function(payload,...){
  called[[length(called)+1L]]<<-payload$case$case
  list(error='planning probe only',warnings=character())})
 x<-ctmaOptimizeFit(ctmaFitFit=f,gambling=FALSE,candidateStudies='auto',referenceCandidates='auto',
  checkpointDirectory=dir,verbose=0)
 expect_identical(x$automaticSelection$candidateStudies,'11')
 expect_true(all(vapply(called,function(x)x$fullData,logical(1))))
 expect_true(file.exists(file.path(dir,'automaticStudyDiagnostics.csv')))
})
test_that('empty auto nominations do not hide invalid exclusion syntax',{
 c<-.ctmaOptimizeContext(optimize_fixture())
 expect_error(.ctmaOptimizePlan(c,candidates='auto',exclusion=c('none','single')),'exclusion must')
 expect_error(.ctmaOptimizePlan(c,candidates='auto',exclusion=NULL),'exclusion must')
})
test_that('planOnly exposes auto nominations without workers or writes',{
 f<-optimize_fixture();dir<-tempfile()
 local_mocked_bindings(.ctmaOptimizeRun=function(...)stop('must not fit'))
 x<-ctmaOptimizeFit(ctmaFitFit=f,gambling=FALSE,referenceCandidates='auto',
  candidateStudies='auto',exclusion='single',planOnly=TRUE,checkpointDirectory=dir)
 expect_true(x$planOnly);expect_equal(nrow(x$casePlan),3L)
 expect_identical(summary(x),x$casePlan);expect_false(dir.exists(dir))
 expect_length(x$caseDiagnostics,0L)
 expect_error(ctmaOptimizeFit(planOnly=TRUE),'gambling=FALSE')
})
test_that('weak zero is not asserted when curvature was not computed',{
 native<-structure(list(estimate=list(raw=c(1,2),loglik=-10),optim=list(converged=TRUE),
  identifiability=list(nweak=0L),uncertainty=list(certification=list(certified=TRUE))),class='ctJuliaFit')
 a<-.ctmaOptimizeWorkerAssess(native)
 expect_false(a$curvatureAssessed);expect_true(is.na(a$weakDirections));expect_false(a$uncertaintyValidated)
 native$uncertainty$hessian<--diag(2)
 native$estimate$cov<-diag(2)
 attr(native$estimate$cov,'ctOptimCovFromHessian')<-list(usedNullProjection=FALSE,usedNearPD=FALSE,infoNearPD=FALSE,nullDirections=0L)
 a<-.ctmaOptimizeWorkerAssess(native)
 expect_true(a$curvatureAssessed);expect_equal(a$weakDirections,0);expect_true(a$uncertaintyValidated)
 native$uncertainty$hessian[1,1]<-NA_real_
 expect_true(is.na(.ctmaOptimizeWorkerAssess(native)$weakDirections))
})
test_that('legacy summary zeroes without assessment provenance are unknown',{
 x<-structure(list(caseDiagnostics=data.frame(weak_directions=c(0L,2L))),class='ctmaOptimizeSystematic')
 s<-summary(x)
 expect_true(is.na(s$weak_directions[1]));expect_equal(s$weak_directions[2],2)
 expect_true(all(is.na(s$curvature_assessed)))
 expect_identical(x$caseDiagnostics$weak_directions,c(0L,2L))
 x$caseDiagnostics$curvature_assessed<-c(TRUE,FALSE)
 expect_identical(summary(x)$weak_directions,c(0L,2L))
})
