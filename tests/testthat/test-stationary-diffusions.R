test_that("stationary options are strict and retain explicit legacy settings", {
  expect_equal(.ctmaStationaryOptions(FALSE,FALSE),list(stationary=FALSE,determinedDiffusions=FALSE))
  expect_true(.ctmaStationaryOptions(FALSE,FALSE,TRUE,TRUE)$stationary)
  expect_error(.ctmaStationaryOptions(FALSE,FALSE,FALSE,TRUE),"Conflicting")
  for(x in list(NA,1,c(TRUE,FALSE),NULL))expect_error(.ctmaStationaryOptions(x,FALSE),"stationary")
  expect_error(.ctmaStationaryOptions(TRUE,TRUE,optimize=FALSE),"maximum likelihood")
  expect_error(.ctmaStationaryOptions(TRUE,TRUE,priors=TRUE),"maximum likelihood")
  expect_error(.ctmaStationaryOptions(FALSE,TRUE),"requires stationary")
  expect_error(.ctmaStationaryOptions(TRUE,TRUE,inits=1:2),"coordinates")
  expect_error(.ctmaStationaryOptions(TRUE,FALSE,backend="julia"),"bundled")
})

.stationary_test_model <- function(n=2) {
  ctsem::ctModel(type="ct",n.latent=n,n.manifest=n,LAMBDA=diag(n),
    DRIFT=matrix(if(n==1)"a" else c("a","c","b","d"),n,n),
    DIFFUSION=if(n==1)matrix("q") else matrix(c("q1","q12",0,"q2"),2,2),
    T0VAR="auto",T0MEANS=matrix(0,n,1),CINT=matrix(0,n,1),MANIFESTMEANS=matrix(0,n,1),
    MANIFESTVAR=diag(.1,n))
}

test_that("determined diffusion defaults do not change any model field",{
  m<-.stationary_test_model()
  expect_identical(.ctmaDeterminedDiffusionModel(m,FALSE),m)
})


.stationary_transformed_values <- function(m, raw) {
  out<-as.list(raw)
  pars<-m$pars[m$pars$matrix=="PARS",,drop=FALSE]
  for(i in seq_len(nrow(pars))){
    e<-list(param=raw[[pars$param[i]]],log1p_exp=function(x)log1p(exp(x)))
    out[[pars$param[i]]]<-eval(parse(text=pars$transform[i]),e)
  }
  out
}

test_that("derived stationary matrices satisfy Lyapunov identity away from the initial point",{
  set.seed(3)
  for(n in 1:2){
    m<-suppressMessages(.ctmaDeterminedDiffusionModel(.stationary_test_model(n),TRUE))
    expect_equal(m$covmattransform,"cholesky")
    a<-m$determinedDiffusions$auxiliary
    expect_length(a,if(n==1)2 else 7)
    expect_false(anyDuplicated(paste(m$pars$matrix,m$pars$row,m$pars$col))>0)
    for(i in 1:8){
      raw<-setNames(rnorm(length(a),sd=.3),a)
      values<-.stationary_transformed_values(m,raw)
      calc<-function(mat){
        p<-m$pars[m$pars$matrix==mat,];M<-matrix(0,n,n)
        for(j in seq_len(nrow(p)))M[p$row[j],p$col[j]]<-if(is.na(p$param[j]))p$value[j] else eval(parse(text=p$param[j]),values)
        M
      }
      A<-calc("DRIFT");L<-calc("DIFFUSION");LS<-calc("T0VAR");S<-tcrossprod(LS);Q<-tcrossprod(L)
      expect_equal(A%*%S+S%*%t(A)+Q,matrix(0,n,n),tolerance=1e-10)
      expect_gt(min(eigen(Q,symmetric=TRUE)$values),0)
      expect_lt(max(Re(eigen(A)$values)),0)
    }
  }
})

test_that("unmapped model constraints are rejected before rewriting",{
  m<-.stationary_test_model()
  bad<-m;bad$pars$value[bad$pars$matrix=="DRIFT"&bad$pars$row!=bad$pars$col]<-0
  bad$pars$param[bad$pars$matrix=="DRIFT"&bad$pars$row!=bad$pars$col]<-NA
  expect_error(.ctmaDeterminedDiffusionModel(bad,TRUE),"distinct free")
  bad<-m;bad$pars$z_effect<-bad$pars$matrix=="DIFFUSION"
  expect_error(.ctmaDeterminedDiffusionModel(bad,TRUE),"invariance")
  bad<-m;bad$pars$indvarying[bad$pars$matrix=="DRIFT"]<-TRUE
  expect_error(.ctmaDeterminedDiffusionModel(bad,TRUE),"random effects")
  bad<-m;bad$pars$param[bad$pars$matrix=="T0MEANS"]<-"t0free";bad$pars$indvarying[bad$pars$matrix=="T0MEANS"]<-TRUE
  expect_error(.ctmaDeterminedDiffusionModel(bad,TRUE),"indVaryingT0")
  bad<-m;bad$n.latent<-4
  expect_error(.ctmaDeterminedDiffusionModel(bad,TRUE),"one or two")
})

test_that("public arguments are appended and default to FALSE",{
  for(fun in list(ctmaInit,ctmaFit)){
    expect_identical(formals(fun)$stationary,FALSE)
    expect_identical(formals(fun)$determinedDiffusions,FALSE)
  }
})

test_that("stationary result metadata preserves the legacy default result",{
  x<-list(summary=list(minus2ll=4),argumentList=list(fit=TRUE))
  expect_identical(.ctmaStationaryResult(x,FALSE,FALSE),x)
  y<-.ctmaStationaryResult(x,TRUE,TRUE)
  expect_true(y$argumentList$stationary)
  expect_true(y$argumentList$determinedDiffusions)
  expect_match(y$summary$stationary_specification$note,"does not fix")
})

test_that("generated model text preserves the determined covariance construction",{
  m<-suppressMessages(.ctmaDeterminedDiffusionModel(.stationary_test_model(),TRUE))
  e<-new.env(parent=globalenv())
  eval(parse(text=.ctmaModelText(m)),e)
  expect_equal(e$m$pars,m$pars)
  expect_identical(e$m$covmattransform,"cholesky")
})

test_that("native primary-study dispatcher forwards stationary", {
  seen<-NULL
  local_mocked_bindings(
    .ctmaRequireCtsemJuliaLaplace=function()list(available=TRUE),
    .ctmaNativeCtsemJuliaFit=function(args){seen<<-args;list()},
    .ctmaNativeCtsemJuliaSummary=function(...)list(),
    .ctmaCtsemJuliaMode=function(...)"direct",
    .ctmaAnnotateCtsemJuliaResult=function(result,...)result
  )
  .ctmaInitCtsemJuliaStudyFit(data.frame(id=1,time=0),list(pars=data.frame(indvarying=FALSE)),
    same_initial_times=FALSE,backend_args=list(stationary=TRUE))
  expect_true(seen$stationary)
})

test_that("Stan determined covariance fallback is explicit and exact at initialization",{
  m<-suppressMessages(.ctmaDeterminedDiffusionModel(.stationary_test_model(),TRUE))
  a<-suppressMessages(.ctmaStationaryFitArgs(list(model=m,stationary=TRUE,backend="stan")))
  expect_false(a$stationary)
  expect_equal(a$model$pars[a$model$pars$matrix=="T0VAR",],m$pars[m$pars$matrix=="T0VAR",])
  expect_true(all(a$model$pars$value[a$model$pars$matrix=="T0MEANS"]==0))
  bad<-m;bad$pars$value[bad$pars$matrix=="CINT"]<-1
  expect_error(.ctmaStationaryFitArgs(list(model=bad,stationary=TRUE,backend="stan")),"CINT")
  x<-list(stationary=FALSE,model=m)
  expect_identical(.ctmaStationaryFitArgs(x),x)
})

test_that("untested binary and random CINT structures are guarded",{
  m<-.stationary_test_model()
  m$manifesttype<-c(1,0)
  expect_error(.ctmaDeterminedDiffusionModel(m,TRUE),"Gaussian")
  m<-.stationary_test_model()
  m$pars$param[m$pars$matrix=="CINT"]<-c("c1","c2")
  m$pars$value[m$pars$matrix=="CINT"]<-NA
  m$pars$indvarying[m$pars$matrix=="CINT"]<-TRUE
  expect_error(.ctmaDeterminedDiffusionModel(m,TRUE),"random CINT")
})

test_that("original drift transforms and signs are preserved exactly",{
  original<-.stationary_test_model()
  m<-suppressMessages(.ctmaDeterminedDiffusionModel(original,TRUE))
  dr<-original$pars[original$pars$matrix=="DRIFT",]
  at<-match(dr$param,m$pars$param[m$pars$matrix=="PARS"])
  expect_equal(m$pars$transform[m$pars$matrix=="PARS"][at],dr$transform)
  expect_true(m$determinedDiffusions$drift_constraints_preserved)
  set.seed(54)
  for(k in 1:30){
    raw<-setNames(rnorm(length(m$determinedDiffusions$auxiliary),sd=.5),
      m$determinedDiffusions$auxiliary)
    values<-.stationary_transformed_values(m,raw)
    for(i in seq_len(nrow(dr))){
      e<-list(param=raw[[dr$param[i]]],log1p_exp=function(x)log1p(exp(x)))
      expected<-eval(parse(text=dr$transform[i]),e)
      actual<-eval(parse(text=m$pars$param[m$pars$matrix=="DRIFT"&
        m$pars$row==dr$row[i]&m$pars$col==dr$col[i]]),values)
      expect_identical(actual,expected)
      if(dr$row[i]==dr$col[i])expect_lt(actual,-1e-6)
    }
  }
})

test_that("native stationary priors are explicit and proper for sampling", {
  expect_message(x <- .ctmaStationaryOptions(TRUE, TRUE, backend="ctsemjulia", optimize=FALSE, priors=TRUE), "stationary raw coordinates")
  expect_true(x$determinedDiffusions)
  expect_error(.ctmaStationaryOptions(TRUE, TRUE, backend="ctsemjulia", optimize=FALSE, priors=FALSE), "priors=TRUE")
  expect_message(.ctmaStationaryOptions(TRUE, TRUE, backend="ctsemjulia", optimize=TRUE, priors=TRUE), "stationary raw coordinates")
})


test_that("manifest RI covariance construction is retained independently of process coordinates", {
  original<-.stationary_test_model()
  ix<-original$pars$matrix=="MANIFESTMEANS"
  original$pars$param[ix]<-c("ri1","ri2")
  original$pars$value[ix]<-NA_real_
  original$pars$indvarying[ix]<-TRUE
  original$covmattransform<-"z"
  m<-suppressMessages(.ctmaDeterminedDiffusionModel(original,TRUE))
  expect_identical(m$covmattransform,"z")
  expect_true(m$determinedDiffusions$random_covariance_construction_preserved)
  actual<-m$pars[m$pars$matrix=="MANIFESTMEANS",];expected<-original$pars[ix,]
  rownames(actual)<-rownames(expected)<-NULL
  expect_equal(actual,expected)
  raw<-m$determinedDiffusions$suggested_raw
  raw[c("b","c")]<-c(-.1,.08)
  values<-.stationary_transformed_values(m,raw)
  calc<-function(mat){p<-m$pars[m$pars$matrix==mat,];M<-matrix(0,2,2)
    for(j in seq_len(nrow(p)))M[p$row[j],p$col[j]]<-if(is.na(p$param[j]))p$value[j] else eval(parse(text=p$param[j]),values)
    M}
  covariance<-function(B){r<-tanh(B[2,1]);D<-diag(diag(B));D%*%matrix(c(1,r,r,1),2)%*%D}
  A<-calc("DRIFT");Q<-covariance(calc("DIFFUSION"));S<-covariance(calc("T0VAR"))
  expect_equal(A%*%S+S%*%t(A)+Q,matrix(0,2,2),tolerance=1e-10)
  expect_gt(min(eigen(Q,symmetric=TRUE)$values),0)
  bad<-original;bad$covmattransform<-"rawcorr"
  expect_error(.ctmaDeterminedDiffusionModel(bad,TRUE),"safeguards")
})
test_that("saturated study covariance slopes transfer only to auxiliary covariance coordinates", {
 m<-.stationary_test_model()
 m$TIpredNames<-c("TI1","TI2")
 for(v in c("TI1_effect","TI2_effect"))m$pars[[v]]<-m$pars$matrix=="DIFFUSION" & !is.na(m$pars$param)
 design<-rbind(c(0,0),c(1,0),c(0,1))
 out<-suppressMessages(.ctmaDeterminedDiffusionModel(m,TRUE,
   c("TI1_effect","TI2_effect"),design))
 p<-out$pars
 for(v in c("TI1_effect","TI2_effect")) {
  expect_equal(sum(p[[v]]),3)
  expect_true(all(p$matrix[p[[v]]]=="PARS"))
  expect_false(any(p[[v]][p$param %in% c("a","b","c","d")]))
 }
 expect_equal(out$determinedDiffusions$study_covariances,"fully study-specific")
 bad<-m;bad$pars$TI1_effect[bad$pars$param %in% "q1"]<-FALSE
 expect_error(.ctmaDeterminedDiffusionModel(bad,TRUE,c("TI1_effect","TI2_effect"),design),"Partial covariance")
 bad<-m;bad$pars$TI1_effect[bad$pars$param %in% "a"]<-TRUE
 expect_error(.ctmaDeterminedDiffusionModel(bad,TRUE,c("TI1_effect","TI2_effect"),design),"common drift")
 expect_error(.ctmaDeterminedDiffusionModel(m,TRUE,c("TI1_effect","TI2_effect"),matrix(0,3,2)),"full rank")
 expect_error(.ctmaDeterminedDiffusionModel(m,TRUE,"TI1_effect",design[,1,drop=FALSE]),"invariance")
})

test_that("public pooled stationary setup supports dummy and weighted covariance coding", {
 # Setup-only mapping test: no native engine is requested or exercised.
 # Published ctsem constructs the legacy covariance layout. Emulate the
 # modern constructor contract only in this dry-run API test; real legacy
 # layouts remain rejected by the separate runtime safeguard tests.
 native_model <- ctsem::ctModel
 local_mocked_bindings(ctModel=function(...) {
   m <- native_model(...); m$covmattransform <- "z"; m
 }, .package="ctsem")
 local_mocked_bindings(.ctmaRequireCtsemJuliaLaplace=function()list(available=TRUE))
 e<-new.env();utils::data("CoTiMAstudyList_3",package="CoTiMA2",envir=e)
 utils::data("CoTiMAInitFit_3",package="CoTiMA2",envir=e)
 fresh<-suppressMessages(suppressWarnings(ctmaInit(primaryStudies=e$CoTiMAstudyList_3,
  n.latent=2,fit=FALSE,silentOverwrite=TRUE,activeDirectory=paste0(tempdir(),"/"))))
 init<-e$CoTiMAInitFit_3;init$emprawList<-fresh$data
 init$studyFitList<-lapply(fresh$data,function(d)list(empraw=d))
 for(coding in c("dummy","weighted")) {
  out<-suppressMessages(suppressWarnings(ctmaFit(ctmaInitFit=init,fit=FALSE,
   backend="ctsemjulia",indVarying="MANIFEST",indVaryingT0=FALSE,
   stationary=TRUE,determinedDiffusions=TRUE,studyCoding=coding,
   activeDirectory=paste0(tempdir(),"/"))))
  p<-out$ctModel$pars;cols<-out$ctModel$determinedDiffusions$study_effect_columns
  expect_length(cols,init$n.studies-1)
  expect_true(all(vapply(cols,function(v)sum(p[[v]])==3,logical(1))))
  expect_true(out$stationary_specification$stationary)
  expect_equal(length(out$data$id),nrow(out$data))
 }
 expect_error(suppressMessages(suppressWarnings(ctmaFit(ctmaInitFit=init,fit=FALSE,
  backend="stan",indVarying="MANIFEST",indVaryingT0=FALSE,stationary=TRUE,
  determinedDiffusions=TRUE,activeDirectory=paste0(tempdir(),"/")))),
  "only for backend='ctsemjulia'")

})

test_that("determined Bayesian study DRIFT mapping preserves covariance contrasts", {
 local_mocked_bindings(.ctmaRequireCtsemJuliaLaplace=function()list(available=TRUE))
 e<-new.env();utils::data("CoTiMAstudyList_3",package="CoTiMA2",envir=e)
 utils::data("CoTiMAInitFit_3",package="CoTiMA2",envir=e)
 fresh<-suppressMessages(suppressWarnings(ctmaInit(primaryStudies=e$CoTiMAstudyList_3,
  n.latent=2,fit=FALSE,silentOverwrite=TRUE,activeDirectory=paste0(tempdir(),"/"))))
 init<-e$CoTiMAInitFit_3;init$emprawList<-fresh$data
 init$studyFitList<-lapply(fresh$data,function(d)list(empraw=d))
 setup<-suppressMessages(suppressWarnings(ctmaFit(ctmaInitFit=init,fit=FALSE,
  backend="ctsemjulia",stationary=TRUE,determinedDiffusions=TRUE,
  optimize=FALSE,priors=TRUE,studyVarying=c("DIFFUSION","T0VAR"),
  backendArgs=list(study_random_effects="DRIFT",intoverpop="none"),
  activeDirectory=paste0(tempdir(),"/"))))
 prep<-suppressMessages(.ctmaPrepareCtsemJuliaLaplace(setup,study_random_effects="DRIFT"))
 p<-prep$model$pars;cols<-setup$ctModel$determinedDiffusions$study_effect_columns
 expect_equal(sum(p$indvarying_study),4)
 expect_true(all(p$matrix[p$indvarying_study]=="PARS"))
 expect_setequal(p$param[p$indvarying_study],setup$ctModel$determinedDiffusions$drift_parameters)
 expect_true(all(vapply(cols,function(v)sum(p[[v]])==3,logical(1))))
 expect_match(setup$stationary_specification$prior_coordinates,"Not equivalent")
 expect_error(.ctmaPrepareCtsemJuliaLaplace(setup,study_random_effects="DIFFUSION"),"derived")
 expect_error(.ctmaPrepareCtsemJuliaLaplace(setup,subject_random_effects="DRIFT"),"study level only")
})
