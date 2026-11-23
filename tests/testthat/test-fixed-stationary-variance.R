fixed_variance_model <- function() {
  m <- ctsem::ctModel(type='ct', LAMBDA=matrix(1), DRIFT=matrix('a'),
    DIFFUSION=matrix('q'), T0VAR=matrix(1), T0MEANS=matrix(0),
    MANIFESTMEANS=matrix(0), MANIFESTVAR=matrix(0), CINT=matrix(0),
    TIpredNames='TI1', n.TIpred=1, silent=TRUE)
  m$pars$indvarying <- FALSE
  m$pars$TI1_effect <- m$pars$matrix=='DRIFT'
  m
}
test_that('fixed stationary variance keeps study drift and derives diffusion', {
  m <- fixed_variance_model()
  a <- CoTiMA2:::.ctmaDeterminedDiffusionModel(m, TRUE, stationaryVariance=1)
  p <- a$pars
  expect_equal(sum(!is.na(p$param) & p$matrix=='PARS'),1)
  expect_true(p$TI1_effect[p$matrix=='PARS'])
  expect_false(any(p$TI1_effect[p$matrix!='PARS']))
  expect_identical(p$transform[p$matrix=='PARS'],m$pars$transform[m$pars$matrix=='DRIFT'])
  for (drift in c(-.01,-.5,-3)) {
    env <- list(a=drift)
    q <- eval(parse(text=p$param[p$matrix=='DIFFUSION']),env)^2
    expect_equal(q/(-2*drift),1)
  }
  expect_equal(p$value[p$matrix=='T0VAR'],1)
})
test_that('fixed variance rejects incompatible settings without silent changes', {
  expect_error(CoTiMA2:::.ctmaStationaryOptions(TRUE,TRUE,backend='stan',stationaryVariance=1),'ctsemjulia')
  expect_error(CoTiMA2:::.ctmaStationaryOptions(TRUE,FALSE,backend='ctsemjulia',stationaryVariance=1),'determinedDiffusions')
  expect_error(CoTiMA2:::.ctmaStationaryOptions(TRUE,TRUE,backend='ctsemjulia',stationaryVariance=0),'positive')
  m<-fixed_variance_model();m$pars$TI1_effect[m$pars$matrix=='DIFFUSION']<-TRUE
  expect_error(CoTiMA2:::.ctmaDeterminedDiffusionModel(m,TRUE,stationaryVariance=1),'separate study effects')
  m<-fixed_variance_model();m$pars$value[m$pars$matrix=='T0VAR']<-2
  expect_error(CoTiMA2:::.ctmaDeterminedDiffusionModel(m,TRUE,stationaryVariance=1),'inconsistent')
  expect_identical(CoTiMA2:::.ctmaDeterminedDiffusionModel(m,FALSE),m)
})
test_that('native mode covariance retains all contrast covariance terms and time scaling', {
  raw<-c(-.3,.2);v<-matrix(c(.04,.01,.01,.09),2)
  X<-cbind(1,c(1,-2))
  fake<-structure(list(args=list(optimize=TRUE),estimate=list(raw=raw,cov=v),
    model_spec=list(tipred_data=matrix(c(1,-2),ncol=1,dimnames=list(NULL,'TI1'))),
    argumentList=list(scaleTime=2)),class='ctJuliaFit')
  local_mocked_bindings(.ctmaNativeStudyExtract=function(fit,profiles,cores,raw=NULL){
    vals <- -exp(raw%*%t(X));list(subj_DRIFT=array(vals,c(nrow(raw),2,1,1)))
  },.package='CoTiMA2')
  out<-ctmaStudyEstimates(fake)
  mu<--exp(X%*%raw);J<-diag(as.vector(mu))%*%X*2
  expect_equal(out$studyEstimates$estimate,as.vector(mu)*2,tolerance=1e-9)
  expect_equal(unname(out$covariance),unname(J%*%v%*%t(J)),tolerance=1e-8)
  expect_true(abs(out$covariance[1,2])>1e-3)
  expect_identical(out$covariance, t(out$covariance))
  fake$args <- list(resolved=list(optimize=TRUE))
  expect_equal(ctmaStudyEstimates(fake)$covariance, out$covariance)
  expect_error(ctmaStudyEstimates(fake,derivativeStep=0),'positive')
})
