point_table <- function() data.frame(matrix="DRIFT",row=1,col=1,
  Mean=-.8,sd=.2,`2.5%`=-1.2,`50%`=-.75,`97.5%`=-.4,mode=-.7,check.names=FALSE)

test_that("mode ordering preserves all draw summaries and diagnostics", {
  x<-point_table()
  answer<-.ctmaPointSummary(list(parmatrices=x),list(),optimize=TRUE,priors=FALSE)
  expect_identical(names(answer$parmatrices)[4:6],c("mode","Mean","sd"))
  expect_identical(answer$parmatrices$Mean,x$Mean)
  expect_identical(answer$parmatrices$sd,x$sd)
  expect_identical(answer$pointEstimate$selected,"mode")
  expect_match(answer$pointEstimateNote,"ML point estimates",fixed=TRUE)
})

test_that("Bayesian reporting uses draws and never invents a mode", {
  x<-point_table();x$mode<-NULL
  a<-.ctmaPointSummary(list(parmatrices=x),list(sample=list(chains=4)))
  expect_identical(a$pointEstimate$selected,"mean")
  expect_false("mode"%in%names(a$parmatrices))
  expect_match(a$pointEstimateNote,"posterior draws",fixed=TRUE)
})

test_that("legacy mean ratios remain explicit and no mode Wald is constructed", {
  x<-point_table()
  expect_identical(.ctmaPointColumn(x),"mode")
  expect_identical(.ctmaPointColumn(x,"mean"),"Mean")
  expect_true(is.na(.ctmaPointQuotient(x)))
  expect_equal(.ctmaPointQuotient(x,"mean"),-4)
  expect_error(.ctmaPointColumn(x[,-9],"mode"),"No optimized")
})

test_that("old objects receive a package-name-independent migration notice", {
  old<-structure(list(argumentList=list(optimize=TRUE,priors=FALSE)),class="CoTiMAFit")
  expect_match(.ctmaLegacyPointNotice(old),"pointEstimate = 'mean'",fixed=TRUE)
  old$argumentList$pointEstimate<-"auto"
  expect_null(.ctmaLegacyPointNotice(old))
})

test_that("effective backend priors take priority over caller fallbacks", {
  f<-list(args=list(resolved=list(optimize=TRUE,priors=TRUE)))
  expect_identical(.ctmaPointMethod(f,optimize=TRUE,priors=FALSE),"MAP")
})

test_that("mixed ML and Bayesian tables label the shared point estimate honestly", {
  a<-structure(list(summary=list(parmatrices=point_table()),argumentList=list(
    optimize=TRUE,priors=FALSE,pointEstimate="auto")),class="CoTiMAFit")
  b<-a;b$summary$parmatrices$mode<-NULL;b$argumentList$optimize<-FALSE
  out<-ctmaTables(list(a,b),parameters="drift",statistics="estimate",
    modelNames=c("ML","Bayes"),format="none",interactive=FALSE)
  expect_true("estimate"%in%names(out$data))
  expect_equal(out$data$estimate,c(-.7,-.8))
})

test_that("moderator profiles evaluate ML points and retain draw intervals", {
  testthat::local_mocked_bindings(.ctBackendPopArrays=function(fit,samples,tipreds) {
    list(pop_DRIFT=array(-exp(samples[,1]+tipreds[1]),c(nrow(samples),1,1)))
  },.package="ctsem")
  model<-list(n.latent=1L,TIpredNames="TI1",pars=data.frame(matrix="DRIFT",row=1,col=1,param="d"))
  fit<-structure(list(estimate=list(raw=log(.4),rawposterior=matrix(log(c(.5,.7,1)),ncol=1)),
    model_spec=list(model=model),argumentList=list(optimize=TRUE,priors=FALSE,pointEstimate="auto")),class="ctJuliaFit")
  profiles<-matrix(c(0,1),ncol=1,dimnames=list(c("A","B"),"TI1"))
  a<-ctmaTables(fit,parameters="drift",moderatorProfiles=profiles,format="none")
  expect_equal(a$data$estimate,c(-.4,-.4*exp(1)))
  b<-ctmaTables(fit,parameters="drift",moderatorProfiles=profiles,statistics="mean",format="none")
  expect_equal(b$data$estimate[1],-mean(c(.5,.7,1)))
  expect_equal(a$data$lower,b$data$lower)
})
