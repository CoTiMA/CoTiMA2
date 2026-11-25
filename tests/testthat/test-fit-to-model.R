test_that("model arguments retain explicit NULL and FALSE", {
  model <- list(argumentList=list(scaleTime=.1, indVarying="MANIFEST"))
  args <- .ctmaModelCallArguments(model,
    quote(ctmaInit(model=model, scaleTime=NULL, indVarying=FALSE)),environment(),ctmaInit)
  expect_true("scaleTime" %in% names(args))
  expect_null(args$scaleTime)
  expect_identical(args$indVarying,FALSE)
})

test_that("invalid conversion targets fail before reading study data", {
  expect_error(ctmaFitToModel(list(),"other"),"Unknown model target")
  expect_error(ctmaFitToModel(list(),NA_character_),"model must be")
})

test_that("original pooled models do not require redundant wide caches", {
  f <- list(studyList=list(list(originalStudyNo=1,empcov=diag(2),delta_t=12,sampleSize=20)),
    n.latent=1L, argumentList=list(scaleTime=1/12),
    data=data.frame(id=c(1,1),time=c(0,1),V1=c(.2,.3)),
    studyFitList=structure(list(modelbase=list(timeName="time",pars=data.frame(
      matrix="DRIFT",row=1,col=1,param="d",value=NA_real_))),class="ctJuliaFit"))
  m <- ctmaFitToModel(f)
  expect_identical(m$data$time,c(0,12))
  expect_identical(m$data$V1,f$data$V1)
  expect_true(m$preserve)
})
