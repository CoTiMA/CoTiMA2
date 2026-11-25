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
