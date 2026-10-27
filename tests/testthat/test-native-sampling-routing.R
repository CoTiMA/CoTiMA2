test_that("native sampling controls retain automatic upstream choices", {
  defaults <- ctmaBackendDefaults("ctsemjulia", "Bayes", "ctmaFit")
  expect_length(.ctmaCtsemJuliaSampleControl(defaults$sampleControl), 0L)
  for (sampler in c("nuts", "saem", "map", "smc")) {
    control <- .ctmaCtsemJuliaSampleControl(defaults$sampleControl,
      list(sampler = sampler))
    expect_identical(control$sampler, sampler)
    expect_null(control$placement)
    expect_null(control$chains)
    expect_null(control$warmup)
  }
  smc <- .ctmaCtsemJuliaSampleControl(defaults$sampleControl,
    list(sampler = "smc", particles = 300, coreset = FALSE))
  expect_identical(smc$particles, 300)
  expect_false(smc$coreset)
  expect_error(.ctmaCtsemJuliaSampleControl(defaults$sampleControl,
    list(sampler = "smc", chains = 4)), "particles")
})

test_that("initial study sampling passes through binary models and explicit targets", {
  captured <- NULL
  local_mocked_bindings(
    .ctmaRequireCtsemJuliaLaplace = function() NULL,
    .ctmaNativeCtsemJuliaFit = function(args) {
      captured <<- args
      list(args = list(resolved = list(intoverpop = args$intoverpop)))
    },
    .ctmaNativeCtsemJuliaSummary = function(...) list(),
    .package = "CoTiMA2"
  )
  model <- list(manifesttype = c(1, 0),
    pars = data.frame(indvarying = c(FALSE, FALSE)))
  defaults <- ctmaBackendDefaults("ctsemjulia", "Bayes", "ctmaInit")
  .ctmaInitCtsemJuliaStudyFit(data.frame(id=1, time=0), model, FALSE,
    defaults, estimation="Bayes", priors=TRUE)
  expect_false(captured$optimize)
  expect_true(captured$priors)
  expect_identical(captured$model$manifesttype, c(1, 0))
  expect_identical(captured$intoverpop, "auto")
  expect_null(captured$sampleControl$sampler)
  model$pars$indvarying[1] <- TRUE
  defaults$intoverpop <- "none"
  defaults$sampleControl <- list(sampler="map")
  result <- .ctmaInitCtsemJuliaStudyFit(data.frame(id=1, time=0), model,
    FALSE, defaults, estimation="Bayes", priors=TRUE)
  expect_identical(captured$sampleControl$sampler, "map")
  expect_null(captured$sampleControl$placement)
  expect_identical(result$ctsemjulia_mode, "joint")
})

test_that("pooled model shapes retain the same explicit integration target", {
  calls <- list()
  local_mocked_bindings(
    .ctmaPrepareCtsemJuliaLaplace = function(object, ...) structure(
      list(model=object$ctModel, data=object$data),
      class="ctmaCtsemJuliaLaplacePrep"),
    .ctmaFitCtsemJuliaLaplace = function(prepared, ...) {
      calls[[length(calls)+1L]] <<- c(list(model=prepared$model), list(...))
      list(args=list(resolved=list(intoverpop="none")))
    },
    .ctmaAnnotateCtsemJuliaResult = function(result, ...) result,
    .package="CoTiMA2"
  )
  for (shape in c("person", "study", "both")) {
    setup <- structure(list(
      data=data.frame(id=1, time=0), argumentList=list(),
      ctModel=list(manifesttype=c(1,0), pars=data.frame(
        param="cint", indvarying=shape != "study"))
    ), class="CoTiMAFit")
    args <- ctmaBackendDefaults("ctsemjulia", "Bayes", "ctmaFit")
    args$intoverpop <- "none"
    args$sampleControl <- list(sampler="map")
    if (shape != "person") args$study_random_effects <- "CINT"
    result <- .ctmaFitCtsemJuliaBackend(setup, args, 2,
      estimation="Bayes", priors=TRUE)
    expect_identical(result$ctsemjulia_mode, "joint")
    expect_identical(tail(calls,1)[[1]]$intoverpop, "none")
    expect_identical(tail(calls,1)[[1]]$model$manifesttype, c(1,0))
  }
})

test_that("sampling targets do not change when a kernel changes", {
  args <- list(intoverpop = "auto", sampleControl = list())
  expect_identical(.ctmaCtsemJuliaRoute(args, "Bayes", FALSE), "auto")
  for (effects in c(TRUE, FALSE)) {
    args$sampleControl <- list(sampler = "nuts")
    expect_identical(.ctmaCtsemJuliaRoute(args, "Bayes", effects), "auto")
  }
  args$sampleControl <- list(sampler = "saem", placement = "saem")
  expect_error(.ctmaCtsemJuliaRoute(args, "Bayes", FALSE), "random effects")
  args$intoverpop <- "laplace"
  expect_error(.ctmaCtsemJuliaRoute(args, "Bayes", TRUE), "joint sampling")
  args$intoverpop <- "none"
  expect_identical(.ctmaCtsemJuliaRoute(args, "Bayes", TRUE), "none")
  args$sampleControl <- list(sampler = "map", placement = "saem")
  expect_error(.ctmaCtsemJuliaRoute(args, "Bayes", TRUE), "placement = 'fit'")
  args$sampleControl <- list(sampler = "smc")
  expect_error(.ctmaCtsemJuliaRoute(args, "Bayes", TRUE), "marginal target")
  args$intoverpop <- "laplace"
  expect_identical(.ctmaCtsemJuliaRoute(args, "Bayes", TRUE), "laplace")
})

test_that("native modes report the target actually resolved by ctsem", {
  result <- list(args = list(resolved = list(intoverpop = "none")))
  expect_identical(.ctmaCtsemJuliaMode(result, "auto", TRUE, TRUE), "joint")
  result$args$resolved$intoverpop <- "laplace"
  expect_identical(.ctmaCtsemJuliaMode(result, "auto", TRUE, TRUE), "laplace")
  expect_identical(.ctmaCtsemJuliaMode(result, "auto", TRUE, FALSE), "direct")
})
