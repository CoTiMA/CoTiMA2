test_that("optimizer convergence does not hide identification or saturation", {
  x <- list(optim=list(converged=TRUE,gradient_norm=1e-8,saturated=TRUE,
    saturated_parameters="q2"),identifiability=list(nweak=1L,parameters="q2"),
    estimate=list(raw=c(0,1)))
  d <- ctmaFitDiagnostics(x,print=FALSE)
  expect_true(d$optimizer_converged)
  expect_identical(d$unidentified_directions,1L)
  expect_true(d$saturated)
  output <- capture.output(ctmaFitDiagnostics(x))
  expect_match(paste(output,collapse="\n"),"does not establish identification")
  expect_match(paste(output,collapse="\n"),"Parameters: q2")
})

test_that("missing diagnostics stay unknown and do not trigger Julia", {
  x <- list()
  d <- ctmaFitDiagnostics(x,print=FALSE)
  expect_true(is.na(d$optimizer_converged))
  expect_true(is.na(d$unidentified_directions))
  expect_false(d$scale$assessed)
  expect_output(ctmaFitDiagnostics(x),"not assessed")
})

test_that("binary scale heuristic uses predictor covariance rather than raw Q", {
  make <- function(lambda,v) list(modelbase=list(manifesttype=1),
    ctsemjulia_diagnostics=list(population_matrices=list(
      LAMBDA=matrix(lambda),T0cov=matrix(v))))
  d <- ctmaFitDiagnostics(make(1,121),print=FALSE)
  expect_length(d$scale$flagged,1)
  expect_true(d$scale$assessed)
  # Reciprocal changes of latent units leave the predictor variance invariant.
  expect_identical(d$scale$flagged,ctmaFitDiagnostics(make(.1,12100),print=FALSE)$scale$flagged)
  expect_length(ctmaFitDiagnostics(make(1,4),print=FALSE)$scale$flagged,0)
  expect_output(ctmaFitDiagnostics(make(1,121)),"heuristic")
})

test_that("posterior diagnostics are distinct from placement convergence", {
  x <- list(optim=list(converged=TRUE),sample=list(sampler="nuts",converged=FALSE,
    divergent=3L,rhat=c(a=1.2,b=Inf)))
  d <- ctmaFitDiagnostics(x,print=FALSE)
  expect_true(d$sampling)
  expect_false(d$posterior$reported_converged)
  expect_true(d$posterior$nonfinite_rhat)
  output <- paste(capture.output(ctmaFitDiagnostics(x)),collapse="\n")
  expect_match(output,"Placement optimizer")
  expect_match(output,"divergences=3")
  expect_match(output,"Non-finite R-hat")
})

test_that("study diagnostics remain associated with primary study numbers", {
  object <- list(studyFitList=list(list(optim=list(converged=FALSE)),
    list(optim=list(converged=TRUE))),primaryStudyList=list(studyNumbers=c(4,17)))
  d <- ctmaFitDiagnostics(object,print=FALSE)
  expect_named(d,c("Study 4","Study 17"))
  expect_false(d[[1]]$optimizer_converged)
  expect_true(d[[2]]$optimizer_converged)
})


test_that("bundled NUTS diagnostics use retained raw posterior summaries", {
  x <- list(estimation="bayes",bayes=list(divergences=2L,
    summary=data.frame(variable=c("a","b"),rhat=c(1.01,1.4),ess_bulk=c(100,3))))
  d <- ctmaFitDiagnostics(x,print=FALSE)
  expect_true(d$sampling)
  expect_identical(d$posterior$sampler,"nuts")
  expect_equal(d$posterior$rhat_raw,1.4)
  expect_equal(d$posterior$ess_raw,3)
  expect_equal(d$posterior$divergences,2L)
  expect_match(d$posterior$worst_rhat_parameter,"b")
  expect_true(is.na(d$posterior$reported_converged))
})

test_that("Laplace accuracy remains separate from certified optimization", {
  x <- list(optim=list(converged=TRUE),laplace=list(nrandom=1L,
    correction=list(status="checked",applied=FALSE,gap=.0532,
      check=list(status="unresolved",max_shift=.1006,noise=.0515,bar=.1))))
  d <- ctmaFitDiagnostics(x,print=FALSE)
  expect_true(d$optimizer_converged)
  expect_identical(d$laplace_accuracy$status,"unresolved")
  expect_false(d$laplace_accuracy$correction_applied)
  expect_equal(d$laplace_accuracy$gap,.0532)
  text <- paste(capture.output(ctmaFitDiagnostics(x)),collapse="\n")
  expect_match(text,"Laplace approximation accuracy: unresolved",fixed=TRUE)
  expect_match(text,"does not certify integration accuracy",fixed=TRUE)
  notes <- .ctmaTablesConvergenceNotes(x,list(),"ML")
  expect_match(paste(notes,collapse=" "),"accuracy: unresolved",fixed=TRUE)
  expect_false(any(grepl("accuracy: checked",notes,fixed=TRUE)))
})

test_that("missing and operational Laplace statuses never imply a pass", {
  expect_false(ctmaFitDiagnostics(list(),print=FALSE)$laplace_accuracy$applicable)
  x <- list(laplace=list(nrandom=1L))
  expect_identical(ctmaFitDiagnostics(x,print=FALSE)$laplace_accuracy$status,
    "not assessed in stored diagnostics")
  x$laplace$correction <- list(status="checked")
  expect_match(paste(.ctmaTablesConvergenceNotes(x,list(),"ML"),collapse=" "),
    "no verdict is stored",fixed=TRUE)
  x$laplace$correction$check <- list(status="failed")
  expect_identical(ctmaFitDiagnostics(x,print=FALSE)$laplace_accuracy$status,"failed")
  x$laplace$correction$check$status <- "ok"
  expect_identical(ctmaFitDiagnostics(x,print=FALSE)$laplace_accuracy$status,"ok")
  x$laplace$correction$gap <- c(.1,.2)
  expect_true(is.na(ctmaFitDiagnostics(x,print=FALSE)$laplace_accuracy$gap))
})

test_that("pooled wrappers and primary rosters preserve Laplace status", {
  native <- structure(list(optim=list(converged=TRUE),laplace=list(nrandom=1,
    correction=list(check=list(status="unresolved")))),class="ctJuliaFit")
  pooled <- list(studyFitList=native)
  expect_identical(ctmaFitDiagnostics(pooled,print=FALSE)$laplace_accuracy$status,"unresolved")
  notes <- .ctmaTablesConvergenceNotes(pooled,list(),"ML")
  expect_length(notes,1L)
  roster <- list(studyFitList=list(native,list()),primaryStudyList=list(studyNumbers=c(2,9)))
  d <- ctmaFitDiagnostics(roster,print=FALSE)
  expect_named(d,c("Study 2","Study 9"))
  expect_true(d[[1]]$laplace_accuracy$applicable)
  expect_false(d[[2]]$laplace_accuracy$applicable)
})
