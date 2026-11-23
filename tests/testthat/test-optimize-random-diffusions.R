random_diffusion_fixture <- function() list(argumentList=list(stationary=TRUE,optimize=TRUE,
 priors=FALSE,backend="ctsemjulia",scaleTime=1/12,indVaryingT0=FALSE,
 inits=1:7,backendArgs=list(inits=1:7)))

test_that("new search is opt-in and rejects changed model assumptions", {
 expect_identical(formals(ctmaOptimizeFit)$randomDeterminedDiffusions,FALSE)
 expect_error(ctmaOptimizeFit(gambling=FALSE,randomDeterminedDiffusions=TRUE),"gambling=TRUE")
 expect_error(ctmaOptimizeFit(randomDeterminedDiffusions=TRUE),"pooled")
 f<-random_diffusion_fixture();f$argumentList$stationary<-FALSE
 expect_error(.ctmaOptimizeRandomDetermined(f,list(),2,1),"stationary=TRUE")
 f<-random_diffusion_fixture();f$argumentList$priors<-TRUE
 expect_error(.ctmaOptimizeRandomDetermined(f,list(),2,1),"ML only")
 expect_error(ctmaOptimizeFit(ctmaFitFit=random_diffusion_fixture(),randomDeterminedDiffusions=TRUE,randomRI=TRUE),"only covariance")
})

test_that("choices are reproducible, both preparations checked, starts discarded and failures retained", {
 seen<-list();preflight<-logical()
 mock<-function()NULL
 formals(mock)<-formals(ctmaFit)
 body(mock)<-quote({
  if(!fit){preflight<<-c(preflight,determinedDiffusions);runif(1);return(list())}
  seen[[length(seen)+1L]]<<-list(variant=determinedDiffusions,stationary=stationary,
    inits=inits,backendArgs=backendArgs,scaleTime=scaleTime,indVaryingT0=indVaryingT0)
  list(summary=list(minus2ll=if(determinedDiffusions)10 else 12))
 })
 environment(mock)<-environment()
 local_mocked_bindings(ctmaFit=mock)
 set.seed(123);old<-.Random.seed
 x<-.ctmaOptimizeRandomDetermined(random_diffusion_fixture(),list(),12,42)
 expect_identical(.Random.seed,old)
 expect_identical(preflight,c(FALSE,TRUE))
 expect_setequal(x$all_determinedDiffusions,c(FALSE,TRUE))
 expect_equal(nrow(x$attempts),12)
 expect_true(all(vapply(seen,function(z)isTRUE(z$stationary)&&is.null(z$inits)&&is.null(z$backendArgs$inits)&&identical(z$scaleTime,1/12)&&identical(z$indVaryingT0,FALSE),logical(1))))
 y<-.ctmaOptimizeRandomDetermined(random_diffusion_fixture(),list(),12,42)
 expect_identical(x$attempts,y$attempts)
 expect_equal(x$bestFit$summary$minus2ll,10)
})

test_that("unsupported preparation stops before any fit or RNG change", {
 calls<-0L;mock<-function(...) {calls<<-calls+1L;runif(1);stop("fixed process")}
 local_mocked_bindings(ctmaFit=mock)
 set.seed(2);old<-.Random.seed
 expect_error(.ctmaOptimizeRandomDetermined(random_diffusion_fixture(),list(),2,1),"Cannot search both")
 expect_equal(calls,1L);expect_identical(.Random.seed,old)
})
