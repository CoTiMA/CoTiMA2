test_that('public core resolution respects constrained availability and validates inputs', {
 local_mocked_bindings(.ctmaAvailableCores=function(...) 2L)
 expect_identical(.ctmaResolveCores(63L, announce=FALSE)$resolved, 2L)
 expect_identical(.ctmaResolveCores(-1L, announce=FALSE)$resolved, 1L)
 for (x in list(NULL, NA_real_, Inf, 1.5, c(1,2), '2'))
  expect_error(.ctmaResolveCores(x,announce=FALSE), 'one integer')
 expect_error(.ctmaResolveJuliaThreads(3L), 'budget')
 expect_error(.ctmaCheckJuliaBudget(2L, 2L, 2L, 1L), 'budget')
 expect_error(.ctmaCheckJuliaBudget(2L, 1L, 2L, 2L), 'budget')
 expect_equal(.ctmaCheckJuliaBudget(1L,2L,1L,1L),2)
})

test_that('Julia allocation uses the resolved budget even when machine is larger', {
 defaults <- ctmaBackendDefaults('julia','ML','ctmaFit')
 a <- .ctmaInitJuliaResources(2L,'ML',list(),defaults,announce=FALSE)$backend_args
 expect_identical(a$julia_threads,2L)
 expect_identical(a$subject_threads,2L)
 expect_error(.ctmaInitJuliaResources(2L,'ML',list(julia_threads=3L),
   utils::modifyList(defaults,list(julia_threads=3L)),announce=FALSE),'budget')
 defaults <- ctmaBackendDefaults('julia','Bayes','ctmaFit')
 defaults$chains <- 1L
 a <- .ctmaInitJuliaResources(2L,'Bayes',list(),defaults,announce=FALSE)$backend_args
 expect_identical(a$parallel_chains,1L)
 expect_identical(a$julia_threads,2L)
 expect_identical(a$threads_per_chain,2L)
 defaults$chains <- 4L
 a <- .ctmaInitJuliaResources(2L,'Bayes',list(),defaults,announce=FALSE)$backend_args
 expect_identical(a$parallel_chains,2L)
 expect_identical(a$threads_per_chain,1L)
})

test_that('owned subprocess limits restore caller environment including absent variables', {
 withr::local_envvar(c(OMP_NUM_THREADS='7', OPENBLAS_NUM_THREADS=NA))
 before <- Sys.getenv(c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS'),unset=NA)
 f <- function(fail) {restore <- .ctmaLimitThreadEnvironment();on.exit(restore())
   expect_equal(Sys.getenv('OMP_NUM_THREADS'),'1')
   expect_equal(Sys.getenv('OPENBLAS_NUM_THREADS'),'1')
   if(fail)stop('test exit')}
 f(FALSE);expect_identical(Sys.getenv(names(before),unset=NA),before)
 expect_error(f(TRUE),'test exit')
 expect_identical(Sys.getenv(names(before),unset=NA),before)
})

test_that('managed project copies are content-keyed and preserve installed templates', {
 root <- tempfile();dir.create(root);on.exit(unlink(root,recursive=TRUE))
 src <- file.path(root,'installed');dir.create(src)
 for(f in c('Project.toml','Manifest.toml','CtsemJulia.jl')) writeLines(f,file.path(src,f))
 before <- tools::md5sum(list.files(src,full.names=TRUE))
 Sys.chmod(names(before), mode='0444')
 cache <- file.path(root,'user-cache')
 plan <- .ctmaJuliaManagedProject(FALSE,src,cache)
 expect_false(dir.exists(cache))
 target <- .ctmaJuliaManagedProject(TRUE,src,cache)
 expect_identical(normalizePath(plan),target)
 expect_true(all(file.access(file.path(target,c('Project.toml','Manifest.toml')),2L)==0L))
 expect_true(all(file.exists(file.path(target,basename(names(before))))))
 writeLines('resolved manifest',file.path(target,'Manifest.toml'))
 expect_identical(.ctmaJuliaManagedProject(TRUE,src,cache),target)
 expect_identical(readLines(file.path(target,'Manifest.toml')),'resolved manifest')
 expect_identical(tools::md5sum(names(before)),before)
 Sys.chmod(names(before), mode='0644')
 writeLines('new template',file.path(src,'Project.toml'))
 expect_false(identical(.ctmaJuliaManagedProject(FALSE,src,cache),plan))
 writeLines('modified source',file.path(target,'CtsemJulia.jl'))
 writeLines('Project.toml',file.path(src,'Project.toml'))
 expect_error(.ctmaJuliaManagedProject(TRUE,src,cache),'source was changed')
})

test_that('noninteractive setup requires consent before project writes or subprocesses', {
 created <- FALSE
 local_mocked_bindings(.ctmaJuliaManagedProject=function(create=FALSE,...) {
   if(create) created <<- TRUE
   tempfile('planned-project-')
 })
 expect_error(ctmaJuliaSetup(julia=file.path(R.home('bin'),'R'),
   verbose=FALSE), 'explicit consent')
 expect_false(created)
 plan <- ctmaJuliaSetup(julia=file.path(R.home('bin'),'R'),dry_run=TRUE,verbose=FALSE)
 expect_false(created)
 expect_s3_class(plan,'ctmaJuliaSetup')
})


test_that('nested native and legacy Stan core controls cannot exceed public budget', {
 local_mocked_bindings(.ctmaAvailableCores=function(...) 2L)
 x <- .ctmaBoundCtFitResources(list(cores=1L,optimcontrol=list(cores=63L)))
 expect_identical(x$cores,1L)
 expect_identical(x$optimcontrol$cores,1L)
 x <- .ctmaBoundCtFitResources(list(cores=63L,optimcontrol=list(cores=63L)))
 expect_identical(x$cores,2L)
 expect_identical(x$optimcontrol$cores,2L)
})


test_that('configuration pilots cannot reset a resolved smaller pool', {
 seen <- integer()
 local_mocked_bindings(.ctmaAvailableCores=function(...) 4L,
   ctmaFit=function(ctmaInitFit=NULL, fit=TRUE, backend=NULL, backendArgs=list(),
     coresToUse=2L, ...) list(),
   ctmaJuliaFit=function(object, ...) {
     seen <<- c(seen,list(...)$julia_threads)
     stop('captured pilot boundary')
   })
 expect_error(ctmaFitConfig(list(),coresToUse=2L,verbose=FALSE), 'captured pilot boundary')
 expect_error(ctmaFitConfig(list(),coresToUse=-1L,verbose=FALSE), 'captured pilot boundary')
 expect_identical(seen,c(2L,3L))
})


test_that('coreset is never partially matched as a core count', {
 local_mocked_bindings(.ctmaAvailableCores=function(...) 2L)
 for(value in c(FALSE,TRUE)) {
   control <- list(coreset=value)
   x <- .ctmaBoundCtFitResources(list(cores=2L,optimcontrol=control))
   expect_identical(x$optimcontrol,control)
   expect_identical(x$cores,2L)
 }
 x <- .ctmaBoundCtFitResources(list(cores=2L,optimcontrol=list(coreset=FALSE,cores=63L)))
 expect_identical(x$optimcontrol$cores,2L)
 expect_false(x$optimcontrol$coreset)
 expect_identical(.ctmaBoundCtFitResources(list(coreset=FALSE)),list(coreset=FALSE))
})
