safety_prepared <- function() list(
  init_template=list(rawpopmeans=c(0,0),tipredeffectparams=numeric(),studyre_sd=.2,
    studyre_effect=matrix(0,2,1)),standata=list(nstudyre=1L,nstudies=2L),
  parameterization='centered',covariance_blocks=data.frame(block_id=integer(),
    structure=character(),size=integer()))

test_that('automatic RE initializers respect seeds without changing caller RNG', {
  p<-safety_prepared();set.seed(11);before<-.Random.seed
  a<-.ctma_init_function(p,seed=123);b<-.ctma_init_function(p,seed=456)
  expect_identical(a(1),a(1));expect_false(identical(a(1),a(2)))
  expect_false(identical(a(1),b(1)));expect_identical(.Random.seed,before)
  set.seed(11);x<-.ctma_init_function(p)(1);expect_identical(.Random.seed,before)
  set.seed(999);y<-.ctma_init_function(p)(1);expect_false(identical(x,y))
  expect_error(a(NA_real_),'chain_id');expect_error(.ctma_init_function(p,seed=Inf),'seed')
  expect_true(all(is.finite(unlist(.ctma_init_function(p,seed=.Machine$integer.max)(2)))))
  expect_false(identical(.ctma_init_function(p,seed=0)(1),
    .ctma_init_function(p,seed=.Machine$integer.max)(1)))
})

test_that('initializer restores absent RNG state and state on error', {
  old<-if(exists('.Random.seed',.GlobalEnv,inherits=FALSE))get('.Random.seed',.GlobalEnv)else NULL
  on.exit(if(is.null(old)){if(exists('.Random.seed',.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)}else assign('.Random.seed',old,.GlobalEnv),add=TRUE)
  if(exists('.Random.seed',.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)
  .ctma_init_function(safety_prepared(),seed=3)(1)
  expect_false(exists('.Random.seed',.GlobalEnv,inherits=FALSE))
  set.seed(22);before<-.Random.seed;p<-safety_prepared();p$standata$nstudyre<-NA_integer_
  expect_error(.ctma_init_function(p,seed=3)(1));expect_identical(.Random.seed,before)
  expect_identical(.ctma_init_function(list(init_template=NULL),seed=1),'random')
})

test_that('ctmaReFit forwards the effective fitting seed to automatic initialization', {
  seen<-NULL
  local_mocked_bindings(.ctma_init_function=function(prepared,jitter,seed){seen<<-seed;stop('captured initializer')},.package='CoTiMA2')
  p<-list(init='auto',init_jitter=.02)
  expect_error(ctmaReFit(p,seed=42,verbose=FALSE),'captured initializer');expect_equal(seen,42)
  expect_error(ctmaReFit(p,seed=42,fit_args=list(seed=99),verbose=FALSE),'captured initializer');expect_equal(seen,99)
})

test_that('retry execution preserves an existing foreach backend on success and errors', {
  doParallel::registerDoParallel(cores=1)
  on.exit(foreach::registerDoSEQ(),add=TRUE)
  name<-foreach::getDoParName();workers<-foreach::getDoParWorkers()
  expect_equal(.ctmaOptimizeInitApply(1:2,function(i)i*2,1),list(2,4))
  expect_error(.ctmaOptimizeInitApply(1:2,function(i)stop('retry failed'),1),'retry failed')
  expect_error(ctmaOptimizeInit(),'arguments are missing')
  expect_identical(foreach::getDoParName(),name);expect_identical(foreach::getDoParWorkers(),workers)
  if(.Platform$OS.type=='unix'){
    expect_equal(.ctmaOptimizeInitApply(1:2,function(i)i+1,2),list(2,3))
    expect_error(.ctmaOptimizeInitApply(1:2,function(i)stop('owned worker failure'),2),'owned worker failure')
    expect_identical(foreach::getDoParName(),name)
  }
})

test_that('batch saves join paths correctly and never prompt or silently replace files', {
  d<-tempfile();dir.create(d);on.exit(unlink(d,recursive=TRUE),add=TRUE)
  local_mocked_bindings(.ctmaSaveFileInteractive=function()FALSE,
    .ctmaSaveFileRead=function(...)stop('unexpected prompt'),.package='CoTiMA2')
  ctmaSaveFile(d,1,'example.rds','nested')
  path<-file.path(d,'nested','example.rds');expect_equal(readRDS(path),1)
  expect_error(ctmaSaveFile(d,2,'example.rds','nested'),'silentOverwrite = TRUE')
  expect_equal(readRDS(path),1)
  ctmaSaveFile(d,2,'example.rds','nested',silentOverwrite=TRUE)
  expect_equal(readRDS(path),2)
  expect_error(ctmaSaveFile(d,3,'example.rds','nested',silentOverwrite=NA),'TRUE or FALSE')
})

test_that('interactive save choices skip, cancel and recheck renamed targets', {
  d<-tempfile();dir.create(d);on.exit(unlink(d,recursive=TRUE),add=TRUE)
  saveRDS(1,file.path(d,'a.rds'));saveRDS(2,file.path(d,'b.rds'))
  answers<-c('c','b','n')
  local_mocked_bindings(.ctmaSaveFileInteractive=function()TRUE,
    .ctmaSaveFileRead=function(...){out<-answers[1];answers<<-answers[-1];out},.package='CoTiMA2')
  ctmaSaveFile(d,3,'a.rds','');expect_equal(readRDS(file.path(d,'b.rds')),2)
  answers<-'';expect_error(ctmaSaveFile(d,3,'a.rds',''),'cancelled')
  answers<-c('c','new');ctmaSaveFile(d,3,'a.rds','')
  expect_equal(readRDS(file.path(d,'new.rds')),3);expect_equal(readRDS(file.path(d,'a.rds')),1)
})

test_that('data generation retains option restoration after its cluster is registered', {
  old<-getOption('scipen');on.exit(options(scipen=old),add=TRUE)
  on.exit(foreach::registerDoSEQ(),add=TRUE)
  options(scipen=7)
  expect_error(suppressMessages(ctmaGenData(activeDirectory=tempdir(),
    drift=list(matrix(-.5)),diff=list(matrix(1)),T0var=list(matrix(1)),
    lambda=list(matrix(1)),n.latent=1,n.manifest=1,ctmaExtract=TRUE,useRawData=TRUE,envir=NULL,
    coresToUse=1,doPar=FALSE)), 'environment')
  expect_equal(getOption('scipen'),7)
})

test_that('interrupts propagate without replacing the foreach backend', {
  doParallel::registerDoParallel(cores=1);on.exit(foreach::registerDoSEQ(),add=TRUE)
  before<-foreach::getDoParName()
  result<-tryCatch(.ctmaOptimizeInitApply(1,function(i)
    stop(structure(list(message='user interrupt',call=NULL),class=c('interrupt','condition'))),1),
    interrupt=function(e)'interrupted')
  expect_identical(result,'interrupted');expect_identical(foreach::getDoParName(),before)
})


test_that('public retry workflow preserves backend when its fitting callback fails', {
  doParallel::registerDoParallel(cores=1);on.exit(foreach::registerDoSEQ(),add=TRUE)
  before<-foreach::getDoParName();count<-0L
  local_mocked_bindings(ctmaInit=function(...) {
    count<<-count+1L
    list(summary=list(minus2ll=10-count),
      studyFitList=list(list(resultsSummary=list())))
  },.package='CoTiMA2')
  primary<-list(studyNumbers=list(1),deltas=list(1),sampleSizes=list(10),n.studies=1)
  args<-list(primaryStudies=primary,problemStudy=1,reFits=2,n.latent=1,
    activeDirectory=tempdir(),coresToUse=1)
  result<-do.call(ctmaOptimizeInit,args)
  expect_equal(result$summary$minus2ll,8);expect_equal(count,2)
  expect_identical(foreach::getDoParName(),before)
  local_mocked_bindings(ctmaInit=function(...)stop('mocked fit failure'),.package='CoTiMA2')
  expect_error(do.call(ctmaOptimizeInit,args),'mocked fit failure')
  expect_identical(foreach::getDoParName(),before)
})
