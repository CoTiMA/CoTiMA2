.libPaths(c("/home/cdormann/R/ctsem-validation-20261027",
  path.expand("~/R/ctsem-juliaFit-lib"),.libPaths()))
Sys.unsetenv(c("JULIA_HOME","JULIA_BINDIR"))
library(CoTiMA2)
ctsem::ctJuliaSetup(threads=4L,agree=TRUE)
root <- "/home/cdormann/SynologyDrive/R/CoTiMA2-Benchmark/e2e-20261027"
initial <- readRDS(file.path(root,"init-ML.rds"))
warnings <- character()
timing <- system.time(fit <- withCallingHandlers(ctmaFit(ctmaInitFit=initial,
  backend="ctsemjulia",binaries=c(0,1),coresToUse=4L,optimize=FALSE,
  priors=TRUE,manifestMeans="auto",backendArgs=list(
    optimcontrol=list(maxiter=150L,uncertainty="hessian",finishsamples=100L),
    sampleControl=list(sampler="nuts",placement="fit",chains=2L,
      warmup=250L,draws=100L,maxDraws=200L,minESS=20,rhatTarget=1.1,
      target_accept=.95,processes=FALSE,seed=20261028L))),
  warning=function(w) warnings <<- c(warnings,conditionMessage(w))))
saveRDS(fit,file.path(root,"fit-Bayes-refined.rds"))
saveRDS(list(timing=timing,warnings=warnings),file.path(root,"refined-result.rds"))
summary(fit)
