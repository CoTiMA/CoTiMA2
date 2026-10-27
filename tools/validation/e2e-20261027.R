root <- "/home/cdormann/SynologyDrive/R/CoTiMA2-Benchmark/e2e-20261027"
dir.create(root, recursive=TRUE, showWarnings=FALSE)
.libPaths(c("/home/cdormann/R/ctsem-validation-20261027",
  path.expand("~/R/ctsem-juliaFit-lib"), .libPaths()))
Sys.unsetenv(c("JULIA_HOME", "JULIA_BINDIR"))
library(CoTiMA2)
ctsem::ctJuliaSetup(threads=4L, agree=TRUE)
set.seed(20261027)
for (s in 1:3) {
  n <- 64; waves <- 6
  x <- y <- matrix(0,n,waves)
  x[,1] <- rnorm(n); y[,1] <- rnorm(n)
  for (j in 2:waves) {
    x[,j] <- .6*x[,j-1]+rnorm(n,0,.8)
    y[,j] <- .7*y[,j-1]+rnorm(n,0,.6)
  }
  binary <- matrix(rbinom(n*waves,1,plogis(y-.5)),n,waves)
  wide <- matrix(0,n,2*waves+waves-1)
  wide[,seq(1,2*waves,2)] <- x
  wide[,seq(2,2*waves,2)] <- binary
  wide[,(2*waves+1):ncol(wide)] <- 1
  path <- file.path(root,paste0("study",s,".txt"))
  write.table(wide,path,row.names=FALSE,col.names=FALSE)
  assign(paste0("rawData",s),list(fileName=path,studyNumbers=s,
    missingValues=-99,standardize=FALSE,header=FALSE,dec=".",sep=" "))
  assign(paste0("delta_t",s),rep(NA_real_,waves-1))
  assign(paste0("sampleSize",s),n)
}
studies <- ctmaPrep(selectedStudies=1:3,activeDirectory=paste0(root,"/"),
  n.latent=2,summary=FALSE)
saveRDS(studies,file.path(root,"studies.rds"))
controls <- list(optimcontrol=list(maxiter=150L,uncertainty="hessian",
  finishsamples=100L))
results <- list()
run <- function(label,expr) {
  writeLines(label,file.path(root,"CURRENT"))
  warnings <- character(); start <- proc.time()[[3]]
  answer <- tryCatch(withCallingHandlers(force(expr),
    warning=function(w) warnings <<- c(warnings,conditionMessage(w))),
    error=function(e) structure(list(error=conditionMessage(e)),class="workflow_error"))
  saveRDS(answer,file.path(root,paste0(label,".rds")))
  results[[label]] <<- list(status=if(inherits(answer,"workflow_error")) "error" else "completed",
    seconds=proc.time()[[3]]-start,warnings=warnings,error=answer$error)
  saveRDS(results,file.path(root,"results.rds"))
  if(inherits(answer,"workflow_error")) stop(label,": ",answer$error)
  answer
}
initial <- run("init-ML",ctmaInit(primaryStudies=studies,
  activeDirectory=paste0(root,"/"),silentOverwrite=TRUE,
  backend="ctsemjulia",binaries=c(0,1),n.latent=2,coresToUse=4L,
  optimize=TRUE,priors=FALSE,manifestMeans="auto",
  drift=c("V1toV1",0,0,"V2toV2"),
  diff=c("V1withV1",0,0,.3),backendArgs=controls))
fit <- run("fit-ML",ctmaFit(ctmaInitFit=initial,backend="ctsemjulia",
  binaries=c(0,1),coresToUse=4L,optimize=TRUE,priors=FALSE,
  manifestMeans="auto",backendArgs=controls))
run("summary-ML",summary(fit))
sampling <- controls
sampling$sampleControl <- list(sampler="nuts",placement="fit",chains=2L,
  warmup=50L,draws=50L,maxDraws=100L,minESS=10,rhatTarget=1.2,
  processes=FALSE,seed=20261027L)
bayes <- run("fit-Bayes",ctmaFit(ctmaInitFit=initial,backend="ctsemjulia",
  binaries=c(0,1),coresToUse=4L,optimize=FALSE,priors=TRUE,
  manifestMeans="auto",backendArgs=sampling))
run("summary-Bayes",summary(bayes))
writeLines(capture.output(sessionInfo()),file.path(root,"sessionInfo.txt"))
writeLines(as.character(Sys.time()),file.path(root,"WORKFLOW_FINISHED"))
