.ctmaOptimizeRandomDetermined <- function(fit,init,reFits,seed,backend=NULL,
    backendArgs=list(),legacyArgs=NULL,cores=2,iter=5000,finishsamples=NULL,
    activeDirectory=NULL,saveModelFits=FALSE,verbose=1) {
  if (is.null(init)) stop("Supply the ctmaInitFit used for this pooled fit.",call.=FALSE)
  a<-fit$argumentList
  if (!is.list(a) || !isTRUE(a$stationary))
    stop("The baseline must already specify stationary=TRUE; stationarity is never selected randomly.",call.=FALSE)
  if (!isTRUE(a$optimize) || isTRUE(a$priors)) stop("This search supports ML only (optimize=TRUE, priors=FALSE).",call.=FALSE)
  .ctmaOptimizeScalar(reFits,"reFits",integer=TRUE,minimum=1)
  .ctmaOptimizeScalar(seed,"seed",integer=TRUE,minimum=0)
  if(seed > .Machine$integer.max)stop("seed exceeds the supported integer range.",call.=FALSE)
  if(is.null(backend))backend<-a$backend
  if(is.null(backend))backend<-"stan"
  backend<-match.arg(tolower(backend),c("stan","ctsemjulia"))
  if(length(backendArgs) && !is.null(legacyArgs))stop("Supply backendArgs or CoTiMAStanctArgs, not both.",call.=FALSE)
  # Retain the complete recorded recipe, not a manually selected subset.
  args<-a[!duplicated(names(a)) & names(a)%in%names(formals(ctmaFit))]
  args$ctmaInitFit<-init
  args$backend<-backend
  saved<-if(identical(backend,a$backend)) a$backendArgs else NULL
  args$backendArgs<-if(length(backendArgs))utils::modifyList(if(is.null(saved))list()else saved,backendArgs) else if(is.null(saved))list()else saved
  args$CoTiMAStanctArgs<-if(backend=="stan")legacyArgs else NULL
  if(!is.null(legacyArgs))args$backendArgs<-list()
  # Coordinate-dependent starts must never leak from either the recipe or controls.
  args$inits<-NULL;args$backendArgs$inits<-NULL;args$CoTiMAStanctArgs$inits<-NULL
  if(!is.null(args$backendArgs$optimcontrol))args$backendArgs$optimcontrol$init<-NULL
  if(!is.null(args$CoTiMAStanctArgs$optimcontrol))args$CoTiMAStanctArgs$optimcontrol$init<-NULL
  args$stationary<-TRUE;args$optimize<-TRUE;args$priors<-FALSE
  args$coresToUse<-cores;args$iter<-iter;args$verbose<-verbose
  if(!is.null(finishsamples))args$finishsamples<-finishsamples
  if(!is.null(activeDirectory))args$activeDirectory<-activeDirectory
  .ctma_print_message("Random determined-diffusion search retains the stationary model and discards coordinate-dependent raw starts.")
  had_rng<-exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  old_rng<-if(had_rng)get(".Random.seed",envir=.GlobalEnv)else NULL
  on.exit(if(had_rng)assign(".Random.seed",old_rng,envir=.GlobalEnv)else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  # Preparation validates actual generated constraints before consuming RNG or fitting.
  for(v in c(FALSE,TRUE)) {
    check<-args;check$fit<-FALSE;check$determinedDiffusions<-v
    tryCatch(do.call(ctmaFit,check),error=function(e)stop(
      "Cannot search both covariance parametrizations: ",conditionMessage(e),call.=FALSE))
  }
  set.seed(as.integer(seed))
  choices<-sample(c(FALSE,TRUE),reFits,replace=TRUE)
  seeds<-sample.int(.Machine$integer.max,reFits,replace=TRUE)
  attempts<-data.frame(attempt=seq_len(reFits),seed=seeds,stationary=TRUE,
    determinedDiffusions=choices,minus2ll=NA_real_,converged=NA,maximum_gradient=NA_real_,status="failed",error="",stringsAsFactors=FALSE)
  warns<-errs<-vector("list",reFits);best<-NULL;best_value<-Inf
  for(i in seq_len(reFits)) {
    set.seed(seeds[i]);run<-args;run$fit<-TRUE;run$determinedDiffusions<-choices[i]
    message("Fit attempt ",i,"/",reFits,": stationary=TRUE, determinedDiffusions=",choices[i],"; fitting the full model.")
    f<-tryCatch(withCallingHandlers(do.call(ctmaFit,run),warning=function(w){
      warns[[i]]<<-c(warns[[i]],conditionMessage(w));invokeRestart("muffleWarning")}),
      error=function(e){errs[[i]]<<-conditionMessage(e);attempts$error[i]<<-conditionMessage(e);NULL})
    value<-.ctmaOptimizeFitValue(f)
    if(is.finite(value)) {
      attempts$minus2ll[i]<-value;attempts$status[i]<-"completed"
      diagnostics<-tryCatch(.ctmaOptimizeWorkerAssess(.ctmaOptimizeNative(f)),error=function(e)NULL)
      if(!is.null(diagnostics)){
        attempts$converged[i]<-diagnostics$converged
        attempts$maximum_gradient[i]<-diagnostics$maxGradient
      }
      if(value<best_value){best<-f;best_value<-value}
      if(!identical(saveModelFits,FALSE)) {
        if(!is.character(saveModelFits)||length(saveModelFits)!=1L||is.null(args$activeDirectory))
          stop("saveModelFits requires a filename prefix and activeDirectory.",call.=FALSE)
        saveRDS(f,file.path(args$activeDirectory,paste0(saveModelFits,"_",i,".rds")))
      }
    }
  }
  if(is.null(best))stop("All parametrization attempts failed: ",paste(unique(unlist(errs)),collapse="; "),call.=FALSE)
  structure(list(bestFit=best,summary=best$summary,all_minus2ll=attempts$minus2ll,
    all_determinedDiffusions=choices,attempts=attempts,warnings=warns,errors=errs,
    randomDeterminedDiffusions=TRUE,gambling=TRUE,seed=seed,backend=backend,
    note="Finite objective values rank attempts; completion does not certify convergence, identification or uncertainty."),class="CoTiMAFit")
}
