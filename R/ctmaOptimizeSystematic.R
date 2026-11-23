# Systematic ML pilots/refits. Each case is an owned, budgeted R subprocess.
.ctmaOptimizeControls <- function(x) {
  defaults<-list(maxiter=100L,maxSeconds=60,fullMaxSeconds=600,maxFullChecks=0L,jitter=.01)
  if(!is.list(x)||any(!names(x)%in%names(defaults))|| (length(x)&&(is.null(names(x))||anyDuplicated(names(x))))) .ctmaOptimizeStop('Unknown or unnamed quickControl setting.')
  x<-utils::modifyList(defaults,x)
  for(n in c('maxiter','maxFullChecks')) .ctmaOptimizeScalar(x[[n]],paste0('quickControl$',n),TRUE,if(n=='maxiter')1 else 0)
  for(n in c('maxSeconds','fullMaxSeconds')) .ctmaOptimizeScalar(x[[n]],paste0('quickControl$',n),FALSE,1)
  .ctmaOptimizeScalar(x$jitter,'quickControl$jitter',FALSE,0)
  x
}
.ctmaOptimizeSetIterations <- function(control,backend,maxiter,optimizer=NULL) {
  if(backend=='stan') {
    if(is.null(optimizer))optimizer<-tryCatch(getFromNamespace('stanoptimis','ctsem'),error=function(e)NULL)
    if(is.function(optimizer)&&!any(c('maxiter','...')%in%names(formals(optimizer)))) {
      control$maxiter<-NULL
      warning('This ctsem Stan optimizer exposes no iteration cap; the subprocess wall-time cap remains active.',call.=FALSE)
      return(control)
    }
  }
  control$maxiter<-as.integer(maxiter)
  control
}
.ctmaOptimizePreparedLayout <- function(prepared,backend) {
  if(backend=='ctsemjulia') {
    s<-unclass(prepared)
    used<-c(s$parameter_table$parnumber,s$ti_effects$coefficient,
      s$laplace$sd_index,s$laplace$cor_index)
    n<-max(c(0,used),na.rm=TRUE)
    fake<-structure(list(model_spec=s,estimate=list(raw=numeric(n))),class='ctJuliaFit')
  } else {
    d<-prepared$standata
    n<-d$nparams+d$nindvarying+d$nindvaryingoffdiagonals+max(c(0,d$TIPREDEFFECTsetup))
    fake<-prepared;fake$stanfit<-list(rawest=numeric(n))
    class(fake)<-'ctStanFit'
  }
  .ctmaOptimizeLayout(fake)
}
.ctmaOptimizeWorkerAssess <- function(native) {
  if(inherits(native,'ctJuliaFit')) {
    gradient<-native$optim$gradient_norm
    converged<-isTRUE(native$optim$converged)
    certification<-native$uncertainty$certification
    stationary<-if(is.null(certification$certified))NA else isTRUE(certification$certified)
    weak<-native$identifiability$nweak
    if(is.null(weak))weak<-NA_integer_
    repair<-attr(native$estimate$cov,'ctOptimCovFromHessian')
    bad<-isTRUE(native$optim$stopped_by_divergence)||isTRUE(native$optim$overshot)||isTRUE(native$optim$stalled)
    ll<-native$estimate$loglik
  } else {
    # Convergence codes differ across ctsem releases. Preserve unknown status;
    # never infer success merely from a finite likelihood or small LP change.
    op<-native$stanfit$optimfit
    converged<-if(is.logical(op$converged)&&length(op$converged)==1L)op$converged else NA
    gradient<-if(is.numeric(op$gradient))max(abs(op$gradient)) else NA_real_
    stationary<-NA;weak<-NA_integer_;bad<-FALSE
    repair<-attr(native$stanfit$cov,'ctOptimCovFromHessian')
    ll<-native$stanfit$transformedparsfull$ll
    if(is.null(ll))ll<-tryCatch(as.numeric(stats::logLik(native)),error=function(e)NA_real_)
  }
  hessian<-if(inherits(native,'ctJuliaFit'))native$uncertainty$hessian else native$stanfit$hessian
  raw<-if(inherits(native,'ctJuliaFit'))native$estimate$raw else native$stanfit$rawest
  curvatureAssessed<-is.matrix(hessian)&&is.numeric(hessian)&&length(raw)>0L&&
    identical(dim(hessian),c(length(raw),length(raw)))&&all(is.finite(hessian))
  if(is.null(gradient))gradient<-NA_real_
  if(length(ll)!=1L||!is.finite(ll))ll<-NA_real_
  repaired<-if(is.null(repair))NA else isTRUE(repair$usedNullProjection)||isTRUE(repair$usedNearPD)||isTRUE(repair$infoNearPD)
  if(!is.null(repair$nullDirections))weak<-max(c(weak,repair$nullDirections),na.rm=TRUE)
  # The engine initializes nweak=0 even when estonly skipped curvature.
  # Absence of evidence must not be displayed as a successful zero count.
  if(!curvatureAssessed)weak<-NA_integer_
  list(loglik=as.numeric(ll),converged=converged,maxGradient=as.numeric(gradient),
    curvatureAssessed=curvatureAssessed,
    stationary=stationary,weakDirections=weak,repaired=repaired,bad=bad,
    uncertaintyValidated=isTRUE(curvatureAssessed)&&isTRUE(converged)&&isTRUE(stationary)&&isFALSE(repaired)&&identical(as.integer(weak),0L))
}
.ctmaOptimizeCaseWorker <- function(payload) {
  Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',
    JULIA_NUM_THREADS=as.character(payload$cores),JULIA_NUM_PRECOMPILE_TASKS=as.character(payload$cores))
  Sys.unsetenv('JULIACONNECTOR_SERVER')
  saveRDS(list(pid=Sys.getpid(),tempdir=tempdir(),started=as.numeric(Sys.time())),payload$ownerFile)
  set.seed(payload$seed)
  warnings<-character()
  run<-function() {
    backend<-payload$case$backend
    if(backend=='ctsemjulia') {
      .ctmaRequireCtsemJuliaLaplace()
      on.exit({try(getExportedValue("ctsem", "ctJuliaKill")('workers',ask=FALSE),silent=TRUE);try(getExportedValue("ctsem", "ctJuliaKill")('session',ask=FALSE),silent=TRUE)},add=TRUE)
      getExportedValue("ctsem", "ctJuliaSetup")(threads=as.integer(payload$cores),agree=FALSE)
      # Record only this fresh subprocess's engine PID, for safe timeout cleanup.
      pid<-tryCatch(as.integer(getFromNamespace('.ctJuliaEval','ctsem')('getpid()')),error=function(e)NA_integer_)
      if(is.finite(pid)&&requireNamespace('ps',quietly=TRUE)) {
        h<-ps::ps_handle(pid)
        saveRDS(list(pid=pid,created=ps::ps_create_time(h)),payload$engineFile)
      }
    }
    allowed<-c('intoverstates','binomial','intoverpop','sameInitialTimes','stationary',
      'savescores','savesubjectmatrices','saveComplexPars','nlcontrol')
    args<-payload$case$resolved[intersect(allowed,names(payload$case$resolved))]
    if(is.null(args$intoverpop))args$intoverpop<-'auto'
    args<-utils::modifyList(args,payload$backendArgs)
    # No backend override may alter target, model, data or priors in this search.
    forbidden<-intersect(names(payload$backendArgs),c('model','datalong','backend','optimize','priors','fit','inits','intoverpop','stationary','intoverstates','sampleControl'))
    if(length(forbidden)) .ctmaOptimizeStop('Systematic backendArgs cannot change the saved model/target: ',paste(forbidden,collapse=', '))
    args$datalong<-payload$case$data;args$model<-payload$case$model
    args$backend<-if(backend=='ctsemjulia')'julia' else 'stan'
    # Published ctsem 3.x has no backend formal; its default engine is Stan.
    if(backend=='stan'&&!('backend'%in%names(formals(ctsem::ctFit))))args$backend<-NULL
    args$optimize<-TRUE;args$priors<-FALSE;args$cores<-payload$cores
    args$verbose<-payload$verbose;args$fit<-FALSE
    prepared<-do.call(ctsem::ctFit,args)
    layout<-.ctmaOptimizePreparedLayout(prepared,backend)
    keys<-layout$keys
    if(any(!keys%in%names(payload$case$targets))) .ctmaOptimizeStop('A new free coordinate has no scientifically mapped starting value.')
    raw<-vapply(keys,function(k)payload$case$targets[[k]],numeric(1))
    if(!is.null(payload$point)) {
      if(!identical(names(payload$point),keys)) .ctmaOptimizeStop('Pilot-to-full parameter layout changed.')
      raw<-as.numeric(payload$point)
    }
    startLL<-NA_real_
    if(backend=='ctsemjulia') {
      v<-getExportedValue("ctsem", "ctJuliaEvaluate")(prepared,pars=raw,gradient=FALSE,contributions=FALSE,cores=payload$cores)
      startLL<-as.numeric(v$value)
      if(payload$case$case$fullData && is.finite(payload$expectedLL) &&
          is.null(payload$point)&&abs(startLL-payload$expectedLL)>1e-5)
        .ctmaOptimizeStop('Full-data reference/start transformation failed likelihood invariance; no optimization was run.')
    }
    if(payload$attempt>1L)raw<-raw+stats::rnorm(length(raw),sd=payload$jitter)
    args$fit<-TRUE;args$inits<-raw
    control<-args$optimcontrol
    if(is.null(control))control<-list()
    control<-.ctmaOptimizeSetIterations(control,backend,payload$maxiter)
    control$estonly<-isTRUE(payload$pilot)
    control$finishsamples<-if(payload$pilot)0L else payload$finishsamples
    if(payload$pilot&&backend=='stan'&&!('backend'%in%names(formals(ctsem::ctFit)))) {
      # ctsem 3.x's estonly return omits standata that its TI postprocessor
      # subsequently requires. Use its complete return path, still budgeted.
      control$estonly<-FALSE
      control$finishsamples<-20L
      warning('This published ctsem Stan API needs its standard uncertainty postprocessing for a complete TI fit; the pilot retains only 20 approximate draws and remains inconclusive.',call.=FALSE)
    }
    if(backend=='ctsemjulia') {
      control$batch<-FALSE;control$coreset<-FALSE;control$carefulfit<-FALSE
      if(!payload$pilot)control$uncertainty<-'hessian'
    }
    args$optimcontrol<-control
    native<-do.call(ctsem::ctFit,args)
    after<-.ctmaOptimizeLayout(native)
    assessment<-.ctmaOptimizeWorkerAssess(native)
    if(backend=='ctsemjulia') {
      v<-getExportedValue("ctsem", "ctJuliaEvaluate")(native,pars=after$raw,gradient=TRUE,contributions=FALSE,cores=payload$cores)
      assessment$loglik<-as.numeric(v$value);assessment$maxGradient<-max(abs(v$gradient))
    }
    setup<-list(data=payload$case$data,ctModel=payload$case$model,
      ctModelText=.ctmaModelText(payload$case$model),studyList=payload$case$studies)
    wrapped<-.ctmaWrapPooledJuliaFit(native,setup)
    if(backend=='stan')class(wrapped)<-'CoTiMAFit'
    wrapped$argumentList<-payload$argumentList
    wrapped$argumentList$backend<-backend
    wrapped$argumentList$scaleTime<-payload$case$scaleTime
    wrapped$summary<-list(loglik=assessment$loglik,minus2ll=-2*assessment$loglik,n.parameters=length(after$raw))
    wrapped$optimizer_controls<-control
    list(fit=wrapped,point=after$raw,assessment=assessment,startLoglik=startLL,
      pilot=payload$pilot,case=payload$case$case,
      parameterMap=data.frame(key=keys,start=raw,stringsAsFactors=FALSE))
  }
  result<-tryCatch(withCallingHandlers(run(),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)list(error=conditionMessage(e),errorCall=paste(deparse(conditionCall(e)),collapse=' ')))
  result$warnings<-warnings
  result
}
.ctmaOptimizeCleanEngine <- function(engineFile,ownerFile) {
  if(!requireNamespace('ps',quietly=TRUE))return(invisible(NULL))
  if(file.exists(engineFile)) {
    x<-tryCatch(readRDS(engineFile),error=function(e)NULL)
    if(!is.null(x))try({h<-ps::ps_handle(x$pid)
      if(abs(ps::ps_create_time(h)-x$created)<.01&&grepl('julia',ps::ps_name(h),ignore.case=TRUE))ps::ps_terminate(h)},silent=TRUE)
  } else if(file.exists(ownerFile)) {
    x<-tryCatch(readRDS(ownerFile),error=function(e)NULL)
    if(!is.null(x))for(pid in ps::ps_pids())try({h<-ps::ps_handle(pid)
      # Strict ownership marker: this worker's unique temporary directory AND
      # its own Julia port-file prefix. Never match a generic Julia process.
      marker<-file.path(x$tempdir,paste0('juliaPort',x$pid))
      if(ps::ps_create_time(h)>=x$started-2&&grepl('julia',ps::ps_name(h),ignore.case=TRUE)&&
          any(grepl(marker,ps::ps_cmdline(h),fixed=TRUE)))ps::ps_terminate(h)
    },silent=TRUE)
  }
  invisible(NULL)
}
.ctmaOptimizeRun <- function(payload,timeout,logfile) {
  if(!requireNamespace('callr',quietly=TRUE)||!requireNamespace('ps',quietly=TRUE))
    .ctmaOptimizeStop('gambling=FALSE numerical checks need optional packages callr and ps. Install them explicitly first.')
  # Pin the namespaces actually loaded by the caller, not a different copy
  # that happens to appear first on .libPaths() after a library switch.
  payload$ctsemLibrary<-dirname(getNamespaceInfo(asNamespace('ctsem'),'path'))
  payload$cotimaLibrary<-dirname(getNamespaceInfo(asNamespace('CoTiMA2'),'path'))
  for(package in c('ctsem','CoTiMA2')) {
    library<-if(package=='ctsem')payload$ctsemLibrary else payload$cotimaLibrary
    if(!file.exists(file.path(library,package,'Meta','package.rds')))
      .ctmaOptimizeStop('Systematic numerical checks require an installed ',package,' package; load_all() source sessions are not supported.')
  }
  payload$ownerFile<-paste0(logfile,'.owner.rds');payload$engineFile<-paste0(logfile,'.engine.rds')
  unlink(c(payload$ownerFile,payload$engineFile))
  start<-proc.time()[['elapsed']]
  p<-callr::r_bg(function(payload,libs){
    .libPaths(unique(c(payload$cotimaLibrary,payload$ctsemLibrary,libs)))
    Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
    loadNamespace('ctsem',lib.loc=payload$ctsemLibrary)
    library('CoTiMA2',lib.loc=payload$cotimaLibrary,character.only=TRUE)
    worker<-get('.ctmaOptimizeCaseWorker',envir=asNamespace('CoTiMA2'),inherits=FALSE)
    worker(payload)
  },
    args=list(payload=payload,libs=.libPaths()),libpath=.libPaths(),stdout=logfile,
    stderr=paste0(logfile,'.err'),supervise=TRUE)
  on.exit({if(p$is_alive())try(p$kill_tree(),silent=TRUE);.ctmaOptimizeCleanEngine(payload$engineFile,payload$ownerFile)},add=TRUE)
  while(p$is_alive()) {
    if(proc.time()[['elapsed']]-start>=timeout) {
      try(p$kill_tree(),silent=TRUE)
      return(list(error='Pilot/full-check wall-time budget reached.',timedOut=TRUE,warnings=character()))
    }
    p$wait(timeout=250)
  }
  tryCatch(p$get_result(),error=function(e)list(error=conditionMessage(e),warnings=character()))
}
.ctmaOptimizeStatus <- function(result,pilot) {
  if(!is.null(result$error))return(if(isTRUE(result$timedOut))'timeout' else 'failed')
  a<-result$assessment
  if(!is.finite(a$loglik)||isTRUE(a$bad))return('unpromising')
  if(pilot)return('screened_inconclusive')
  if(isTRUE(a$uncertaintyValidated))'numerically_checked' else if(isTRUE(a$converged))'converged_uncertainty_unverified' else 'not_certified'
}
.ctmaOptimizeSystematic <- function(fit,init,backend,backendArgs,candidateStudies,
 exclusion,studyBlocks,referenceStudy,referenceCandidates,quickCheck,quickControl,
 timeBudget,maxCases,seed,checkpointDirectory,resume,coresToUse,reFits,iter,verbose,activeDirectory,finishsamples,planOnly=FALSE) {
  began<-proc.time()[['elapsed']]
  .ctmaOptimizeScalar(timeBudget,'timeBudget',FALSE,1)
  .ctmaOptimizeScalar(seed,'seed',TRUE,0);.ctmaOptimizeScalar(reFits,'reFits',TRUE,1)
  .ctmaOptimizeScalar(iter,'iter',TRUE,1)
  .ctmaOptimizeScalar(finishsamples,'finishsamples',TRUE,0)
  if(!is.logical(quickCheck)||length(quickCheck)!=1L||is.na(quickCheck)||
      !is.logical(resume)||length(resume)!=1L||is.na(resume)) .ctmaOptimizeStop('quickCheck and resume must be TRUE or FALSE.')
  if(!is.list(backendArgs)||(length(backendArgs)&&(is.null(names(backendArgs))||anyDuplicated(names(backendArgs)))))
    .ctmaOptimizeStop('backendArgs must be a named list without duplicate keys.')
  if(seed>2147483646) .ctmaOptimizeStop('seed must be <= 2147483646.')
  if(!is.logical(planOnly)||length(planOnly)!=1L||is.na(planOnly)) .ctmaOptimizeStop('planOnly must be TRUE or FALSE.')
  qc<-.ctmaOptimizeControls(quickControl)
  forbidden<-intersect(names(backendArgs),c('model','datalong','backend','optimize','priors','fit','inits','intoverpop','stationary','intoverstates','sampleControl'))
  if(length(forbidden)) .ctmaOptimizeStop('Systematic backendArgs cannot change model/target: ',paste(forbidden,collapse=', '))
  context<-.ctmaOptimizeContext(fit,init,backend)
  plan<-.ctmaOptimizePlan(context,candidateStudies,exclusion,studyBlocks,referenceStudy,referenceCandidates,maxCases)
  casePlan<-do.call(rbind,lapply(plan,function(x)data.frame(case=x$case,
    reference=x$reference,excluded=paste(x$excluded,collapse=', '),
    nstudies=length(x$included),nsubjects=x$nsubjects,stringsAsFactors=FALSE)))
  if(planOnly) {
    result<-list(bestFit=fit,bestConvergedFit=NULL,all_minus2ll=numeric(),
      summary=fit$summary,gambling=FALSE,backend=context$backend,planOnly=TRUE,
      casePlan=casePlan,automaticSelection=attr(plan,'automaticSelection'),
      sensitivity=list(plan=plan,cases=list()),caseDiagnostics=data.frame(),
      warnings=list(),errors=list(),checkpointDirectory=NULL,
      notes=c('Plan only: no numerical checks, exclusions or checkpoint writes were performed.',
        'Nominations are descriptive, reference-dependent sensitivity candidates, not grounds for deleting studies.'))
    class(result)<-c('ctmaOptimizeSystematic','CoTiMAFit')
    return(invisible(result))
  }
  cores<-.ctmaResolveCores(coresToUse,backend=context$backend,workflow='ctmaOptimizeFit')$resolved
  if(is.null(checkpointDirectory)) {
    if(resume) .ctmaOptimizeStop('resume=TRUE needs the original checkpointDirectory.')
    if(is.null(activeDirectory))activeDirectory<-if(!is.null(fit$activeDirectory))fit$activeDirectory else getwd()
    checkpointDirectory<-file.path(activeDirectory,paste0('ctmaOptimizeFit-systematic-',format(Sys.time(),'%Y%m%d-%H%M%S'),'-',Sys.getpid()))
  }
  if(!dir.exists(checkpointDirectory))dir.create(checkpointDirectory,recursive=TRUE,showWarnings=FALSE)
  if(!dir.exists(checkpointDirectory)||file.access(checkpointDirectory,2)!=0) .ctmaOptimizeStop('checkpointDirectory is not writable.')
  automatic<-attr(plan,'automaticSelection')
  manifest<-list(version=3L,context=context$fingerprint,plan=plan,backend=context$backend,
    cores=cores,ctsemRevision=utils::packageDescription('ctsem')$RemoteSha,
    implementation=.ctmaOptimizeHash(list(body(.ctmaOptimizeCase),body(.ctmaOptimizeCaseWorker),body(.ctmaOptimizeAutoCandidates),body(.ctmaOptimizeWorkerAssess))),
    controls=qc,quickCheck=quickCheck,backendArgs=backendArgs,reFits=reFits,iter=iter,seed=seed,finishsamples=finishsamples,
    CoTiMA2=as.character(utils::packageVersion('CoTiMA2')),ctsem=as.character(utils::packageVersion('ctsem')))
  path<-file.path(checkpointDirectory,'manifest.rds')
  if(resume) {
    if(!file.exists(path)||!identical(readRDS(path),manifest)) .ctmaOptimizeStop('Resume manifest differs: data, plan, versions or settings changed.')
  } else {
    if(file.exists(path)) .ctmaOptimizeStop('Existing checkpoint manifest: use resume=TRUE or a new directory.')
    saveRDS(manifest,path)
  }
  if(!is.null(automatic)) {
    utils::write.csv(automatic$studyDiagnostics,file.path(checkpointDirectory,'automaticStudyDiagnostics.csv'),row.names=FALSE)
    message('Automatic nominations: references [',paste(automatic$referenceCandidates,collapse=', '),
      ']; exclusion candidates [',paste(automatic$candidateStudies,collapse=', '),']. Reasons are in automaticSelection$studyDiagnostics.')
    for(note in automatic$notes)message('Note: ',note)
  }
  baseline<-.ctmaOptimizeWorkerAssess(context$native)
  expected<-baseline$loglik
  best<-fit;bestValue<-if(is.finite(expected))expected else -Inf
  bestConverged<-if(isTRUE(baseline$converged)&&!isFALSE(baseline$stationary))fit else NULL
  bestConvergedValue<-if(!is.null(bestConverged)&&is.finite(expected))expected else -Inf
  cases<-list();rows<-list();fullCriteria<-numeric();warnings<-errors<-list()
  deadline<-began+timeBudget
  execute<-function(prepared,caseNumber,attempt,pilot,point=NULL) {
    left<-deadline-proc.time()[['elapsed']]
    if(left<1)return(list(error='Global timeBudget exhausted.',timedOut=TRUE,warnings=character()))
    stem<-paste0(sprintf('%03d',caseNumber),'-',if(pilot)'pilot' else 'full','-',attempt)
    checkpoint<-file.path(checkpointDirectory,paste0(stem,'.rds'))
    if(resume&&file.exists(checkpoint))return(readRDS(checkpoint))
    # A case/start seed is stable under resume; original R RNG is untouched.
    key<-.ctmaOptimizeHash(list(prepared$case,attempt,pilot))
    caseSeed<-as.integer((as.double(strtoi(substr(key,1,7),16L))+seed)%%2147483646)+1L
    payload<-list(case=prepared,cores=cores,seed=caseSeed,attempt=attempt,
      maxiter=if(pilot)qc$maxiter else iter,pilot=pilot,jitter=qc$jitter,finishsamples=finishsamples,
      backendArgs=backendArgs,point=point,expectedLL=expected,verbose=verbose,
      argumentList=context$args)
    message('CoTiMA2 systematic check: ',prepared$case$case,'; ',if(pilot)'short pilot' else 'full numerical check',
      '; reference ',prepared$case$reference,'; ',length(prepared$case$included),' studies retained.')
    before<-proc.time()[['elapsed']]
    result<-.ctmaOptimizeRun(payload,min(left,if(pilot)qc$maxSeconds else qc$fullMaxSeconds),file.path(checkpointDirectory,paste0(stem,'.log')))
    result$elapsed<-proc.time()[['elapsed']]-before;result$seed<-caseSeed
    result$status<-.ctmaOptimizeStatus(result,pilot)
    saveRDS(result,checkpoint)
    result
  }
  record<-function(result,case,pilot,attempt) {
    if(is.null(result$status))result$status<-.ctmaOptimizeStatus(result,pilot)
    if(is.null(result$elapsed))result$elapsed<-0
    if(is.null(result$startLoglik))result$startLoglik<-NA_real_
    a<-result$assessment
    if(is.null(a))a<-list(loglik=NA_real_,maxGradient=NA_real_,converged=NA,weakDirections=NA,repaired=NA,uncertaintyValidated=FALSE)
    gain<-if(is.finite(result$startLoglik)&&is.finite(a$loglik))a$loglik-result$startLoglik else NA_real_
    rows[[length(rows)+1L]]<<-data.frame(case=case$case,stage=if(pilot)'pilot' else 'full',attempt=attempt,
      excluded=paste(case$excluded,collapse=', '),reference=case$reference,nstudies=length(case$included),nsubjects=case$nsubjects,
      loglik=a$loglik,within_case_gain=gain,max_gradient=a$maxGradient,
      converged=a$converged,curvature_assessed=if(is.null(a$curvatureAssessed))NA else a$curvatureAssessed,
      weak_directions=a$weakDirections,repaired=a$repaired,
      uncertainty_validated=isTRUE(a$uncertaintyValidated),status=result$status,
      elapsed=result$elapsed,stringsAsFactors=FALSE)
    if(length(result$warnings))warnings[[length(warnings)+1L]]<<-list(case=case$case,stage=if(pilot)'pilot' else 'full',messages=result$warnings)
    if(!is.null(result$error))errors[[length(errors)+1L]]<<-list(case=case$case,message=result$error)
    if(case$fullData&&!pilot&&is.finite(a$loglik)&&!isTRUE(a$bad)) {
      fullCriteria<<-c(fullCriteria,-2*a$loglik)
      if(a$loglik>bestValue){best<<-result$fit;bestValue<<-a$loglik}
      if(isTRUE(a$converged)&&!isFALSE(a$stationary)&&a$loglik>bestConvergedValue) {
        bestConverged<<-result$fit;bestConvergedValue<<-a$loglik
      }
    }
    utils::write.csv(do.call(rbind,rows),file.path(checkpointDirectory,'caseDiagnostics.csv'),row.names=FALSE)
  }
  promoted<-0L
  for(i in seq_along(plan)) {
    if(proc.time()[['elapsed']]>=deadline)break
    prepared<-tryCatch(.ctmaOptimizeCase(context,plan[[i]]),error=function(e)e)
    if(inherits(prepared,'condition')) {
      result<-list(error=conditionMessage(prepared),warnings=character(),status='failed',elapsed=0)
      record(result,plan[[i]],quickCheck,1L)
      cases[[plan[[i]]$case]]<-list(plan=plan[[i]],attempts=list(result));next
    }
    attempts<-list()
    for(attempt in seq_len(reFits)) {
      if(proc.time()[['elapsed']]>=deadline)break
      result<-execute(prepared,i,attempt,quickCheck)
      record(result,plan[[i]],quickCheck,attempt)
      attempts[[length(attempts)+1L]]<-result

    }
    cases[[plan[[i]]$case]]<-list(plan=plan[[i]],dataFingerprint=prepared$dataFingerprint,attempts=attempts)
  }
  # Finish the entire screening stage before any explicitly requested
  # promotion. Eligibility/order is numerical and prespecified, never based
  # on desired effect signs or across-dataset likelihood rankings.
  if(quickCheck&&qc$maxFullChecks>0L)for(i in seq_along(plan)) {
    if(promoted>=qc$maxFullChecks||proc.time()[['elapsed']]>=deadline)break
    saved<-cases[[plan[[i]]$case]]$attempts
    eligible<-which(vapply(saved,function(x)is.null(x$error)&&
      is.finite(x$assessment$loglik)&&!isTRUE(x$assessment$bad),logical(1)))
    if(!length(eligible))next
    chosen<-eligible[which.max(vapply(saved[eligible],function(x)x$assessment$loglik,numeric(1)))]
    prepared<-tryCatch(.ctmaOptimizeCase(context,plan[[i]]),error=function(e)NULL)
    if(is.null(prepared))next
    confirmed<-execute(prepared,i,chosen,FALSE,point=saved[[chosen]]$point)
    record(confirmed,plan[[i]],FALSE,chosen)
    cases[[plan[[i]]$case]]$attempts[[length(saved)+1L]]<-confirmed
    promoted<-promoted+1L
  }
  diagnostics<-if(length(rows))do.call(rbind,rows) else data.frame()
  result<-list(bestFit=best,bestConvergedFit=bestConverged,all_minus2ll=fullCriteria,
    summary=best$summary,usedStudyList=if(is.null(best$studyList))context$studies else best$studyList,usedTimeScale=context$scaleTime,
    usedScaleTI=context$args$scaleTI,gambling=FALSE,backend=context$backend,backendArgs=backendArgs,
    sensitivity=list(plan=plan,cases=cases),casePlan=casePlan,planOnly=FALSE,automaticSelection=automatic,caseDiagnostics=diagnostics,warnings=warnings,errors=errors,
    checkpointDirectory=normalizePath(checkpointDirectory),manifest=manifest,
    elapsed=proc.time()[['elapsed']]-began,
    notes=c('Full-data and excluded-study likelihoods are not ranked against one another.',
      'Pilots are screening calculations, not final estimates, confidence intervals, or calibrated convergence probabilities.',
      'Finite likelihood, small gradient or disappeared warning does not establish identification.',
      'bestFit remains the best full-data attempt; exclusions never silently replace the primary fit.'))
  class(result)<-c('ctmaOptimizeSystematic','CoTiMAFit')
  saveRDS(result,file.path(checkpointDirectory,'result.rds'))
  invisible(result)
}
.ctmaOptimizeDiagnosticTable <- function(object) {
  tab<-object$caseDiagnostics
  if(is.data.frame(tab)&&'weak_directions'%in%names(tab)) {
    if(!('curvature_assessed'%in%names(tab)))tab$curvature_assessed<-NA
    unknown<-is.na(tab$curvature_assessed)|!tab$curvature_assessed
    tab$weak_directions[unknown&!is.na(tab$weak_directions)&tab$weak_directions==0]<-NA_integer_
  }
  tab
}
#' @rdname ctmaOptimizeFit
#' @export
print.ctmaOptimizeSystematic <- function(x,...) {
  cat('CoTiMA2 systematic optimization/sensitivity checks\n')
  print(if(isTRUE(x$planOnly))x$casePlan else .ctmaOptimizeDiagnosticTable(x),row.names=FALSE,...)
  cat(paste(x$notes,collapse='\n'),'\n')
  invisible(x)
}
#' @rdname ctmaOptimizeFit
#' @export
summary.ctmaOptimizeSystematic <- function(object,...) {
  if(isTRUE(object$planOnly))object$casePlan else .ctmaOptimizeDiagnosticTable(object)
}
