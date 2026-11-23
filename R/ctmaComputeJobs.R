.ctmaComputeCoordinatorPlan <- function(compute, jobdir) {
  compute$execution <- "foreground"
  # Freeze inherited budgets before changing the coordinator's local budget.
  compute$remote <- lapply(compute$remote, function(h) {
    if (is.null(h$threads)) h$threads <- compute$threads
    h
  })
  co <- compute$coordinator
  if (!is.null(co)) {
    same <- vapply(compute$remote, function(h)
      identical(h$host, co$host) && identical(h$user, co$user) &&
      identical(h$port, co$port), logical(1))
    compute$local <- sum(vapply(compute$remote[same], function(h) h$workers, integer(1)))
    if (any(same)) compute$threads <- compute$remote[same][[1]]$threads
    compute$remote <- compute$remote[!same]
  }
  compute$coordinator <- NULL
  compute$directory <- file.path(jobdir, "checkpoints")
  compute
}
# Durable background jobs. The coordinator owns the R process and its workers.
# Collection reconnects to a job; it never assumes a dead Julia socket survived.
.ctmaComputeAtomic <- function(value,file) {
  tmp <- paste0(file,".tmp-",Sys.getpid());saveRDS(value,tmp)
  if(!file.rename(tmp,file)) stop("Cannot atomically save ",file)
}
.ctmaComputePackage <- function(compute) {
  source <- getOption("cotima.source")
  if(is.null(source)||!file.exists(file.path(source,"DESCRIPTION")))
    stop("Background jobs require options(cotima.source='/path/to/CoTiMA2') for an exact source snapshot.")
  root <- .ctmaComputeRoot(compute);stage <- tempfile("package-",tmpdir=root)
  dir.create(stage);dir.create(file.path(stage,"R"));dir.create(file.path(stage,"inst"));dir.create(file.path(stage,"inst","julia"))
  for(f in c("DESCRIPTION","NAMESPACE")) file.copy(file.path(source,f),file.path(stage,f))
  file.copy(list.files(file.path(source,"R"),full.names=TRUE,pattern="\\.[Rr]$"),file.path(stage,"R"))
  file.copy(list.files(file.path(source,"inst","julia"),full.names=TRUE),file.path(stage,"inst","julia"),recursive=TRUE)
  if (dir.exists(file.path(source,"data"))) {
    if (!file.copy(file.path(source,"data"),stage,recursive=TRUE))
      stop("Could not copy package data into the job snapshot.")
  }
  .ctmaComputeRun("git",c("-C",stage,"init","-q"));.ctmaComputeRun("git",c("-C",stage,"add","."))
  .ctmaComputeRun("git",c("-C",stage,"-c","user.name=CoTiMA","-c","user.email=compute@localhost","commit","-qm","Background package snapshot"))
  bundle <- paste0(stage,".bundle");.ctmaComputeRun("git",c("-C",stage,"bundle","create",bundle,"HEAD"))
  list(source=stage,bundle=bundle)
}
.ctmaComputeJobRunner <- function() c(
  'args <- commandArgs(TRUE); jobdir <- normalizePath(args[1],mustWork=TRUE)',
  'atomic <- function(x,n) { p<-file.path(jobdir,n); t<-paste0(p,".tmp-",Sys.getpid());saveRDS(x,t);if(!file.rename(t,p))stop("Atomic write failed") }',
  'lock <- file.path(jobdir,"lock"); if(!dir.create(lock,showWarnings=FALSE))stop("Job already locked; inspect prior process before resume")',
  'writeLines(as.character(Sys.getpid()),file.path(jobdir,"pid"))',
  'atomic(list(state="running",pid=Sys.getpid(),started=Sys.time()),"status.rds")',
  'tryCatch({',
  '  request <- readRDS(file.path(jobdir,"request.rds"))',
  '  .libPaths(c(file.path(jobdir,"library"),.libPaths()));library(CoTiMA2)',
  '  options(cotima.source=file.path(jobdir,"source"))',
  '  if (identical(request$args$backend, "ctsemjulia")) {',
  '    request$args$compute <- NULL',
  '    request$args$activeDirectory <- paste0(jobdir, "/")',
  '  }',
  '  if(!is.null(request$args$compute)) {',
  '    request$args$compute <- CoTiMA2:::.ctmaComputeCoordinatorPlan(request$args$compute,jobdir)',
  '  }',
  '  result <- do.call(getExportedValue("CoTiMA2",request$fun),request$args)',
  '  atomic(result,"result.rds");atomic(list(state="completed",finished=Sys.time()),"status.rds")',
  '},error=function(e) atomic(list(state="failed",error=conditionMessage(e),finished=Sys.time()),"status.rds"))',
  'unlink(lock,recursive=TRUE)')
.ctmaComputeSubmit <- function(fun,args,compute) {
  if(!is.null(compute$coordinator) && compute$local>0L)
    stop("Remote background coordinator requires local=0; list all compute hosts explicitly. The submitting laptop must not be an implicit dependency.")
  root <- .ctmaComputeRoot(compute);dir <- tempfile("job-",tmpdir=root);dir.create(dir,mode="0700")
  id <- basename(dir);snapshot <- .ctmaComputePackage(compute)
  saveRDS(list(fun=fun,args=args),file.path(dir,"request.rds"))
  writeLines(.ctmaComputeJobRunner(),file.path(dir,"run.R"))
  h <- compute$coordinator
  if(is.null(h)) {
    file.copy(snapshot$source,file.path(dir,"source"),recursive=TRUE)
    dir.create(file.path(dir,"library"))
    .ctmaComputeRun(file.path(R.home("bin"),"R"),c("CMD","INSTALL","--no-byte-compile","--no-help","--no-test-load","-l",file.path(dir,"library"),file.path(dir,"source")))
    if(!requireNamespace("callr",quietly=TRUE)) stop("Install callr for local background jobs")
    p <- callr::rscript_process$new(script=file.path(dir,"run.R"),cmdargs=dir,
      stdout=file.path(dir,"run.log"),stderr="2>&1",supervise=FALSE)
    job <- structure(list(id=id,local=dir,remote=NULL,coordinator=NULL,pid=p$get_pid()),class="ctmaJob")
    saveRDS(job,file.path(dir,"job.rds"));return(job)
  }
  h <- .ctmaComputeConnect(h);on.exit(.ctmaComputeClose(h),add=TRUE)
  target <- paste0(h$directory,"/",id)
  .ctmaComputeSsh(h,paste("umask 077 && mkdir -p",shQuote(target)))
  for(pair in list(c(snapshot$bundle,"source.bundle"),c(file.path(dir,"request.rds"),"request.rds"),c(file.path(dir,"run.R"),"run.R")))
    .ctmaComputeRun("scp",c("-o",paste0("ControlPath=",h$socket),"-o","BatchMode=yes","-P",as.character(h$port),pair[1],paste0(.ctmaComputeDestination(h),":",target,"/",pair[2])))
  .ctmaComputeSsh(h,paste("git -c advice.detachedHead=false clone -q",shQuote(paste0(target,"/source.bundle")),shQuote(paste0(target,"/source")),
    "&& mkdir -p",shQuote(paste0(target,"/library")),"&&", shQuote(h$rscript),
    "--vanilla -e",shQuote("quit(status = system2(file.path(R.home(\"bin\"), \"R\"), c(\"CMD\", \"INSTALL\", \"--no-byte-compile\", \"--no-help\", \"--no-test-load\", shQuote(commandArgs(TRUE))), stdout = \"\", stderr = \"\"))"),
    shQuote(paste0("--library=",target,"/library")),shQuote(paste0(target,"/source"))))
  command <- paste(shQuote(h$rscript),"--vanilla",shQuote(paste0(target,"/run.R")),shQuote(target),">",shQuote(paste0(target,"/run.log")),"2>&1")
  .ctmaComputeSsh(h,paste("tmux new-session -d -s",shQuote(id),shQuote(command)))
  raw <- compute$coordinator
  job <- structure(list(id=id,local=dir,remote=target,coordinator=raw),class="ctmaJob")
  saveRDS(job,file.path(dir,"job.rds"));job
}
.ctmaComputeJobFetch <- function(job,name) {
  stopifnot(inherits(job,"ctmaJob"),name %in% c("status.rds","result.rds","run.log","pid"))
  path <- file.path(job$local,name)
  if(is.null(job$coordinator)) return(if(file.exists(path)) path else NULL)
  h <- .ctmaComputeConnect(job$coordinator);on.exit(.ctmaComputeClose(h),add=TRUE)
  exists <- .ctmaComputeSsh(h,paste("if test -f",shQuote(paste0(job$remote,"/",name)),"; then echo yes; else echo no; fi"))
  if(!identical(tail(exists,1),"yes")) return(NULL)
  tmp <- tempfile(tmpdir=job$local)
  .ctmaComputeRun("scp",c("-o",paste0("ControlPath=",h$socket),"-o","BatchMode=yes","-P",as.character(h$port),
    paste0(.ctmaComputeDestination(h),":",job$remote,"/",name),tmp))
  file.rename(tmp,path);path
}
#' @rdname ctmaCompute
#' @export
ctmaJobStatus <- function(job) {
  path <- .ctmaComputeJobFetch(job,"status.rds")
  if(is.null(path)) return(list(state="starting",id=job$id))
  status <- readRDS(path);status$id <- job$id
  if(identical(status$state,"running") && !.ctmaComputeJobAlive(job,status$pid)) status$state <- "interrupted"
  status
}
#' @rdname ctmaCompute
#' @export
ctmaCollect <- function(job) {
  status <- ctmaJobStatus(job)
  if(status$state!="completed") stop("Job ",job$id," is ",status$state,
    if(!is.null(status$error)) paste0(": ",status$error) else "",call.=FALSE)
  path <- .ctmaComputeJobFetch(job,"result.rds")
  if(is.null(path)) stop("Completed job has no result; inspect its run.log")
  readRDS(path)
}
#' @rdname ctmaCompute
#' @export
ctmaResume <- function(job) {
  status <- ctmaJobStatus(job)
  if(status$state=="completed") return(job)
  if(status$state %in% c("running","starting")) stop("Job may still be running; resume will not duplicate it.")
  if(is.null(job$coordinator)) {
    if(.ctmaComputeJobAlive(job,status$pid)) stop("Job still running")
    unlink(file.path(job$local,"lock"),recursive=TRUE)
    p <- callr::rscript_process$new(script=file.path(job$local,"run.R"),cmdargs=job$local,
      stdout=file.path(job$local,"run.log"),stderr="2>&1",supervise=FALSE)
    job$pid <- p$get_pid()
  } else {
    h <- .ctmaComputeConnect(job$coordinator);on.exit(.ctmaComputeClose(h),add=TRUE)
    if(.ctmaComputeJobAlive(job,status$pid)) stop("Job still running")
    .ctmaComputeSsh(h,paste("rmdir",shQuote(paste0(job$remote,"/lock")),"2>/dev/null || true"))
    cmd <- paste(shQuote(h$rscript),"--vanilla",shQuote(paste0(job$remote,"/run.R")),shQuote(job$remote),">>",shQuote(paste0(job$remote,"/run.log")),"2>&1")
    .ctmaComputeSsh(h,paste("test ! -d",shQuote(paste0(job$remote,"/lock")),"&& tmux new-session -d -s",shQuote(job$id),shQuote(cmd)))
  }
  job
}
.ctmaComputeEntry <- function(compute,fun,call,environment,backend) {
  if(is.null(compute)) return(NULL)
  if(!inherits(compute,"ctmaCompute")) stop("Invalid compute plan")
  if (identical(compute$tasks, "studies") ||
      (identical(compute$tasks, "likelihood") && identical(tolower(backend[1L]), "ctsemjulia")))
    return(NULL)
  if (identical(tolower(backend[1L]), "ctsemjulia")) {
    if (!identical(compute$execution, "background") ||
        is.null(compute$coordinator) || compute$local != 0L ||
        length(compute$remote) != 1L || compute$remote[[1L]]$workers != 1L ||
        !identical(compute$tasks, "auto") ||
        !identical(.ctmaComputeDestination(compute$remote[[1L]]),
                   .ctmaComputeDestination(compute$coordinator)) ||
        !identical(compute$remote[[1L]]$port, compute$coordinator$port)) {
      stop("ctsemjulia compute requires one coordinator with one worker, local=0, tasks='auto' and execution='background'. Native fits do not split across hosts.",
           call. = FALSE)
    }
    args <- lapply(as.list(call)[-1L], eval, envir = environment)
    args$backend <- "ctsemjulia"
    host_threads <- compute$remote[[1L]]$threads
    args$coresToUse <- if (is.null(host_threads)) compute$threads else host_threads
    return(.ctmaComputeSubmit(fun, args, compute))
  }
  if(compute$execution!="background") return(NULL)
  args <- lapply(as.list(call)[-1L],eval,envir=environment)
  .ctmaComputeSubmit(fun,args,compute)
}

.ctmaComputeJobAlive <- function(job,pid=NULL) {
  if(is.null(pid)||!is.numeric(pid)||length(pid)!=1L||is.na(pid)||pid<1) return(TRUE)
  if(is.null(job$coordinator)) {
    if(!requireNamespace("ps",quietly=TRUE)) return(TRUE)
    return(tryCatch(ps::ps_is_running(ps::ps_handle(pid)),error=function(e) FALSE))
  }
  h <- .ctmaComputeConnect(job$coordinator);on.exit(.ctmaComputeClose(h),add=TRUE)
  result <- .ctmaComputeSsh(h,paste("if kill -0",as.integer(pid),"2>/dev/null; then echo yes; else echo no; fi"))
  identical(tail(result,1),"yes")
}
#' @export
print.ctmaJob <- function(x,...) {
  cat("CoTiMA background job:",x$id,"\nCoordinator:",if(is.null(x$coordinator)) "local" else x$coordinator$host,
      "\nSaved handle:",file.path(x$local,"job.rds"),"\nUse ctmaJobStatus(), ctmaCollect(), or ctmaResume().\n")
  invisible(x)
}
