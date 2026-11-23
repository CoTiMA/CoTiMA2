# Connection and execution plans contain configuration only, never passwords.
.ctmaComputeInt <- function(x, name, zero=FALSE) {
  if (!is.numeric(x) || length(x)!=1L || is.na(x) || !is.finite(x) ||
      x!=floor(x) || x < if(zero) 0 else 1 || x>.Machine$integer.max)
    stop(name," must be a whole number >= ",if(zero) 0 else 1,call.=FALSE)
  as.integer(x)
}
.ctmaComputeString <- function(x, name, pattern=NULL) {
  if (!is.character(x)||length(x)!=1L||is.na(x)||!nzchar(x)||grepl("[\r\n]",x) ||
      (!is.null(pattern)&&!grepl(pattern,x))) stop("Invalid ",name,call.=FALSE)
  x
}
#' @rdname ctmaCompute
#' @export
ctmaHost <- function(host, user=NULL, workers=1L, port=22L,
    authentication=c("auto","ask"), julia="julia", project=NULL,
    directory="~/CoTiMA2-Benchmark", rscript="Rscript", depot=NULL, threads=NULL) {
  host <- .ctmaComputeString(host,"host","^[A-Za-z0-9][A-Za-z0-9._:-]*$")
  if(!is.null(user)) user <- .ctmaComputeString(user,"user","^[A-Za-z0-9_][A-Za-z0-9_.-]*$")
  port <- .ctmaComputeInt(port,"port"); if(port>65535) stop("Invalid port")
  for(n in c("julia","directory","rscript")) .ctmaComputeString(get(n),n)
  .ctmaComputeString(directory,"directory","^[A-Za-z0-9_~./-]+$")
  for(n in c("project","depot")) if(!is.null(get(n))) .ctmaComputeString(get(n),n)
  if(!is.null(threads)) threads <- .ctmaComputeInt(threads,"threads")
  structure(list(threads=threads,host=host,user=user,workers=.ctmaComputeInt(workers,"workers"),
    port=port,authentication=match.arg(authentication),julia=julia,project=project,
    directory=directory,rscript=rscript,depot=depot),class="ctmaHost")
}
#' @rdname ctmaCompute
#' @export
ctmaCompute <- function(workers=c(local=1L), threads=1L, local=NULL, remote=NULL,
    tasks=c("auto","likelihood","starts","studies"), tunnel="ssh",
    reconnect=TRUE, checkpoint=TRUE, execution=c("foreground","background"),
    coordinator=NULL, retries=1L, directory=NULL) {
  if(!identical(tunnel,"ssh")) stop("Only tunnel='ssh' is supported.",call.=FALSE)
  if(is.null(names(workers))||any(!nzchar(names(workers)))||anyDuplicated(names(workers)))
    stop("workers must have unique host names, including 'local' if used.")
  hosts <- list(); local_count <- 0L
  for(n in names(workers)) {
    k <- .ctmaComputeInt(workers[[n]],paste("workers",n),zero=TRUE)
    if(n=="local") local_count <- k else if(k>0L) {
      profile <- getOption("cotima.hosts",list())[[n]]
      if(is.null(profile)) profile <- ctmaHost(n,workers=k)
      if(!inherits(profile,"ctmaHost")) stop("Invalid cotima.hosts profile: ",n)
      profile$workers <- k; hosts[[length(hosts)+1L]] <- profile
    }
  }
  if(!is.null(local)) local_count <- .ctmaComputeInt(local,"local",zero=TRUE)
  if(!is.null(remote)) {
    if(inherits(remote,"ctmaHost")) remote <- list(remote)
    if(!is.list(remote)||!all(vapply(remote,inherits,logical(1),"ctmaHost"))) stop("remote must contain ctmaHost objects")
    hosts <- c(hosts,remote)
  }
  keys <- vapply(hosts,function(h) paste(h$host,h$user,h$port),character(1))
  if(anyDuplicated(keys)) stop("Duplicate remote host; combine worker counts in one profile")
  if(local_count+sum(vapply(hosts,`[[`,integer(1),"workers"))<1L) stop("At least one worker is required")
  for(n in c("reconnect","checkpoint")) if(!is.logical(get(n))||length(get(n))!=1L||is.na(get(n))) stop(n," must be TRUE or FALSE")
  if(!is.null(coordinator)&&!inherits(coordinator,"ctmaHost")) stop("coordinator must be a ctmaHost")
  structure(list(local=local_count,remote=hosts,threads=.ctmaComputeInt(threads,"threads"),
    tasks=match.arg(tasks),tunnel=tunnel,reconnect=reconnect,checkpoint=checkpoint,
    execution=match.arg(execution),coordinator=coordinator,retries=.ctmaComputeInt(retries,"retries",TRUE),
    directory=directory),class="ctmaCompute")
}
.ctmaComputeRun <- function(exe,args, output=TRUE) {
  ans <- suppressWarnings(system2(exe, args=vapply(args,shQuote,character(1)),
    stdout=if(output) TRUE else "",stderr=if(output) TRUE else ""))
  code <- if(output) attr(ans,"status") else ans
  if(!is.null(code)&&code!=0L) stop("Command failed (",basename(exe),"): ",
    if(output) paste(ans,collapse="\n") else code,call.=FALSE)
  ans
}
.ctmaComputeDestination <- function(h) paste0(if(is.null(h$user)) "" else paste0(h$user,"@"),h$host)
.ctmaComputeSsh <- function(h,command) {
  args <- c("-S",h$socket,"-p",as.character(h$port),"-o","BatchMode=yes",
    "-o","ConnectTimeout=10",.ctmaComputeDestination(h),command)
  .ctmaComputeRun("ssh",args)
}
.ctmaComputeConnect <- function(h) {
  if(.Platform$OS.type=="windows") stop("SSH compute currently requires macOS or Linux.")
  h$socket <- tempfile("ctma-ssh-",tmpdir="/tmp")
  args <- c("-M","-S",h$socket,"-f","-N","-T","-p",as.character(h$port),
    "-o","ControlPersist=600","-o","ServerAliveInterval=15","-o","ServerAliveCountMax=3",
    "-o","ConnectTimeout=10","-o","StrictHostKeyChecking=ask")
  if(h$authentication=="auto") args <- c(args,"-o","BatchMode=yes") else {
    if(!interactive()) stop("authentication='ask' requires an interactive R session; use an unlocked SSH key for unattended jobs.")
    args <- c(args,"-o","BatchMode=no","-o","PreferredAuthentications=keyboard-interactive,password")
  }
  # OpenSSH owns the password dialogue; R never receives or serializes the secret.
  ask <- NULL
  old <- Sys.getenv(c("SSH_ASKPASS","SSH_ASKPASS_REQUIRE","DISPLAY"),unset=NA)
  on.exit({for(n in names(old)) if(is.na(old[[n]])) Sys.unsetenv(n) else do.call(Sys.setenv,setNames(list(old[[n]]),n));if(!is.null(ask)) unlink(ask)},add=TRUE)
  if(h$authentication=="ask" && Sys.info()[["sysname"]]=="Darwin") {
    ask <- tempfile("ctma-askpass-");writeLines(c("#!/bin/sh",
      "exec /usr/bin/osascript - \"$1\" <<'APPLESCRIPT'",
      "on run argv", "set promptText to item 1 of argv",
      "set reply to display dialog promptText default answer \"\" with hidden answer buttons {\"Cancel\", \"OK\"} default button \"OK\"",
      "return text returned of reply", "end run", "APPLESCRIPT"),ask)
    Sys.chmod(ask,"0700");Sys.setenv(SSH_ASKPASS=ask,SSH_ASKPASS_REQUIRE="force",DISPLAY=":0")
  }
  .ctmaComputeRun("ssh",c(args,.ctmaComputeDestination(h)),output=FALSE)
  connected <- FALSE
  on.exit(if(!connected) .ctmaComputeClose(h),add=TRUE)
  home <- tail(.ctmaComputeSsh(h,"printf '%s\\n' \"$HOME\""),1)
  for(n in c("directory","project","julia","depot")) if(!is.null(h[[n]]))
    h[[n]] <- gsub("~/",paste0(home,"/"),h[[n]],fixed=TRUE)
  h$home <- home
  connected <- TRUE
  h
}
.ctmaComputeClose <- function(h) {
  if(!is.null(h$socket)) try(suppressWarnings(system2("ssh",vapply(c("-S",h$socket,"-O","exit",
    "-p",as.character(h$port),.ctmaComputeDestination(h)),shQuote,character(1)),stdout=FALSE,stderr=FALSE)),silent=TRUE)
}
#' @rdname ctmaCompute
#' @export
ctmaCheckCompute <- function(compute) {
  stopifnot(inherits(compute,"ctmaCompute"))
  checks <- list()
  for(h in compute$remote) {
    h <- .ctmaComputeConnect(h)
    info <- tryCatch(.ctmaComputeSsh(h,paste(shQuote(h$julia),"--startup-file=no --version")),finally=.ctmaComputeClose(h))
    checks[[length(checks)+1L]] <- data.frame(host=h$host,workers=h$workers,threads=if(is.null(h$threads)) compute$threads else h$threads,julia=paste(info,collapse=" "))
  }
  if(compute$local>0L) checks <- c(list(data.frame(host="local",workers=compute$local,threads=compute$threads,
    julia=paste(.ctmaComputeRun("julia",c("--startup-file=no","--version")),collapse=" "))),checks)
  do.call(rbind,checks)
}
.ctmaComputeBind <- function(compute, backend, backendArgs, workflow) {
  if(is.null(compute)) return(backendArgs)
  if(!inherits(compute,"ctmaCompute")) stop("compute must be created by ctmaCompute()")
  if (identical(compute$tasks, "likelihood") && tolower(backend[1L]) == "ctsemjulia") {
    if (workflow != "ctmaFit" || compute$execution != "foreground" || !is.null(compute$coordinator))
      stop("Native distributed likelihood requires foreground ctmaFit().")
    return(backendArgs)
  }
  if (identical(compute$tasks, "studies")) {
    if (workflow != "ctmaInit" || tolower(backend[1L]) != "ctsemjulia" ||
        compute$execution != "foreground" || !is.null(compute$coordinator))
      stop("tasks='studies' requires foreground ctmaInit(backend='ctsemjulia').")
    return(backendArgs)
  }
  if(tolower(backend[1L])!="julia") stop("compute currently supports backend='julia' only")
  if(!is.null(backendArgs$compute)) stop("Supply compute only once")
  if(compute$tasks=="auto") compute$tasks <- "starts"
  if (compute$tasks == "chains") stop("Distributed chain scheduling is not supported.")
  for(n in c("parallel_starts","parallel_chains","julia_threads","subject_threads","blas_threads","threads_per_start"))
    if(!is.null(backendArgs[[n]]) && !identical(as.numeric(backendArgs[[n]]),1)) stop("compute owns resource limits; remove backendArgs$",n)
  backendArgs$compute <- compute
  backendArgs$julia_threads <- 1L;backendArgs$subject_threads <- 1L;backendArgs$blas_threads <- 1L
  backendArgs$parallel_starts <- 1L;backendArgs$threads_per_start <- 1L
  backendArgs
}
.ctmaComputeRoot <- function(compute) {
  root <- compute$directory
  if(is.null(root)) root <- file.path(tools::R_user_dir("CoTiMA2","cache"),"compute")
  dir.create(root,recursive=TRUE,showWarnings=FALSE,mode="0700")
  normalizePath(root,mustWork=TRUE)
}
.ctmaComputeStage <- function(compute,julia_file) {
  root <- .ctmaComputeRoot(compute)
  helper <- file.path(dirname(julia_file),"CtsemCompute.jl")
  if(!file.exists(helper)) helper <- system.file("julia","CtsemCompute.jl",package="CoTiMA2")
  if(!nzchar(helper)||!file.exists(helper)) stop("CtsemCompute.jl is missing")
  identity <- tempfile();writeLines(unname(tools::md5sum(c(julia_file,helper))),identity)
  id <- unname(tools::md5sum(identity));unlink(identity)
  stage <- file.path(root,paste0("source-",id))
  if(!dir.exists(file.path(stage,".git"))) {
    dir.create(stage,recursive=TRUE,showWarnings=FALSE)
    file.copy(c(julia_file,helper),file.path(stage,c("CtsemJulia.jl","CtsemCompute.jl")),overwrite=TRUE)
    .ctmaComputeRun("git",c("-C",stage,"init","-q"))
    .ctmaComputeRun("git",c("-C",stage,"add","CtsemJulia.jl","CtsemCompute.jl"))
    .ctmaComputeRun("git",c("-C",stage,"-c","user.name=CoTiMA","-c","user.email=compute@localhost","commit","-qm","Compute snapshot"))
  }
  bundle <- file.path(root,paste0(id,".bundle"))
  .ctmaComputeRun("git",c("-C",stage,"bundle","create",bundle,"HEAD"))
  list(id=id,source=stage,bundle=bundle)
}
.ctmaComputePrepare <- function(compute,julia_file) {
  if(compute$local*compute$threads > .ctmaAvailableCores()) stop("Local compute budget exceeds available cores")
  stage <- .ctmaComputeStage(compute,julia_file);connections <- list();config <- list()
  success <- FALSE
  on.exit(if(!success) lapply(connections,.ctmaComputeClose),add=TRUE)
  if(compute$local>0L) config[[1L]] <- list(host="local",workers=compute$local,
    threads=compute$threads,project=normalizePath(dirname(julia_file)),source=stage$source,
    julia="",port=22L,socket="",depot="")
  for(raw in compute$remote) {
    h <- .ctmaComputeConnect(raw);connections[[length(connections)+1L]] <- h
    if(is.null(h$project)) {
      code <- 'cat(system.file("julia",package="CoTiMA2"))'
      h$project <- tail(.ctmaComputeSsh(h,paste(shQuote(h$rscript),"--vanilla -e",shQuote(code))),1)
      if(!length(h$project)||!nzchar(h$project)) stop("Specify a Julia project for host ",h$host)
    }
    target <- paste0(h$directory,"/compute-source-",stage$id)
    .ctmaComputeSsh(h,paste("mkdir -p",shQuote(h$directory)))
    remote_bundle <- paste0(target,".bundle")
    .ctmaComputeRun("scp",c("-o",paste0("ControlPath=",h$socket),"-o","BatchMode=yes","-P",as.character(h$port),
      stage$bundle,paste0(.ctmaComputeDestination(h),":",remote_bundle)))
    .ctmaComputeSsh(h,paste("test -d",shQuote(paste0(target,"/.git")),"|| git -c advice.detachedHead=false clone -q",shQuote(remote_bundle),shQuote(target)))
    config[[length(config)+1L]] <- list(host=.ctmaComputeDestination(h),workers=h$workers,
      threads=if(is.null(h$threads)) compute$threads else h$threads,project=h$project,
      source=target,julia=h$julia,port=h$port,socket=h$socket,depot=if(is.null(h$depot)) "" else h$depot)
  }
  success <- TRUE
  list(config=config,connections=connections,source=stage$source)
}
.ctmaComputeOptimize <- function(compute,julia_file,inits,maxiters,chunksize,
    gradient_tolerance,optimizer_initial_step,multistart,multistart_jitter,seed,
    function_relative_tolerance,parameter_relative_tolerance,show_trace) {
  if(!inherits(compute,"ctmaCompute")) stop("Invalid compute plan")
  if(compute$execution!="foreground") stop("Use a background job through ctmaFit()/ctmaInit(), not ctmaJuliaFit() directly")
  if(compute$tasks=="auto") compute$tasks <- "starts"
  p <- .ctmaComputePrepare(compute,julia_file)
  on.exit(lapply(p$connections,.ctmaComputeClose),add=TRUE)
  JuliaCall::julia_assign("ctma_compute_source",p$source)
  JuliaCall::julia_eval('using SHA; path=joinpath(ctma_compute_source,"CtsemCompute.jl"); if !isdefined(Main,:CtsemCompute) || CtsemCompute.CtsemJulia !== CtsemJulia || CtsemCompute.SOURCE_HASH != bytes2hex(SHA.sha256(read(path))); include(path); end')
  hashes <- JuliaCall::julia_eval('[bytes2hex(SHA.sha256(read(joinpath(ctma_compute_source,f)))) for f in ["CtsemJulia.jl","CtsemCompute.jl"]]')
  for(i in seq_along(p$config)) p$config[[i]]$hashes <- as.character(hashes)
  JuliaCall::julia_assign("ctma_compute_config",p$config)
  JuliaCall::julia_assign("ctma_compute_start",as.numeric(inits))
  # Checkpoint names identify data/start/settings; Julia additionally checks source hashes.
  keyfile <- tempfile();saveRDS(list(JuliaCall::julia_eval("ctma_sd"),inits,maxiters,chunksize,
    gradient_tolerance,optimizer_initial_step,multistart,multistart_jitter,seed,compute$tasks,
    function_relative_tolerance,parameter_relative_tolerance,basename(p$source)),keyfile)
  key <- unname(tools::md5sum(keyfile));unlink(keyfile)
  checkpoint <- if(compute$checkpoint) file.path(.ctmaComputeRoot(compute),paste0("fit-",key,".jls")) else ""
  JuliaCall::julia_assign("ctma_compute_checkpoint",checkpoint)
  text <- sprintf(paste0('CtsemCompute.optimize(ctma_compute_start,ctma_sd,ctma_compute_config;',
    'tasks="%s",starts=%d,jitter=%.17g,seed=%d,maxiters=%d,chunksize=%d,',
    'gradient_tolerance=%.17g,initial_raw_step=%.17g,show_trace=%s,checkpoint=ctma_compute_checkpoint,',
    'retries=%d,function_relative_tolerance=%.17g,parameter_relative_tolerance=%.17g)'),
    compute$tasks,as.integer(multistart),multistart_jitter,if(is.null(seed)) 1L else as.integer(seed),
    as.integer(maxiters),as.integer(chunksize),gradient_tolerance,optimizer_initial_step,
    if(show_trace) "true" else "false",if(compute$reconnect) compute$retries else 0L,
    function_relative_tolerance,parameter_relative_tolerance)
  JuliaCall::julia_eval(text)
}

#' @export
print.ctmaCompute <- function(x,...) {
  cat("CoTiMA compute plan:",x$tasks,"/",x$execution,"via",x$tunnel,"\n")
  cat("Local workers:",x$local,"with",x$threads,"thread(s) each\n")
  for(h in x$remote) cat(h$host,":",h$workers,"workers x",if(is.null(h$threads)) x$threads else h$threads,"thread(s)\n")
  if(!is.null(x$coordinator)) cat("Coordinator:",x$coordinator$host,"\n")
  invisible(x)
}
