# Reporting policy: preserve draw summaries and explicitly select point estimates.
.ctmaPointChoice <- function(pointEstimate) {
  match.arg(tolower(pointEstimate[1L]),c("auto","mean","mode"))
}

.ctmaPointMethod <- function(fit, optimize = NULL, priors = NULL) {
  native <- fit
  if(inherits(fit,"CoTiMAFit")) {
    candidate <- .ctmaTablesNative(fit)
    if(!is.null(candidate)) native <- candidate
  }
  if (!is.null(native$sample) || identical(native$estimation,"Bayes")) return("Bayesian")
  if (identical(native$estimation,"ML")) return("ML")
  candidates <- list(native$args$resolved,native$args$input,native$argumentList,
    native$argumentListOfCtmaInit,fit$argumentList,native$args,
    list(priors=native$standata$priors))
  pick <- function(key, fallback) {
    for(a in candidates) {
      value <- a[[key]]
      if(length(value)==1L && (is.logical(value)||is.numeric(value)) &&
          !is.na(value) && value%in%c(0,1)) return(as.logical(value))
    }
    fallback
  }
  optimize <- pick("optimize",optimize)
  priors <- pick("priors",priors)
  if (isFALSE(optimize)) return("Bayesian")
  if (isTRUE(optimize)) return(if(isTRUE(priors)) "MAP" else if(isFALSE(priors)) "ML" else "optimized")
  "unknown"
}

.ctmaStatColumn <- function(table, statistic) {
  index <- match(tolower(statistic),tolower(colnames(table)))
  if (is.na(index)) stop("Missing statistic column: ",statistic,call.=FALSE)
  colnames(table)[index]
}

.ctmaPointColumn <- function(table, pointEstimate="auto", optimize=TRUE) {
  pointEstimate <- .ctmaPointChoice(pointEstimate)
  mode <- which(tolower(colnames(table))=="mode")
  if (pointEstimate=="mode" && !length(mode)) stop(
    "No optimized point estimates are available for this table. Use pointEstimate='mean' or a fit retaining its optimized values.",call.=FALSE)
  if (length(mode) && (pointEstimate=="mode" || (pointEstimate=="auto" && isTRUE(optimize))))
    return(colnames(table)[mode])
  .ctmaStatColumn(table,"mean")
}

.ctmaPointQuotient <- function(table, pointEstimate="auto") {
  # Explicit mean selects the historical mean/SD ratio, never mode/SD.
  if (.ctmaPointChoice(pointEstimate)!="mean") return(rep(NA_real_,nrow(table)))
  if(is.matrix(table)) table<-as.data.frame(table,check.names=FALSE)
  mu<-table[[.ctmaStatColumn(table,"mean")]]
  sd<-table[[.ctmaStatColumn(table,"sd")]]
  ifelse(is.finite(sd)&sd>0,mu/sd,NA_real_)
}

.ctmaPointNote <- function(method, available=TRUE, legacy=FALSE) {
  if (method=="Bayesian") return(paste0(
    "Mean, SD and quantiles summarize posterior draws. SD measures posterior uncertainty. ",
    "Mean/SD, where displayed, is a descriptive ratio and is not a general Bayesian significance test."))
  if (!method %in% c("ML","MAP","optimized")) return(NULL)
  point <- if(available) paste0("The mode column reports ",method,
    " point estimates. Derived quantities are evaluated at the optimized parameter vector. ") else
    "Optimized values are unavailable in this stored summary; point estimates fall back to draw means. "
  paste0(point,"Mean, SD and quantiles summarize stored uncertainty draws. ",
    "SD measures draw spread. A mode/SD Wald test requires an appropriate standard error and regularity conditions. ",
    if(legacy) "Historical Mean/SD ratios were explicitly requested." else "")
}

.ctmaOrderMode <- function(table, legacy=FALSE) {
  if (!is.data.frame(table) && !is.matrix(table)) return(table)
  columns<-colnames(table)
  mode<-which(tolower(columns)=="mode")
  if(length(mode)) {
    stat<-which(tolower(columns)%in%c("mode","mean","sd","se","2.5%","50%","97.5%","z","n_eff","rhat"))
    first<-min(stat)
    columns<-append(columns[-mode],columns[mode],after=first-1L)
    table<-table[,columns,drop=FALSE]
  }
  if(!legacy) colnames(table)[tolower(colnames(table))=="z"]<-"mean_over_sd"
  if(legacy) colnames(table)[tolower(colnames(table))=="mean_over_sd"]<-"z"
  if("Tvalues"%in%colnames(table)) table[,"Tvalues"] <-
    .ctmaPointQuotient(table,if(legacy) "mean" else "auto")
  table
}

.ctmaPointSummary <- function(answer, fit, pointEstimate="auto", optimize=NULL, priors=NULL) {
  if(!is.list(answer) || is.data.frame(answer)) return(answer)
  choice<-.ctmaPointChoice(pointEstimate)
  method<-.ctmaPointMethod(fit,optimize,priors)
  if(choice=="mode" && method=="Bayesian")
    stop("pointEstimate='mode' requires an optimized ML/MAP fit; Bayesian sampling summaries use posterior means or medians.",call.=FALSE)
  table_name <- if(is.data.frame(answer$parmatrices)) "parmatrices" else
    if(is.data.frame(answer$estimates)) "estimates" else NULL
  native <- fit
  if(inherits(fit,"CoTiMAFit")) {
    candidate <- .ctmaTablesNative(fit)
    if(!is.null(candidate)) native <- candidate
  }
  if(!is.null(table_name)) {
    result <- .ctmaEnsurePointMode(list(parmatrices=answer[[table_name]]),native,
      method %in% c("ML","MAP","optimized") && choice!="mean")
    answer[[table_name]] <- result$parmatrices
  }
  tables<-which(vapply(answer,function(x) (is.data.frame(x)||is.matrix(x)) &&
    any(tolower(colnames(x))=="mean"),logical(1)))
  available<-any(vapply(answer[tables],function(x) any(tolower(colnames(x))=="mode"),logical(1)))
  if(choice=="mode" && length(tables) && !available)
    stop("No optimized mode values are available in this summary.",call.=FALSE)
  if(length(answer$parameterTables)) available <- available ||
    any(vapply(answer$parameterTables,function(t) "mode"%in%colnames(t),logical(1)))
  for(i in tables) answer[[i]]<-.ctmaOrderMode(answer[[i]],legacy=choice=="mean")
  note<-.ctmaPointNote(method,available,choice=="mean")
  answer$legacyPointEstimateNote<-.ctmaLegacyPointNotice(fit,choice)
  answer$pointEstimateNote<-note
  answer$pointEstimate<-list(requested=choice,method=method,
    selected=if(choice=="mode" || (choice=="auto"&&available&&method%in%c("ML","MAP","optimized"))) "mode" else "mean")
  answer
}

.ctmaPrintPointNote <- function(answer) {
  if(length(answer$legacyPointEstimateNote)) .ctma_print_message(answer$legacyPointEstimateNote)
  if(length(answer$pointEstimateNote)) .ctma_print_message(answer$pointEstimateNote)
  invisible(NULL)
}

.ctmaLegacyPointNotice <- function(fit, choice="auto") {
  if (!is.null(fit$argumentList$pointEstimate)) return(NULL)
  cotima <- inherits(fit,c("CoTiMAFit","ctmaCtsemJuliaFit","ctmaJuliaFit")) ||
    !is.null(fit$argumentListOfCtmaInit) || !is.null(fit$cotima_call) || !is.null(fit$cotima_study)
  if (!cotima) return(NULL)
  paste0("This saved object predates the pointEstimate reporting policy and used draw means in the earlier CoTiMA summaries. ",
    "The current package selects optimized values ('mode') for ML/MAP when available. ",
    "Use pointEstimate = 'mean' to select the earlier point estimates. ",
    "Bayesian point estimates continue to use posterior means by default.")
}

.ctmaEnsurePointMode <- function(answer, fit, optimize=TRUE) {
  if(inherits(fit,"ctmaJuliaFit") && !is.null(fit$ctsem_fit)) {
    bridge<-fit$ctsem_fit
    bridge$stanfit$rawest<-fit$coefficients$raw
    fit<-bridge
  }
  table<-answer$parmatrices
  if (!isTRUE(optimize)||!is.data.frame(table)||"mode"%in%names(table)||!is.null(fit$sample)) return(answer)
  point<-tryCatch({
    if(inherits(fit,"ctJuliaFit")) {
      f<-getFromNamespace(".ctBackendSummary","ctsem")
      result<-f(fit,digits=12,residualcov=FALSE)
      d<-result$parmatrices
      if(!"mode"%in%names(d) && is.null(fit$estimate$rawposterior))
        d$mode<-d[[.ctmaStatColumn(d,"mean")]]
      d
    } else if(inherits(fit,"ctStanFit") && !is.null(fit$stanfit$rawest) &&
        !is.null(fit$stanmodel) && !is.null(fit$standata)) {
      raw<-matrix(as.numeric(fit$stanfit$rawest),nrow=1L)
      transform<-getFromNamespace("stan_constrainsamples","ctsem")
      draws<-transform(sm=fit$stanmodel,standata=fit$standata,samples=raw,
        cores=1L,savescores=FALSE,savesubjectmatrices=FALSE,dokalman=FALSE)
      clone<-fit
      clone$stanfit$rawposterior<-raw
      clone$stanfit$transformedpars<-draws
      matrices<-ctsem::ctSummaryMatrices(clone,calcfunc=mean,calcfuncargs=list())
      d<-ctsem::ctModelUnlist(matrices,matnames=names(matrices))
      d$mode<-d$value
      d
    } else NULL
  },error=function(e) NULL)
  if(is.data.frame(point)&&"mode"%in%names(point)&&
      all(c("matrix","row","col")%in%names(table))&&
      all(c("matrix","row","col")%in%names(point))) {
    key<-function(x)paste(x$matrix,x$row,x$col,sep="/")
    index<-match(key(table),key(point))
    if(!anyNA(index)) answer$parmatrices$mode<-point$mode[index]
  }
  answer
}

.ctmaPointInitSummary <- function(object, choice) {
  fits<-object$studyFitList
  if (length(object$n.studies)!=1L) return(object)
  if (!is.list(fits)||inherits(fits,c("ctJuliaFit","ctStanFit","ctmaJuliaFit")) ||
      length(fits)!=object$n.studies) return(object)
  original<-object$summary
  if(!is.matrix(original$estimates) && !is.data.frame(original$estimates)) return(object)
  tables<-vector("list",length(fits))
  selected<-character(length(fits))
  scale<-.ctmaSummaryTimeFactor(object)
  for(i in seq_along(fits)) {
    s<-fits[[i]]$resultsSummary
    if(!is.list(s)||!is.data.frame(s$parmatrices)) next
    s<-.ctmaPointSummary(s,fits[[i]],choice,
      object$argumentList$optimize,object$argumentList$priors)
    t<-s$parmatrices
    tables[[i]]<-t
    selected[i]<-.ctmaPointColumn(t,choice,s$pointEstimate$method%in%c("ML","MAP","optimized"))
    rows<-if("matrix"%in%names(t)) which(t$matrix=="DRIFT") else which(rownames(t)=="DRIFT")
    if(length(rows)!=object$n.latent^2) next
    digits<-object$argumentList$digits
    if(is.null(digits)) digits<-4L
    values<-round(t[rows,selected[i]],digits)
    columns<-seq.int(2L,2L*length(values),by=2L)
    original$estimates[i,columns]<-values
    if(!is.null(original$drift_estimates_original_time_scale))
      original$drift_estimates_original_time_scale[i,columns]<-values*scale
  }
  original$parameterTables<-tables
  original$pointEstimateNote<-.ctmaPointNote(.ctmaPointMethod(object),
    all(selected=="mode"),choice=="mean")
  original$legacyPointEstimateNote<-.ctmaLegacyPointNotice(object,choice)
  original$pointEstimate<-list(requested=choice,selected=selected)
  attr(original,"ctma_time_reporting")<-NULL
  original$originalTime<-NULL
  original$estimates_original_time_scale<-NULL
  object$summary<-original
  object
}
