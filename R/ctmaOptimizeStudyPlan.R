# Pure planning and reparameterization helpers. Never fit a model or consume RNG.
.ctmaOptimizeStop <- function(...) stop(...,call.=FALSE)
.ctmaOptimizeHash <- function(x) {
  path<-tempfile();on.exit(unlink(path));saveRDS(x,path,version=2,compress=FALSE)
  unname(tools::md5sum(path))
}
.ctmaOptimizeScalar <- function(x,name,integer=FALSE,minimum=0) {
  if(!is.numeric(x)||length(x)!=1L||!is.finite(x)||x<minimum||
      (integer && x!=floor(x))) .ctmaOptimizeStop(name," must be one finite ",if(integer)'integer' else 'number'," >= ",minimum,".")
  x
}
.ctmaOptimizeNative <- function(fit) {
  if(inherits(fit,c('ctJuliaFit','ctStanFit','ctmaJuliaFit')))return(fit)
  if(inherits(fit$studyFitList,c('ctJuliaFit','ctStanFit','ctmaJuliaFit')))return(fit$studyFitList)
  .ctmaOptimizeStop('Supply one pooled fit; a ctmaInit studyFitList is not a pooled model.')
}
.ctmaOptimizeLayout <- function(fit) {
  native<-.ctmaOptimizeNative(fit)
  julia<-inherits(native,'ctJuliaFit')
  if(!julia && !inherits(native,'ctStanFit'))
    .ctmaOptimizeStop('Systematic searches currently support stan and native ctsemjulia ML; bundled julia remains available with gambling=TRUE.')
  raw<-if(julia)native$estimate$raw else native$stanfit$rawest
  if(!is.numeric(raw)||!length(raw)||any(!is.finite(raw)))
    .ctmaOptimizeStop('The pooled fit must retain finite raw ML point estimates.')
  keys<-rep(NA_character_,length(raw)); ti<-data.frame(index=integer(),base=character(),predictor=character())
  if(julia) {
    spec<-native$model_spec; p<-spec$parameter_table
    if(!is.data.frame(p)||!all(c('parnumber','param')%in%names(p))) .ctmaOptimizeStop('Native Julia parameter layout is unavailable.')
    good<-!is.na(p$parnumber)&p$parnumber>0&!is.na(p$param)
    for(i in which(good))keys[p$parnumber[i]]<-paste0('base:',p$param[i])
    levels<-spec$laplace$levels
    if(length(levels)>1L).ctmaOptimizeStop('Systematic study exclusions currently support fixed study effects, not hierarchical study random effects.')
    if(length(levels)) {
      lv<-levels[[1L]];k<-length(lv$re_index)
      if(length(lv$load_index)||(!is.null(lv$rank)&&lv$rank!=k)) .ctmaOptimizeStop('Reduced-rank population covariance is not supported by systematic start mapping.')
      labels<-sub('^base:','',keys[lv$re_index])
      for(i in seq_len(k))if(lv$sd_index[i]>0)keys[lv$sd_index[i]]<-paste0('sd:',labels[i])
      pair<-which(lower.tri(diag(k)),arr.ind=TRUE);pair<-pair[order(pair[,2],pair[,1]),,drop=FALSE]
      for(i in seq_len(nrow(pair)))if(lv$cor_index[i]>0)
        keys[lv$cor_index[i]]<-paste0('cor:',labels[pair[i,1]],'|',labels[pair[i,2]])
    }
    e<-spec$ti_effects
    if(is.data.frame(e)&&nrow(e))ti<-data.frame(index=e$coefficient,
      base=sub('^base:','',keys[e$parameter]),predictor=spec$TIpredNames[e$predictor],stringsAsFactors=FALSE)
  } else {
    d<-native$standata;p<-native$ctstanmodel$pars
    if(!is.data.frame(p)||nrow(p)!=nrow(d$matsetup)) .ctmaOptimizeStop('Stan generated parameter layout is unavailable.')
    idx<-d$matsetup[,3L];good<-idx>0&!is.na(p$param)
    for(i in which(good)) keys[idx[i]]<-paste0('base:',p$param[i])
    n<-as.integer(d$nparams);k<-as.integer(d$nindvarying)
    labels<-sub('^base:','',keys[as.integer(d$indvaryingindex)])
    if(k) {
      if(d$nindvaryingoffdiagonals!=k*(k-1)/2) .ctmaOptimizeStop('Nonstandard Stan population covariance layout is not supported.')
      keys[n+seq_len(k)]<-paste0('sd:',labels)
      pair<-which(lower.tri(diag(k)),arr.ind=TRUE);pair<-pair[order(pair[,2],pair[,1]),,drop=FALSE]
      for(i in seq_len(nrow(pair)))keys[n+k+i]<-paste0('cor:',labels[pair[i,1]],'|',labels[pair[i,2]])
    }
    e<-d$TIPREDEFFECTsetup
    if(length(e)&&any(e>0)) {
      cells<-which(e>0,arr.ind=TRUE);offset<-n+k+as.integer(d$nindvaryingoffdiagonals)
      ti<-data.frame(index=offset+e[cells],base=sub('^base:','',keys[cells[,1L]]),
        predictor=native$ctstanmodelbase$TIpredNames[cells[,2L]],stringsAsFactors=FALSE)
    }
  }
  if(nrow(ti)) for(i in seq_len(nrow(ti))) {
    key<-paste0('ti:',ti$base[i],'|',ti$predictor[i]);at<-ti$index[i]
    if(!is.na(keys[at])&&!identical(keys[at],key)) .ctmaOptimizeStop('Tied or inconsistent TI coefficient layouts are not supported.')
    keys[at]<-key
  }
  if(anyNA(keys)||anyDuplicated(keys)) .ctmaOptimizeStop('Raw coordinates do not have a complete, unique scientific parameter map.')
  list(raw=stats::setNames(as.numeric(raw),keys),keys=keys,ti=ti)
}
.ctmaOptimizeContext <- function(fit,init=NULL,backend=NULL) {
  if(is.null(fit)) .ctmaOptimizeStop('gambling=FALSE requires ctmaFitFit.')
  native<-.ctmaOptimizeNative(fit)
  actual<-if(inherits(native,'ctJuliaFit'))'ctsemjulia' else if(inherits(native,'ctStanFit'))'stan' else 'julia'
  if(is.null(backend))backend<-actual
  if(!identical(tolower(backend),actual)) .ctmaOptimizeStop('Systematic searches retain the fitted backend; fit a separate baseline to change engines.')
  a<-fit$argumentList
  if(is.null(a))a<-native$argumentList
  resolved<-native$args$resolved
  optimize<-if(!is.null(a$optimize))a$optimize else resolved$optimize
  priors<-if(!is.null(a$priors))a$priors else resolved$priors
  if(!isTRUE(optimize)||!isFALSE(priors)) .ctmaOptimizeStop('gambling=FALSE currently supports ML only (optimize=TRUE, priors=FALSE); rebasing MAP/Bayesian priors is not automatic.')
  model<-fit$ctModel
  if(is.null(model))model<-native$source_object$ctModel
  if(is.null(model))model<-if(actual=='ctsemjulia')native$modelbase else native$ctstanmodelbase
  data<-fit$data
  if(!is.data.frame(data))data<-native$source_object$data
  if(!is.data.frame(data)) .ctmaOptimizeStop('Retained long data are required; supply the complete pooled CoTiMA fit.')
  if(!is.list(model)||!is.data.frame(model$pars)) .ctmaOptimizeStop('The original model specification was not retained.')
  if(length(model$groupIDnames)) .ctmaOptimizeStop('Hierarchical study random effects are not yet supported by this fixed-study sensitivity workflow.')
  coding<-if(is.null(a$studyCoding))if(isTRUE(a$WEC))'weighted' else 'dummy' else a$studyCoding
  if(!identical(coding,'dummy')) .ctmaOptimizeStop('Reference/exclusion planning currently requires dummy study coding; weighted coding is not silently recoded.')
  studies<-fit$studyList
  if(is.null(studies))studies<-native$source_object$studyList
  if(is.null(studies)&&!is.null(init))studies<-init$studyList
  if(!length(studies)||length(studies)<2L) .ctmaOptimizeStop('Original study metadata for at least two studies are required.')
  ids<-vapply(studies,function(s)as.character(s$originalStudyNo),character(1))
  if(anyNA(ids)||anyDuplicated(ids)) .ctmaOptimizeStop('Original study identifiers must be unique and nonmissing.')
  columns<-paste0('TI',seq_len(length(ids)-1L))
  if(!all(c('id','time',columns)%in%names(data))) .ctmaOptimizeStop('Cannot verify the CoTiMA study-indicator design in the retained long data.')
  x<-as.matrix(data[,columns,drop=FALSE]);storage.mode(x)<-'double'
  if(any(!is.finite(x))) .ctmaOptimizeStop('Study indicators contain missing or nonfinite values.')
  active<-do.call(cbind,lapply(seq_len(ncol(x)),function(j){u<-sort(unique(x[,j]));if(length(u)!=2L).ctmaOptimizeStop('Each study dummy must have two observed levels.');x[,j]==u[2L]}))
  if(any(rowSums(active)>1L)) .ctmaOptimizeStop('Study indicators do not define mutually exclusive groups.')
  group<-ifelse(rowSums(active)==0,length(ids),max.col(active,ties.method='first'))
  if(any(vapply(split(group,data$id),function(g)length(unique(g)),integer(1))!=1L)) .ctmaOptimizeStop('Study membership changes within a participant.')
  design<-do.call(rbind,lapply(seq_along(ids),function(i){z<-unique(x[group==i,,drop=FALSE]);if(nrow(z)!=1L).ctmaOptimizeStop('Study design is not constant within a study.');as.numeric(z)}))
  design<-cbind(Intercept=1,design);colnames(design)<-c('Intercept',columns)
  if(qr(design)$rank!=length(ids)) .ctmaOptimizeStop('The full study-contrast design is rank deficient.')
  layout<-.ctmaOptimizeLayout(fit)
  # Only ordinary Boolean study-effect flags can safely be rewritten here.
  for(nm in paste0(columns,'_effect')) {
    z<-model$pars[[nm]]
    if(is.null(z)||any(!is.na(z)&!as.character(z)%in%c('TRUE','FALSE'))) .ctmaOptimizeStop('Named, tied or fixed study effects require an explicit model-specific transformation.')
  }
  list(fit=fit,init=init,native=native,backend=actual,args=a,resolved=resolved,model=model,data=data,
    studies=studies,ids=ids,group=group,columns=columns,design=design,layout=layout,
    counts=vapply(seq_along(ids),function(i)length(unique(data$id[group==i])),integer(1)),
    scaleTime=if(is.null(a$scaleTime))1 else a$scaleTime,
    fingerprint=.ctmaOptimizeHash(list(model=model,data=data,ids=ids,backend=actual,
      raw=layout$raw,resolved=resolved,scales=a[c('scaleTime','scaleTI','studyCoding')])))
}
.ctmaOptimizePlan <- function(context,candidates=NULL,exclusion='none',blocks=NULL,
 reference='current',references=NULL,maxCases=20L) {
  .ctmaOptimizeScalar(maxCases,'maxCases',TRUE,1)
  if(!is.character(exclusion)||!length(exclusion)||anyNA(exclusion)||any(!exclusion%in%c('none','single','blocks','all'))||('none'%in%exclusion&&length(exclusion)>1)) .ctmaOptimizeStop('exclusion must select none, single, blocks and/or all.')
  ids<-context$ids
  automatic<-.ctmaOptimizeAutoCandidates(context,candidates,references)
  if(!is.null(automatic)) {
    candidates<-automatic$candidateStudies;references<-automatic$referenceCandidates
    if(automatic$autoCandidates&&!length(candidates)) {
      exclusion<-setdiff(exclusion,c('single','all'))
      if(!length(exclusion))exclusion<-'none'
    }
  }
  validate<-function(x,label){z<-as.character(x);if(!length(z)||anyNA(z)||anyDuplicated(z)||any(!z%in%ids)) .ctmaOptimizeStop('Invalid ',label,': use unique original study IDs.');z}
  if(length(candidates))candidates<-validate(candidates,'candidateStudies')
  if(any(c('single','all')%in%exclusion)&&!length(candidates)) .ctmaOptimizeStop('Explicit candidateStudies are required for exclusions; diagnostic flags do not automatically delete data.')
  if('blocks'%in%exclusion) {
    if(!is.list(blocks)||!length(blocks)) .ctmaOptimizeStop('studyBlocks must be a nonempty list of original study IDs.')
    blocks<-lapply(blocks,validate,label='studyBlocks')
  }
  choose<-function(ref,keep){
    if(length(ref)!=1L||is.na(ref)) .ctmaOptimizeStop('referenceStudy must be current, largestN, or one original study ID.')
    ref<-as.character(ref)
    if(ref=='current')ref<-tail(ids,1L)
    if(ref=='largestN'||!ref%in%keep) {
      if(!ref%in%c('largestN',ids)) .ctmaOptimizeStop('Unknown referenceStudy.')
      ref<-keep[which.max(context$counts[match(keep,ids)])]
    }
    if(!ref%in%ids) .ctmaOptimizeStop('Unknown referenceStudy.')
    ref
  }
  rows<-list();seen<-character()
  add<-function(label,excluded=character(),ref=reference){
    keep<-ids[!ids%in%excluded];if(length(keep)<2L) .ctmaOptimizeStop('A planned case leaves fewer than two studies.')
    r<-choose(ref,keep);key<-paste(paste(sort(keep),collapse='|'),r,sep='::')
    if(key%in%seen)return(invisible(NULL));seen<<-c(seen,key)
    rows[[length(rows)+1L]]<<-list(case=label,included=keep,excluded=excluded,reference=r,
      fullData=!length(excluded),nsubjects=sum(context$counts[match(keep,ids)]))
  }
  add('full_current',ref='current')
  add('full_reference',ref=reference)
  if(length(references))for(r in references)add(paste0('reference_',r),ref=r)
  if('single'%in%exclusion)for(id in candidates)add(paste0('omit_',id),id)
  if('blocks'%in%exclusion)for(i in seq_along(blocks))add(paste0('block_',i),blocks[[i]])
  if('all'%in%exclusion)add('omit_all_candidates',candidates)
  if(length(rows)>maxCases) .ctmaOptimizeStop('Study plan exceeds maxCases; reduce the explicit plan.')
  attr(rows,'automaticSelection')<-automatic
  rows
}
.ctmaOptimizeCase <- function(context,case) {
  ids<-context$ids;keep<-case$included
  order<-c(keep[keep!=case$reference],case$reference)
  selected<-context$group%in%match(keep,ids)
  data<-context$data[selected,,drop=FALSE]
  groups<-ids[context$group[selected]]
  cols<-paste0('TI',seq_len(length(order)-1L))
  extras<-setdiff(context$model$TIpredNames,context$columns)
  data[,context$columns]<-NULL
  newx<-sapply(order[-length(order)],function(id)as.numeric(groups==id))
  if(is.null(dim(newx)))newx<-matrix(newx,ncol=1)
  if(isTRUE(context$args$scaleTI)) {
    subject<-!duplicated(data$id)
    means<-colMeans(newx[subject,,drop=FALSE]);sds<-apply(newx[subject,,drop=FALSE],2,stats::sd)
    newx<-sweep(sweep(newx,2,means,'-'),2,sds,'/')
  }
  for(i in seq_along(cols))data[[cols[i]]]<-newx[,i]
  newdesign<-do.call(rbind,lapply(order,function(id)as.numeric(unique(newx[groups==id,,drop=FALSE]))))
  newdesign<-cbind(Intercept=1,newdesign);colnames(newdesign)<-c('Intercept',cols)
  model<-context$model;p<-model$pars
  oldFlags<-do.call(cbind,lapply(paste0(context$columns,'_effect'),function(nm)as.logical(p[[nm]])))
  oldFlags[is.na(oldFlags)]<-FALSE
  rebased<-case$reference!=tail(ids,1L)
  if(rebased&&any(rowSums(oldFlags)>0&rowSums(oldFlags)<ncol(oldFlags)))
    .ctmaOptimizeStop('Changing reference with partially constrained study effects would change the model; choose current reference or an explicit equivalent model.')
  p[,paste0(context$columns,'_effect')]<-NULL
  for(i in seq_along(cols))p[[paste0(cols[i],'_effect')]]<-if(rebased)rowSums(oldFlags)>0 else oldFlags[,match(order[i],ids)]
  model$pars<-p;model$TIpredNames<-c(cols,extras);model$n.TIpred<-length(model$TIpredNames)
  # Baseline plus study contrasts are linear on ctsem's RAW predictor scale.
  # Solving the full study design also handles centered/scaled dummy columns.
  raw<-context$layout$raw
  bases<-unique(context$layout$ti$base[context$layout$ti$predictor%in%context$columns])
  coefficient_values<-list();newbase<-raw[startsWith(names(raw),'base:')]
  for(base in bases) {
    old<-c(raw[[paste0('base:',base)]],vapply(context$columns,function(col){k<-paste0('ti:',base,'|',col);if(k%in%names(raw))raw[[k]] else 0},numeric(1)))
    eta<-as.numeric(context$design[match(order,ids),,drop=FALSE]%*%old)
    if(rebased)beta<-as.numeric(solve(newdesign,eta)) else {
      # Under current reference, drop the removed contrasts without altering
      # constraints. A design solve is valid only for fully free families.
      beta<-as.numeric(solve(newdesign,eta))
    }
    if(max(abs(newdesign%*%beta-eta))>1e-8) .ctmaOptimizeStop('Reference transformation did not preserve study-specific raw predictors.')
    newbase[[paste0('base:',base)]]<-beta[1]
    for(i in seq_along(cols))coefficient_values[[paste0('ti:',base,'|',cols[i])]]<-beta[i+1L]
  }
  # Build the target map explicitly: no recycled logical masks or positional TI mapping.
  targets<-as.list(raw[!startsWith(names(raw),'ti:')])
  for(n in names(newbase))targets[[n]]<-unname(newbase[[n]])
  for(i in seq_len(nrow(context$layout$ti)))if(context$layout$ti$predictor[i]%in%extras){
    k<-names(raw)[context$layout$ti$index[i]];targets[[k]]<-unname(raw[[k]])
  }
  for(n in names(coefficient_values))targets[[n]]<-coefficient_values[[n]]
  list(case=case,model=model,data=data,targets=targets,studies=context$studies[match(order,ids)],
    dataFingerprint=.ctmaOptimizeHash(list(data=data,ids=order)),backend=context$backend,
    resolved=context$resolved,scaleTime=context$scaleTime)
}
