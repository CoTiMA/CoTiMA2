.long_named_stationary_model <- function(k=0L) {
 m<-ctsem::ctModel(type='ct',n.latent=2,n.manifest=2,LAMBDA=diag(2),
  DRIFT=matrix(c('V1toV1','V1toV2','V2toV1','V2toV2'),2),
  DIFFUSION=matrix(c('diff_eta1','diff_eta2_eta1',0,'diff_eta2'),2),
  T0VAR='auto',T0MEANS=matrix(0,2,1),CINT=matrix(0,2,1),
  MANIFESTMEANS=matrix(c('mm_V1|param|TRUE|1','mm_V2|param|TRUE|1'),2,1),
  MANIFESTVAR=diag(.1,2),n.TIpred=k,TIpredNames=if(k)paste0('TI',seq_len(k))else NULL)
 m$covmattransform<-'z'
 for(v in paste0(m$TIpredNames,'_effect'))if(v %in% names(m$pars))
  m$pars[[v]]<-m$pars$matrix=='DIFFUSION' & !is.na(m$pars$param)
 m
}

test_that('realistic parameter labels do not exceed R environment-name limits', {
 for(k in c(0L,2L,44L)) {
  m<-.long_named_stationary_model(k)
  cols<-if(k)paste0(m$TIpredNames,'_effect')else NULL
  X<-if(k)rbind(rep(0,k),diag(k))else NULL
  z<-suppressMessages(.ctmaDeterminedDiffusionModel(m,TRUE,cols,X))
  expect_lt(max(nchar(z$pars$param,type='bytes'),na.rm=TRUE),10000L)
  # Exercise the exact private ctsem setup step named in the user traceback.
  prepare<-getFromNamespace('ctModelTransformsToNum','ctsem')
  transformed<-prepare(z)
  # Exercise R's actual name-size limit under either ctsem API version.
  registry<-new.env(parent=emptyenv())
  expect_no_error(for(label in stats::na.omit(transformed$pars$param))
    assign(paste0("p",label),TRUE,envir=registry))
  ns<-asNamespace("ctsem")
  if(exists(".ctModelMatSetup",envir=ns,inherits=FALSE)) {
    setup<-get(".ctModelMatSetup",envir=ns,inherits=FALSE)
    expect_no_error(result<-setup(transformed))
    expect_true(is.list(result))
  } else {
    # Published ctsem has a different writer API; no native API is claimed.
    expect_true(all(c("matrix","param","transform") %in% names(transformed$pars)))
    expect_gt(length(ls(registry)),0L)
  }
  expect_identical(z$covmattransform,'z')
 }
})

test_that('compact z expressions retain the stationary Lyapunov identity', {
 m<-suppressMessages(.ctmaDeterminedDiffusionModel(.long_named_stationary_model(),TRUE))
 rows<-m$pars[m$pars$matrix=='PARS',,drop=FALSE]
 for(seed in 1:12) {
  set.seed(seed);raw<-stats::setNames(rnorm(nrow(rows),sd=.35),rows$param)
  values<-as.list(raw)
  for(i in seq_len(nrow(rows))) values[[rows$param[i]]]<-
   eval(parse(text=rows$transform[i]),list(param=raw[[rows$param[i]]],log1p_exp=function(x)log1p(exp(x))))
  calc<-function(family) {
   p<-m$pars[m$pars$matrix==family,,drop=FALSE];B<-matrix(0,2,2)
   for(i in seq_len(nrow(p)))B[p$row[i],p$col[i]]<-
    if(is.na(p$param[i]))p$value[i]else eval(parse(text=p$param[i]),values)
   B
  }
  cov<-function(B){s<-diag(B);R<-diag(2);R[1,2]<-R[2,1]<-tanh(B[2,1]);diag(s)%*%R%*%diag(s)}
  A<-calc('DRIFT');Q<-cov(calc('DIFFUSION'));S<-cov(calc('T0VAR'))
  expect_equal(A%*%S+S%*%t(A)+Q,matrix(0,2,2),tolerance=1e-9)
  expect_gt(min(eigen(Q,symmetric=TRUE)$values),0)
  expect_gt(min(eigen(S,symmetric=TRUE)$values),0)
 }
})
