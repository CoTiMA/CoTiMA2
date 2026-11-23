test_that("explicit FALSE removes only initial-mean individual variation", {
 m<-list(pars=data.frame(matrix=c("T0MEANS","MANIFESTMEANS","CINT","T0VAR"),
   param=c("t0","ri","ci","v"),value=NA_real_,indvarying=TRUE),n.latent=1L)
 out<-.ctmaApplyT0MeanVariation(m,FALSE)
 expect_identical(out$pars$indvarying,c(FALSE,TRUE,TRUE,TRUE))
 expect_identical(out$pars$param,m$pars$param)
 expect_identical(.ctmaApplyT0MeanVariation(m,NULL),m)
 expect_identical(.ctmaApplyT0MeanVariation(m,TRUE),m)
 for(x in list(NA,0,c(TRUE,FALSE),"FALSE"))expect_error(.ctmaApplyT0MeanVariation(m,x),"indVaryingT0")
})

test_that("both public constructors apply the explicit choice before determining diffusion", {
 for(fn in list(ctmaInit,ctmaFit)){
  text<-paste(deparse(body(fn)),collapse="\n")
  expect_true(regexpr(".ctmaApplyT0MeanVariation",text,fixed=TRUE)[1] <
    regexpr(".ctmaDeterminedDiffusionModel",text,fixed=TRUE)[1])
 }
})

test_that("ctmaInit fit=FALSE honours the explicit initial-mean switch", {
 e<-new.env();utils::data("CoTiMAstudyList_3",package="CoTiMA2",envir=e)
 s<-e$CoTiMAstudyList_3
 for(nm in names(s))if(is.list(s[[nm]])&&!is.data.frame(s[[nm]])&&!inherits(s[[nm]],"Workbook")&&length(s[[nm]])==4L)s[[nm]]<-s[[nm]][1L]
 s$n.studies<-1L
 p<-suppressMessages(suppressWarnings(ctmaInit(primaryStudies=s,n.latent=2,
   indVarying="MANIFEST",indVaryingT0=FALSE,fit=FALSE,silentOverwrite=TRUE,
   activeDirectory=paste0(tempdir(),"/"))))
 expect_false(any(p$ctModel$pars$indvarying[p$ctModel$pars$matrix=="T0MEANS"]))
 expect_true(all(p$ctModel$pars$indvarying[p$ctModel$pars$matrix=="MANIFESTMEANS"]))
})

test_that("ctmaFit fit=FALSE retains manifest random effects with T0 variation disabled", {
 e<-new.env();utils::data("CoTiMAstudyList_3",package="CoTiMA2",envir=e)
 utils::data("CoTiMAInitFit_3",package="CoTiMA2",envir=e)
 fresh<-suppressMessages(suppressWarnings(ctmaInit(primaryStudies=e$CoTiMAstudyList_3,
  n.latent=2,fit=FALSE,silentOverwrite=TRUE,activeDirectory=paste0(tempdir(),"/"))))
 init<-e$CoTiMAInitFit_3;init$emprawList<-fresh$data
 init$studyFitList<-lapply(fresh$data,function(d)list(empraw=d))
 p<-suppressMessages(suppressWarnings(ctmaFit(ctmaInitFit=init,fit=FALSE,
  indVarying="MANIFEST",indVaryingT0=FALSE,T0means="auto",allInvModel=TRUE,
  activeDirectory=paste0(tempdir(),"/"))))
 expect_false(any(p$ctModel$pars$indvarying[p$ctModel$pars$matrix=="T0MEANS"]))
 expect_true(all(p$ctModel$pars$indvarying[p$ctModel$pars$matrix=="MANIFESTMEANS"]))
})
