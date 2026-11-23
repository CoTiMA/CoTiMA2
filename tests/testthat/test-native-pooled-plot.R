.native_plot_fixture <- function(scale=.25) {
 A<-matrix(c(-.4,.12,-.08,-.6),2)
 model<-list(n.latent=2,pars=data.frame(matrix=rep('DRIFT',4),param=c('V1toV1','V1toV2','V2toV1','V2toV2'),indvarying=FALSE,TI1_effect=FALSE))
 native<-structure(list(model_spec=list(model=model),ctsemjulia_diagnostics=list(population_matrices=list(DRIFT=A))),class='ctJuliaFit')
 out<-structure(list(studyFitList=native,ctModel=model,argumentList=list(scaleTime=scale,activeDirectory=tempdir()),studyList=list(list(originalStudyNo=1,delta_t=4,sampleSize=100))),class=c('CoTiMAPooledJuliaFit','CoTiMAFit'))
 list(fit=out,A=A)
}
test_that('native plotting restores dimensions, matrices and original time units without touching fit',{
 z<-.native_plot_fixture(); before<-z$fit
 x<-suppressMessages(.ctmaPlotCompatible(z$fit))
 expect_identical(x$n.latent,2L);expect_identical(x$n.studies,1L)
 expect_equal(matrix(x$modelResults$DRIFT,2),z$A)
 expect_equal(matrix(x$modelResults$DRIFToriginal_time_scale,2),z$A*.25)
 expect_identical(x$summary$scaledTime,.25)
 expect_identical(z$fit,before)
})
test_that('native plotting does not silently replace predictor or random drift',{
 z<-.native_plot_fixture();z$fit$ctModel$pars$TI1_effect[2]<-TRUE
 expect_error(.ctmaPlotCompatible(z$fit),'conditional curves')
 z<-.native_plot_fixture();z$fit$ctModel$pars$indvarying[1]<-TRUE
 expect_error(.ctmaPlotCompatible(z$fit),'conditional curves')
})
test_that('determined-diffusion auxiliary drift dependencies are checked',{
 z<-.native_plot_fixture(); p<-z$fit$ctModel$pars
 p$matrix<-'PARS'
 q<-p;q$matrix<-'DRIFT';q$param<-paste0('1*',p$param);q$TI1_effect<-FALSE
 z$fit$ctModel$pars<-rbind(p,q)
 expect_silent(suppressMessages(.ctmaPlotCompatible(z$fit)))
 z$fit$ctModel$pars$TI1_effect[1]<-TRUE
 expect_error(.ctmaPlotCompatible(z$fit),'conditional curves')
})
