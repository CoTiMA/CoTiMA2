test_that('removed inert Stan defaults are filtered without changing other controls',{
 checker<-function(controls,backend)stop('Unrecognised optimcontrol name(s): ',names(controls))
 x<-list(optimcontrol=list(parstepsAutoModel=FALSE,groupFreeThreshold=.5,tol=1e-8))
 y<-.ctmaRetiredOptimControls(x,checker)
 expect_identical(y$optimcontrol,list(tol=1e-8));expect_length(x$optimcontrol,3L)
 old<-function(controls,backend)invisible(TRUE)
 expect_identical(.ctmaRetiredOptimControls(x,old),x)
 x$optimcontrol$parstepsAutoModel<-TRUE
 expect_error(.ctmaRetiredOptimControls(x,checker),'not silently discarded')
 x$optimcontrol$parstepsAutoModel<-FALSE;x$optimcontrol$groupFreeThreshold<-.7
 expect_error(.ctmaRetiredOptimControls(x,checker),'not silently discarded')
})
