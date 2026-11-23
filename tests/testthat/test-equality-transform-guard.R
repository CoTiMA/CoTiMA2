test_that('equality rejects different scientific-scale transformations', {
 m<-list(pars=data.frame(transform=c('-exp(param)','param','param','-exp(param)')))
 expect_error(.ctmaEqualCheckTransforms(m,1:4),'different transformations')
 expect_true(.ctmaEqualCheckTransforms(m,2:3))
 expect_true(.ctmaEqualCheckTransforms(m,c(1,4)))
 m$pars$transform[3]<-' param '
 expect_true(.ctmaEqualCheckTransforms(m,2:3))
 m$pars$multiplier<-c(1,1,2,1)
 expect_error(.ctmaEqualCheckTransforms(m,2:3),'different transformations')
})
