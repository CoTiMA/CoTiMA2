test_that('table warnings distinguish posterior failure from placement failure',{
 s<-list(parmatrices=data.frame(matrix='DRIFT',row=1,col=1,Mean=-.3,sd=.02),npars=1)
 fit<-structure(list(summary=s,optim=list(converged=FALSE),
   sample=list(sampler='nuts',converged=TRUE),args=list(resolved=list(optimize=FALSE,priors=TRUE))),class='ctJuliaFit')
 before<-serialize(fit,NULL)
 x<-ctmaTables(fit,parameters='drift',statistics='mean',format='none')
 expect_true(any(grepl('placement optimizer',x$notes,fixed=TRUE)))
 expect_false(any(grepl('WARNING:',x$notes,fixed=TRUE)))
 expect_identical(serialize(fit,NULL),before)
 fit$sample$converged<-FALSE
 x<-ctmaTables(fit,parameters='drift',statistics='mean',format='none')
 expect_equal(sum(grepl('posterior sampling explicitly reports',x$notes,fixed=TRUE)),1)
 fit$sample$converged<-TRUE
 fit$summary$popmeans<-data.frame(Rhat=1.04)
 x<-ctmaTables(fit,parameters='drift',statistics='mean',format='none')
 expect_true(any(grepl('transformed-summary Rhat',x$notes,fixed=TRUE)))
 fit$sample<-NULL;fit$args$resolved$optimize<-TRUE
 x<-ctmaTables(fit,parameters='drift',statistics='mean',format='none')
 expect_true(any(grepl('fit explicitly reports non-convergence',x$notes,fixed=TRUE)))
})
test_that('native cached resultsSummary is read without calling summary',{
 s<-list(parmatrices=data.frame(matrix='DRIFT',row=1,col=1,Mean=-.3,sd=.02),npars=1)
 fit<-structure(list(resultsSummary=s,sample=list(converged=TRUE,sampler='saem')),
   class=c('ctJuliaFit','ctFit'))
 x<-ctmaTables(fit,parameters='drift',statistics='mean',format='none')
 expect_equal(x$data$Mean,-.3)
})
