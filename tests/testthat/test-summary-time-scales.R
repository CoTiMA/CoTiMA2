scale_table <- function() data.frame(matrix=c('DRIFT','DIFFUSIONcov','T0cov','DIFFUSION','PARS'),
 row=1:5,col=1,Mean=c(-1.2,2.4,1,.5,8),sd=c(.12,.24,.1,.05,.8),
 `2.5%`=c(-1.4,2.1,.8,.3,7),`50%`=c(-1.2,2.4,1,.5,8),`97.5%`=c(-1,2.7,1.2,.7,9),
 mode=c(-1.21,2.41,1,.5,8),Variance=c(.0144,.0576,.01,.0025,.64),
 z=1:5,n_eff=101:105,Rhat=rep(1.01,5),check.names=FALSE)
test_that('time conversion acts on transformed rates, not indices or diagnostics', {
 d<-scale_table();f<-list(argumentList=list(scaleTime=1/12));s<-list(parmatrices=d,loglik=-123,npars=5)
 x<-.ctmaAddSummaryTimeScales(s,f,interval=c(1,5))
 expect_identical(x$parmatrices,d);expect_identical(x$estimates_fitted_time_scale,d)
 q<-x$estimates_original_time_scale;expect_equal(q$matrix,c('DRIFT','DIFFUSIONcov'))
 for(col in c('Mean','sd','2.5%','50%','97.5%','mode'))expect_equal(q[[col]],d[[col]][1:2]/12)
 expect_equal(q$Variance,d$Variance[1:2]/144)
 for(col in c('row','col','z','n_eff','Rhat'))expect_identical(q[[col]],d[[col]][1:2])
 expect_equal(x$time_scale_info$interval_original,c(12,60))
 expect_equal(x$loglik,-123);expect_equal(x$npars,5)
 expect_identical(.ctmaAddSummaryTimeScales(x,f),x)
 expect_error(.ctmaAddSummaryTimeScales(x,list(argumentList=list(scaleTime=2))),'different time scale')
})
test_that('initialization and conventional summaries retain legacy original tables', {
 d<-matrix(c('source','-1.2','.12'),1,dimnames=list('Study No 7',c('source','a','SE')))
 orig<-matrix(c('source','-.10001','.01001'),1,dimnames=dimnames(d))
 f<-list(argumentList=list(scaleTime=1/12),summary=list(scaleTime=1/12))
 s<-list(estimates=d,drift_estimates_original_time_scale=orig,minus2ll=10)
 x<-.ctmaAddSummaryTimeScales(s,f,'initialization')
 expect_identical(x$estimates,d);expect_identical(x$estimates_original_time_scale,orig)
 expect_identical(x$drift_estimates_original_time_scale,orig)
 y<-.ctmaAddSummaryTimeScales(list(estimates=scale_table(),estimates_original_time_scale=orig),f,'conventional')
 expect_identical(y$estimates_original_time_scale,orig) # Do not round/reconstruct a stored original table.
 expect_identical(f$summary,list(scaleTime=1/12))
})
test_that('unit scale, absent tables and invalid metadata are explicit', {
 expect_identical(.ctmaSummaryTimeFactor(list(summary='setup only')),1)
 expect_identical(.ctmaSummaryTimeFactor(list(argumentList=list(scaleTimeAuto=TRUE))),1)
 for(s in list(0,-1,NA_real_,Inf,c(1,2),'1'))expect_error(.ctmaSummaryTimeFactor(list(argumentList=list(scaleTime=s))),'positive number')
 expect_error(.ctmaSummaryTimeFactor(list(argumentList=list(scaleTime=.5),summary=list(scaledTime=2))),'Conflicting')
 x<-.ctmaAddSummaryTimeScales(list(npars=3),list())
 expect_null(x$estimates_original_time_scale);expect_null(x$estimates_fitted_time_scale)
 expect_null(.ctmaOriginalRateTable(data.frame(matrix='DIFFUSION',Mean=2),.1))
 expect_identical(.ctmaOriginalRateTable(scale_table(),1),scale_table()[1:2,])
})
test_that('printing labels both scales and does not expose duplicate aliases', {
 x<-.ctmaAddSummaryTimeScales(list(parmatrices=scale_table()),list(argumentList=list(scaleTime=.25)),interval=1)
 out<-capture.output(print(x))
 expect_match(paste(out,collapse='\n'),'original time scale')
 expect_match(paste(out,collapse='\n'),'original \\[4\\]')
 expect_false(any(grepl('^\\$estimates_fitted_time_scale',out)))
 y<-.ctmaAddSummaryTimeScales(list(estimates=scale_table()),list(),'conventional')
 expect_true(any(grepl('^\\$estimates_fitted_time_scale',capture.output(print(y)))))
})

test_that('unavailable uncertainty remains unavailable after conversion', {
  d<-scale_table();d$sd[1]<-NA_real_;d$`2.5%`[1]<-NA_real_
  q<-.ctmaOriginalRateTable(d,2)
  expect_true(is.na(q$sd[1]));expect_true(is.na(q$`2.5%`[1]))
  expect_equal(q$Mean,d$Mean[1:2]*2)
})
