library(testthat)
test_that('native controls forward iter without changing target or precedence',{
 # Exercise control forwarding independently of the optional engine installation.
 local_mocked_bindings(.ctmaRequireCtsemJuliaLaplace=function()list(available=TRUE),
   .package='CoTiMA2')
 x<-.ctmaResolveFitBackend('ctsemjulia',list(sampleControl=list(iter=40L,draws=12L,warmup=20L)),estimation='Bayes',workflow='ctmaInit')
 expect_identical(x$backend_args$sampleControl$iter,40L)
 expect_identical(x$backend_args$sampleControl$draws,12L)
 expect_identical(x$backend_args$sampleControl$warmup,20L)
 expect_identical(.ctmaCtsemJuliaRoute(x$backend_args,'Bayes',FALSE),'auto')
})
test_that('sampled compact status never borrows placement status',{
 for(optimizer in c(TRUE,FALSE)) for(posterior in c(TRUE,FALSE)) {
  x<-list(optim=list(converged=optimizer),sample=list(converged=posterior,sampler='nuts'))
  z<-.ctmaAnnotateCtsemJuliaResult(x)$ctsemjulia_diagnostics
  expect_identical(z$converged,posterior);expect_identical(z$optimizer_converged,optimizer)
  expect_identical(z$convergence_source,'posterior sampling')
 }
 z<-.ctmaAnnotateCtsemJuliaResult(list(sample=list(sampler='saem'),optim=list(converged=TRUE)))
 expect_true(is.na(z$ctsemjulia_diagnostics$converged))
})
test_that('pooled and single-primary storage retain different diagnostic shapes',{
 native<-structure(list(optim=list(converged=FALSE),sample=list(converged=TRUE,sampler='nuts')),class='ctJuliaFit')
 pooled<-structure(list(studyFitList=native,summary=list()),class='CoTiMAFit')
 z<-ctmaFitDiagnostics(pooled,print=FALSE);expect_true(z$sampling);expect_true(z$posterior$reported_converged)
 expect_output(summary(pooled),'Pooled model')
 one<-list(studyFitList=list(native),primaryStudyList=list(studyNumbers=17))
 expect_named(ctmaFitDiagnostics(one,print=FALSE),'Study 17')
})
test_that('stored transformed diagnostics remain visible and SMC does not need Rhat',{
 x<-list(sample=list(sampler='saem',converged=TRUE,rhat=c(a=1.001),ess=c(a=500)),
  resultsSummary=list(rawpopcorr=data.frame(n_eff=178,Rhat=1.04,row.names='rho')))
 z<-ctmaFitDiagnostics(x,print=FALSE)$posterior
 expect_equal(z$rhat,1.04);expect_equal(z$rhat_raw,1.001);expect_equal(z$rhat_transformed,1.04)
 expect_match(z$worst_rhat_parameter,'rho');expect_equal(z$ess_transformed,178)
 expect_output(ctmaFitDiagnostics(x),'R-hat exceeds target')
 x$sample$sampler<-'smc';x$sample$rhat<-NA_real_
 z<-ctmaFitDiagnostics(x,print=FALSE)$posterior
 expect_false(z$nonfinite_rhat);expect_false(z$rhat_applicable);expect_true(is.na(z$rhat))
 expect_output(ctmaFitDiagnostics(x),'not applicable for SMC')
})
test_that('unknown covariance codes fail before Julia initialization',{
 for(code in c(-1L,0L,1L,2L))expect_identical(.ctmaJuliaCheckCovarianceCode(list(choleskymats=code)),code)
 for(code in list(NULL,NA_real_,NaN,Inf,3L,'2',c(1,2)))expect_error(
  .ctmaJuliaCheckCovarianceCode(list(choleskymats=code)),'Unsupported ctsem covariance')
})

test_that('quadrature is an explicit optimized random-effects route', {
 # Availability has dedicated tests; this test checks only route semantics.
 local_mocked_bindings(.ctmaRequireCtsemJuliaLaplace=function()list(available=TRUE),
   .package='CoTiMA2')
 x <- .ctmaResolveFitBackend('ctsemjulia', list(intoverpop='quadrature'),
   estimation='ML', workflow='ctmaFit')
 expect_identical(x$backend_args$intoverpop, 'quadrature')
 expect_identical(.ctmaCtsemJuliaRoute(x$backend_args, 'ML', TRUE), 'quadrature')
 expect_error(.ctmaCtsemJuliaRoute(x$backend_args, 'ML', FALSE), 'random effects')
 expect_error(.ctmaCtsemJuliaRoute(x$backend_args, 'Bayes', TRUE), 'optimized Laplace')
 expect_identical(.ctmaCtsemJuliaMode(list(), 'quadrature', FALSE, TRUE), 'quadrature')
})

test_that('checked likelihood, resolution and accuracy remain separate', {
 x <- list(optim=list(converged=TRUE, gradient_norm=0.01,
   predicted_gain=5e-10, convergence_floor=1e-6, resolution=0),
   estimate=list(loglik=-2502.268, loglik_laplace=-2502.321,
     loglik_se=0.0198, loglik_method='checked'),
   laplace=list(correction=list(status='checked', applied=FALSE,
     gap=0.0532, check=list(status='unresolved', max_shift=0.1006))))
 z <- .ctmaAnnotateCtsemJuliaResult(x, use_laplace=TRUE)
 d <- ctmaFitDiagnostics(z, print=FALSE)
 expect_true(d$optimizer_converged)
 expect_identical(d$laplace_accuracy$status, 'unresolved')
 expect_equal(d$likelihood$loglik_laplace, -2502.321)
 expect_identical(d$likelihood$method, 'checked')
 expect_equal(d$optimizer_resolution$convergence_floor, 1e-6)
 expect_equal(z$ctsemjulia_diagnostics$likelihood, d$likelihood)
 expect_output(ctmaFitDiagnostics(z), 'keep separate from the reported likelihood')
 expect_output(ctmaFitDiagnostics(z), 'accuracy: unresolved')
 old <- ctmaFitDiagnostics(list(optim=list(converged=TRUE)), print=FALSE)
 expect_identical(old$likelihood$method, 'not recorded')
 expect_true(is.na(old$likelihood$loglik))
})
