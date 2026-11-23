test_that('explicit partial invariance is not suppressed by family defaults', {
 base<-c('CINT','DIFFUSION','T0VAR')
 f<-CoTiMA2:::.ctma_resolve_invariant_drift_families
 expect_identical(f(NULL,base,TRUE),base)
 expect_identical(f('all',base,TRUE),base)
 expect_setequal(f(c('V2toV1','V1toV2'),base,TRUE),c(base,'DRIFT'))
 expect_setequal(f('none',base,TRUE),c(base,'DRIFT'))
 expect_error(f('V1toV2',base,FALSE),'studyVarying excludes DRIFT')
 expect_identical(f('V1toV2',c(base,'DRIFT'),FALSE),c(base,'DRIFT'))
})
