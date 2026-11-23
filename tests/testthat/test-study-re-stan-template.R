re_template_fixture <- function(layout = c('legacy', 'additive')) {
  layout <- match.arg(layout)
  t0 <- if (layout == 'legacy')
    'if(intoverpop && nindvarying > 0) T0VAR[intoverpopindvaryingindex, intoverpopindvaryingindex] = rawpopcovbase;' else
    paste0('T0cov[intoverpopindvaryingindex, intoverpopindvaryingindex] = ',
           'T0cov[intoverpopindvaryingindex, intoverpopindvaryingindex] + ',
           'quad_form_diag(rawpopcov, popstatescale);')
  paste(c(
    '  int nindvaryingoffdiagonals; //number of off diagonal parameters needed for popcov matrix',
    '  vector[(nsubsets > 1) ? 1 : 0] subsetpar;',
    '  matrix[nindvarying, nindvarying] rawpopcorr;',
    '  }//end indvarying par setup',
    '  if(si > 0 &&  ntieffects > 0){', t0,
    '  real priormod2 = priormod / nsubsets;'), collapse='\n')
}
re_patch_fixture <- function(text) {
  CoTiMA2:::.ctma_patch_stan(text, 'independent', 'noncentered',
    data.frame(structure='independent',size=4L,start=1L,block_id=1L))
}

test_that('both covariance layouts are recognized and malformed templates fail', {
  for (layout in c('legacy','additive')) {
    text <- re_template_fixture(layout)
    expect_identical(CoTiMA2:::.ctma_stan_t0_layout(text)$name, layout)
    expect_error(re_patch_fixture(paste(text,text)), 'exactly one supported')
  }
  expect_error(re_patch_fixture('unknown template'), 'exactly one supported')
  expect_error(re_patch_fixture(paste(re_template_fixture('legacy'),
    re_template_fixture('additive'))), 'exactly one supported')
})

test_that('additive covariance retains the conditional block and scales once', {
  text <- re_patch_fixture(re_template_fixture('additive'))
  expect_match(text, paste0('if(si > 0) T0cov[intoverpopindvaryingindex, intoverpopindvaryingindex] = ',
    'T0cov[intoverpopindvaryingindex, intoverpopindvaryingindex] + studyre_T0cov[study_id[si]];'),fixed=TRUE)
  expect_match(text, 'else T0cov[', fixed=TRUE)
  expect_match(text, 'studyre_T0cov[s] = quad_form_diag(studyre_T0cov[s], popstatescale);',fixed=TRUE)
  expect_false(grepl('studyre_T0cov[s][q,] *=', text, fixed=TRUE))
  expect_match(text,'studyre_rawpopcovbase[s][r,c] = sqrtpcov[studyre_index[q]] + studyre_effect[s,q];',fixed=TRUE)
  expect_false(grepl('inv_logit(sqrtpcov', text, fixed=TRUE))
})

test_that('legacy covariance construction remains supported', {
  text <- re_patch_fixture(re_template_fixture('legacy'))
  expect_match(text, 'if(si > 0) T0VAR[', fixed=TRUE)
  expect_match(text, 'studyre_T0cov[s][q,] *=', fixed=TRUE)
  expect_match(text, 'inv_logit(sqrtpcov', fixed=TRUE)
  expect_false(grepl('popstatescale',text,fixed=TRUE))
})

test_that('ctsem 4 preparation uses Stan and accepts bivariate manifest variation', {
  skip_if_not_installed('ctsem','4.0.0')
  skip_if_not_installed('rstan')
  set.seed(22)
  dat <- data.frame(id=rep(1:6,each=3),time=rep(0:2,6),
    y1=rnorm(18),y2=rnorm(18),TI1=rep(c(0,0,1,1,0,0),each=3),
    TI2=rep(c(0,0,0,0,1,1),each=3))
  m <- ctsem::ctModel(LAMBDA=diag(2),manifestNames=c('y1','y2'),
    latentNames=c('TP','OCB'),TIpredNames=c('TI1','TI2'),
    MANIFESTMEANS=matrix(c('mu1','mu2'),2),MANIFESTVAR=diag(.1,2))
  m$pars$indvarying <- m$pars$matrix=='MANIFESTMEANS'
  wrapper <- structure(list(data=dat,ctModel=m,argumentList=list()),class='CoTiMAFit')
  setup <- CoTiMA2:::.ctma_generate_ctsem_setup(wrapper,FALSE)
  expect_false(is.null(setup$stanmodeltext))
  prepared <- CoTiMA2::ctmaRePrep(setup,study_dummy_ti=c('TI1','TI2'),
    reDRIFT=matrix(TRUE,2,2),reDIFFUSION=matrix(FALSE,2,2),
    reT0VAR=matrix(FALSE,2,2),covariance='independent',
    parameterization='noncentered',verbose=FALSE)
  expect_equal(nrow(prepared$specification),4)
  expect_equal(prepared$standata$nstudies,3)
  expect_equal(prepared$standata$TIPREDEFFECTsetup[prepared$specification$index,],
    matrix(0L,4,2))
  expect_gt(prepared$standata$nindvarying,0)
  # Syntax/type validation, without C++ compilation or sampling.
  expect_true(rstan::stanc(model_code=prepared$stanmodeltext,
    model_name='cotima_study_re_template_test')$status)
})
