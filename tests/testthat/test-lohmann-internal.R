test_that("Lohmann data keeps compact lower triangles and padding aligned", {
  a <- diag(4); a[1,2] <- a[2,1] <- .2
  b <- diag(6)
  d <- .ctmaLohmannData(list(a,b), c(100,200), c(1,2),list(diag(6)/100,diag(15)/200))
  expect_silent(.ctmaLohmannValidateData(d))
  expect_equal(d$r_i[1,1:6], a[lower.tri(a)])
  expect_equal(d$r_i[1,7:15], rep(-99,9))
  expect_equal(d$V_i[1,1:6,1:6],diag(6)/100)
  d$V_i[1,1,1] <- -1
  expect_error(.ctmaLohmannValidateData(d),"sampling covariance")
  expect_error(.ctmaLohmannData(list(a),100,1,NULL),"samplingCovariances")
  a[1,1] <- 2
  expect_error(.ctmaLohmannData(list(a),100,1,list(diag(6))),"correlation matrices")
})

test_that("Lohmann compiler is parsed without executing arbitrary code", {
  f <- tempfile(fileext=".R"); on.exit(unlink(f))
  writeLines(c('stop("must never run")', 'OS_URE_char = "data {} parameters {} model {}"'), f)
  x <- .ctmaLohmannSource(f,"URE")
  expect_equal(x$code,"data {} parameters {} model {}")
  expect_equal(x$md5,unname(tools::md5sum(f)))
  expect_error(.ctmaLohmannSource(f,"CRE"),"literal Stan")
  d <- .ctmaLohmannData(list(diag(4)),100,1,list(diag(6)/100))
  expect_warning(p <- ctmaLohmannFit(stanData=d,supplementDirectory=f,fit=FALSE),"Literal supplement")
  expect_null(p$fit)
  expect_identical(p$stanData,d)
  expect_error(ctmaLohmannFit(stanData=d,supplementDirectory=f,fit=FALSE,optimize=TRUE),"Unknown")
  expect_error(ctmaLohmannFit(stanData=d,timeLag=1,supplementDirectory=f,fit=FALSE),"not both")
  expect_error(ctmaLohmannFit(reStructure="bad"),"reStructure")
})

test_that("Lohmann accepts authors' list-of-matrices V_i with one study", {
  f <- tempfile(fileext=".R"); on.exit(unlink(f))
  writeLines('OS_URE_char = "data {} parameters {} model {}"',f)
  d <- .ctmaLohmannData(list(diag(4)),100,1,list(diag(6)/100))
  d$V_i <- list(diag(6)/100)
  expect_warning(x <- ctmaLohmannFit(stanData=d,supplementDirectory=f,fit=FALSE),"Literal")
  expect_identical(dim(x$stanData$V_i),c(1L,6L,6L))
  expect_error(ctmaLohmannFit(stanData=d,supplementDirectory=f,fit=FALSE,
    activeDirectory=tempfile()),"activeDirectory")
})

test_that("original split loading loops equal the wave-major identity map", {
  for (w in c(2L,3L,7L)) {
    L <- matrix(0,2*w,2*w+2)
    for (i in seq_len(w)) {
      L[1+2*(i-1),1] <- 1; L[2*i,2] <- 1
      L[i,2+i] <- 1; L[w+i,w+2+i] <- 1
    }
    expect_equal(L[,3:(2*w+2)],diag(2*w))
    expect_equal(L[,1],rep(c(1,0),w))
  }
})

test_that("univariate extension allows two waves and records its distinct provenance", {
  d <- .ctmaLohmannData(list(diag(2)),100,1,list(matrix(.01)),n.latent=1)
  expect_silent(.ctmaLohmannValidateData(d,1))
  expect_warning(p <- ctmaLohmannFit(stanData=d,n.latent=1,fit=FALSE),"Univariate extension")
  expect_equal(p$parameter_order,c("A11_raw","trait_X_logit"))
  expect_true(is.na(p$source$path))
  expect_match(p$modelCode,"diffusion\\[s\\]=-2\\*drift")
  expect_error(ctmaLohmannFit(stanData=d,n.latent=3,fit=FALSE),"n.latent")
  expect_error(.ctmaLohmannValidateData(d,2),"dimensions")
  for (variant in c("FE","URE","ERE","CRE"))
    expect_match(.ctmaLohmannUnivariateCode(variant),"multi_normal")
})
