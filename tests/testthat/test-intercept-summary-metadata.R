test_that('intercept metadata is captured before covariance-summary mutation', {
  assignments <- list(); fields <- list()
  walk <- function(x) {
    if(!is.call(x)) return(invisible(NULL))
    if(identical(x[[1]],as.name('<-')) && identical(x[[2]],as.name('randomInterceptStructure')))
      assignments[[length(assignments)+1L]] <<- x[[3]]
    args <- as.list(x)[-1L]
    if(!is.null(names(args))) for(i in which(names(args)=='randomInterceptStructure'))
      fields[[length(fields)+1L]] <<- args[[i]]
    for(i in seq_along(args)) {
      if(identical(unname(args[i]),unname(alist(x=)))) next
      if(is.call(args[[i]])) walk(args[[i]])
    }
  }
  walk(body(ctmaFit))
  expect_length(assignments,1L)
  expect_length(fields,2L)
  expect_true(all(vapply(fields,identical,logical(1),as.name('randomInterceptStructure'))))
  for(mode in list(FALSE,'CINT','MANIFEST')) for(pooled in c(FALSE,TRUE)) {
    env <- list2env(list(randomIntercepts=mode,pooledRandomIntercepts=pooled))
    env$randomInterceptStructure <- eval(assignments[[1]],env)
    expected <- if(pooled)'pooled' else if(identical(mode,FALSE))'none' else 'study-specific'
    # Simulate the downstream covariance summary that triggered the user's error.
    env$randomIntercepts <- list(popsd=1,popcov_mean=diag(2),popcor_mean=diag(2))
    for(field in fields) expect_identical(eval(field,env),expected)
  }
})
