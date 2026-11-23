test_that('compact identification printing preserves full messages and unknown notes', {
  wave <- paste0('At most two observed waves per person in a separately fitted primary study ',
    'with its own random-intercept covariance: a design limitation. The model is retained unchanged.')
  raw <- c(paste0('Study 20: ',wave),paste0('Study 3: ',wave),
    'Study 20: 3 weak parameter direction(s); optimizer convergence does not establish identification.',
    'Study 20: unavailable uncertainty for a, b',
    'Study 3: random-effect correlation at its boundary.',
    'Study 3: a future warning must remain visible.', 'A note without a study prefix.')
  x <- structure(raw,class=c('ctmaIdentificationNotes','character'))
  out <- capture.output(print(x))
  expect_identical(as.character(x),raw)
  value <- NULL; invisible(capture.output(value <- print(x)))
  expect_identical(value,x)
  expect_equal(sum(grepl('At most two observed waves',out)),1)
  expect_match(paste(out,collapse='\n'),'3 weak|Weak: number')
  expect_true(any(grepl('a future warning',out)))
  expect_true(any(grepl('without a study prefix',out)))
  t <- .ctmaInitNoteFormat(x)$table
  expect_equal(t$Study,c('3','20'));expect_true(is.na(t$Weak[1]))
  expect_equal(t$Unavailable[2],2L);expect_equal(t$Flags,c('T,B,O','T'))
  detail <- capture.output(print(x,details=TRUE))
  expect_equal(sum(grepl('At most two observed waves',detail)),2)
  expect_error(print(x,details=NA),'TRUE or FALSE')
  # Base list printing must dispatch the compact method, as fit$summary does.
  expect_equal(sum(grepl('At most two observed waves',capture.output(print(list(identification_notes=x))))),1)
})
test_that('two-wave warning references are per call, never a global suppression', {
  text <- 'At most two observed waves per person: full design explanation.'
  expect_identical(.ctmaInitWaveMessage(text),paste0('[T] ',text))
  expect_match(.ctmaInitWaveMessage(text,FALSE),'earlier in this ctmaInit')
  expect_false(grepl('full design explanation',.ctmaInitWaveMessage(text,FALSE)))
  expect_identical(.ctmaInitWaveMessage(text),paste0('[T] ',text))
})
