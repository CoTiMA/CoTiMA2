test_that("selection retains legacy prep format and trims only confirmed empty closure", {
  p<-list(studyNumbers=list(4,1,NA),deltas=list(2,1,NA),
    sampleSizes=list(200,100,NA),empcovs=list(diag(4),diag(4),matrix(NA,0,0)),
    rawData=list(list(),list(),list(fileName=NULL,studyNumbers=NULL,missingValues=-99)),
    moderators=list(NA,NA,NA),n.studies=2,
    moderatorValues=list(c('a','b'),c('c','d')))
  original<-p
  expect_equal(.ctmaPrepSelectionIDs(p),c(4,1))
  expect_identical(p,original)
  clean<-p;clean$studyNumbers<-list(4,1)
  expect_equal(.ctmaPrepSelectionIDs(clean),c(4,1))
  bad<-p;bad$sampleSizes[[3]]<-99
  expect_error(.ctmaPrepSelectionIDs(bad),'contains data')
  bad<-p;bad$rawData[[3]]$fileName<-'real-data.csv'
  expect_error(.ctmaPrepSelectionIDs(bad),'raw-data')
  bad<-p;bad$studyNumbers<-list(4,NA,1)
  expect_error(.ctmaPrepSelectionIDs(bad),'invalid missing')
  bad<-p;bad$n.studies<-3
  expect_error(.ctmaPrepSelectionIDs(bad),'invalid missing')
  bad<-p;bad$empcovs<-list(diag(4))
  expect_error(.ctmaPrepSelectionIDs(bad),'misaligned')
})
