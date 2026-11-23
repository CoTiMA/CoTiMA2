test_that("observed occasions count unique person-times with actual outcomes", {
 d<-data.frame(id=c(1,1,1,2,2,2,2),time=c(0,1,1,0,1,2,3),V1=c(1,1,1,1,1,1,NA))
 r<-.ctmaObservedStudyDesign(d,"V1")
 expect_identical(r$nWaves,3L)
 expect_identical(r$minWaves,2L)
 expect_identical(r$nSubjectsTwoWaves,1L)
 expect_identical(r$nSubjectsThreeOrMoreWaves,1L)
 expect_identical(.ctmaObservedStudyDesign(d[FALSE,],"V1")$nSubjects,0L)
})

test_that("partitions preserve ctmaPrep structure, custom study fields and global metadata", {
 s<-structure(list(n.studies=3L,studyNumbers=list(9,31,8),empcovs=list(diag(2),diag(2),diag(2)),
  custom=list("a","b","c"),sampleSizes=list(10,20,30),
  moderatorLabels=c("global1","global2","global3"),
  moderatorValues=list(1,2,3),rawData=list("path1","path2","path3"),
  summary=data.frame(study=c(9,31,8)),reUseEmprawData=TRUE),class="CoTiMAFit")
 out<-.ctmaPrepStudySubset(s,c(3L,1L))
 expect_s3_class(out,"CoTiMAFit")
 expect_identical(out$studyNumbers,list(8,9))
 expect_identical(out$custom,list("c","a"))
 expect_identical(out$rawData,list("path3","path1"))
 expect_identical(out$moderatorLabels,s$moderatorLabels)
 expect_identical(out$moderatorValues,s$moderatorValues)
 expect_identical(out$summary$study,c(8,9))
 expect_identical(s$n.studies,3L)
 expect_null(.ctmaPrepStudySubset(s,integer()))
})

test_that("mixed and incomplete designs remain visible without excluding studies", {
 d<-function(times)data.frame(id=rep(1,length(times)),time=times,V1=rep(1,length(times)))
 s<-list(n.studies=4L,studyNumbers=as.list(c(10,20,30,40)),empcovs=rep(list(diag(2)),4))
 m<-list(manifestNames="V1",pars=data.frame(matrix="MANIFESTMEANS",indvarying=TRUE))
 x<-list(n.studies=4L)
 out<-.ctmaAttachStudyDesign(x,s,list(d(0:1),d(0:2),d(0),d(numeric())),m)
 expect_identical(out$n.studies,4L)
 expect_identical(out$studyLists$twoWaves$studyNumbers,list(10))
 expect_identical(out$studyLists$threeOrMoreWaves$studyNumbers,list(20))
 expect_identical(out$studyDesign$waveGroup,c("twoWaves","threeOrMoreWaves","fewerThanTwoWaves","unknown"))
 expect_identical(out$studyDesign$randomInterceptIdentificationRisk,c(TRUE,FALSE,TRUE,NA))
})

test_that("pooled covariance is not given a separately-fitted two-wave warning", {
 x<-list(studyDesign=data.frame(study=c("6","52"),maxWaves=c(2L,3L)))
 expect_warning(.ctmaWarnStudySpecificTwoWaves(x,"study-specific"),"studies 6")
 expect_no_warning(.ctmaWarnStudySpecificTwoWaves(x,"pooled"))
 expect_no_warning(.ctmaWarnStudySpecificTwoWaves(x,"none"))
})

test_that("legacy wide-data reconstruction counts time 10 rather than lexical positions", {
 f<-list(empraw=data.frame(V1_T0=c(1,NA),V1_T2=c(1,NA),V1_T10=c(1,NA)))
 x<-list(n.studies=1L,studyFitList=list(f))
 d<-.ctmaInitObservedLong(x)[[1]]
 expect_equal(sort(d$time),c(0,2,10))
 expect_identical(.ctmaObservedStudyDesign(d)$maxWaves,3L)
})

test_that("ctmaInit public setup returns reusable prepared subsets", {
 e<-new.env();utils::data("CoTiMAstudyList_3",package="CoTiMA2",envir=e)
 s<-e$CoTiMAstudyList_3
 out<-suppressMessages(suppressWarnings(ctmaInit(primaryStudies=s,n.latent=2,
  fit=FALSE,silentOverwrite=TRUE,activeDirectory=paste0(tempdir(),"/"))))
 expect_equal(nrow(out$studyDesign),s$n.studies)
 expect_equal(length(out$data),s$n.studies)
 for(group in c("twoWaves","threeOrMoreWaves")) {
  subset<-out$studyLists[[group]]
  if(is.null(subset)) next
  setup<-suppressMessages(suppressWarnings(ctmaInit(primaryStudies=subset,n.latent=2,
   fit=FALSE,silentOverwrite=TRUE,activeDirectory=paste0(tempdir(),"/"))))
  expect_equal(nrow(setup$studyDesign),subset$n.studies)
  expect_identical(setup$studyDesign$study,as.character(unlist(subset$studyNumbers)))
 }
})
test_that("public pooled setup distinguishes individual and study-specific RI warnings", {
 e<-new.env();utils::data("CoTiMAstudyList_3",package="CoTiMA2",envir=e)
 utils::data("CoTiMAInitFit_3",package="CoTiMA2",envir=e)
 fresh<-suppressMessages(suppressWarnings(ctmaInit(primaryStudies=e$CoTiMAstudyList_3,
  n.latent=2,fit=FALSE,silentOverwrite=TRUE,activeDirectory=paste0(tempdir(),"/"))))
 init<-e$CoTiMAInitFit_3;init$emprawList<-fresh$data
 # Make the first study genuinely two-wave in its observed wide records.
 columns<-grep("_T[0-9]+$",names(init$emprawList[[1]]),value=TRUE)
 later<-columns[as.integer(sub("^.*_T","",columns))>1L]
 init$emprawList[[1]][,later]<-NA_real_
 init$studyFitList<-lapply(init$emprawList,function(d)list(empraw=d))
 init$studyDesign<-.ctmaStudyDesignTable(e$CoTiMAstudyList_3,.ctmaInitObservedLong(init))
 collect<-function(...) {
  messages<-character()
  setup<-withCallingHandlers(suppressMessages(ctmaFit(ctmaInitFit=init,fit=FALSE,
   activeDirectory=paste0(tempdir(),"/"),...)),warning=function(w){
    messages<<-c(messages,conditionMessage(w));invokeRestart("muffleWarning")})
  list(setup=setup,warnings=messages)
 }
 shared<-collect(indVarying="MANIFEST")
 expect_false(any(grepl("Study-specific random-intercept covariances requested",shared$warnings,fixed=TRUE)))
 separate<-collect(randomIntercepts="MANIFEST")
 expect_true(any(grepl("Study-specific random-intercept covariances requested",separate$warnings,fixed=TRUE)))
 expect_identical(separate$setup$studyList,shared$setup$studyList)
})
