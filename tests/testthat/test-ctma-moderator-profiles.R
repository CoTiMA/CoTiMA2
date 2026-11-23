test_that("moderator profiles use the stored coding and study-level reference", {
  coding <- .ctmaModeratorCoding(matrix(c(1, 1, 3, 5), ncol=1L), list("TI2"),
    c(1, 1, 2, 3), "Age", list(2.5), list(2))
  fit <- list(moderatorCoding=coding)
  automatic <- .ctmaModeratorProfiles(fit, c("TI1","TI2"), "Age", NULL, NULL, c(TI1=1))
  expect_equal(unname(automatic$values[,"TI2"]), (c(1, 3, 5)-2.5)/2)
  expect_equal(unname(automatic$values[,"TI1"]), rep(1,3))
  specified <- .ctmaModeratorProfiles(fit, c("TI1","TI2"), "Age", c(-1,0,1), NULL, NULL)
  expect_equal(unname(specified$values[,"TI2"]), (c(-1,0,1)-2.5)/2)
})

test_that("centered categorical coding preserves both levels", {
  coding <- .ctmaModeratorCoding(matrix(c(1, 1, 2, 2), ncol=1L), list("TI2"),
    c(1, 1, 2, 2), "Group", list(.5), list(1), list(c(1,2)))
  profiles <- .ctmaModeratorProfiles(list(moderatorCoding=coding),c("TI1","TI2"),
    "Group",NULL,NULL,NULL)$values
  expect_equal(unname(profiles[,"TI2"]), c(-.5,.5))
  expect_identical(rownames(profiles), c("1","2"))
})
