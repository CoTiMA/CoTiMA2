test_that("ctmaPub matches common author-name formats by surname", {
  expect_identical(
    CoTiMA2:::.ctmaPubSurnameKey(c(
      "J; de Jonge", "de Jonge, J", "& Bakker, A, B", "Arnold B.; Bakker"
    )),
    c("dejonge", "dejonge", "bakker", "bakker")
  )
})

test_that("ctmaPub computes recency directly from the requested window", {
  frequencies <- c(`2015` = 1, `2016` = 3, `2017` = 7,
    `2018` = 5, `2019` = 9)

  result_mean <- CoTiMA2:::.ctmaPubAuthorCounts(
    frequencies, target_year = 2020, recency = 3,
    years_to_exclude = integer(), aggregate_function = mean)
  expect_equal(result_mean$total, mean(frequencies))
  expect_equal(result_mean$recent, mean(c(7, 5, 9)))

  result_max <- CoTiMA2:::.ctmaPubAuthorCounts(
    frequencies, target_year = 2020, recency = 3,
    years_to_exclude = 2019, aggregate_function = max)
  expect_equal(result_max$total, 7)
  expect_equal(result_max$recent, 7)
})

test_that("ctmaPub extracts authors without treating study labels as authors", {
  source <- c("de Jonge, J", "& Bakker, A, B", "Study 1", "2003")
  result <- CoTiMA2:::.ctmaPubStudyAuthors(source)
  expect_identical(result$authors, c("de Jonge, J", "& Bakker, A, B"))
  expect_identical(result$year, 2003L)
})

test_that("ctmaPub produces deterministic publication scores", {
  publications <- structure(list(
    Authors = c("J; de Jonge", "Arnold B.; Bakker"),
    `Publication Frequencies` = list(
      `J; de Jonge` = c(`1999` = 1, `2000` = 2, `2001` = 3,
        `2002` = 4, `2003` = 100),
      `Arnold B.; Bakker` = c(`1999` = 2, `2000` = 4, `2001` = 6,
        `2002` = 8, `2003` = 100)
    )
  ), class = "CoTiMAFit")
  studies <- structure(list(
    source = list(c("de Jonge, J", "& Bakker, A, B", "Study1", "2003"), NA),
    studyNumbers = list(1L, NA),
    moderators = list(NA, NA),
    n.studies = 1L
  ), class = "CoTiMAFit")

  result <- suppressMessages(ctmaPub(
    publications, studies, recency = 2, indFUN = "sum", colFUN = "mean"))
  expect_equal(result$NEPP, mean(c(10, 20)))
  expect_equal(result$NEPPrecency, mean(c(7, 14)))
  expect_length(result$messages, 0L)

  excluded <- suppressMessages(ctmaPub(
    publications, studies, recency = 2, yearsToExclude = 2002,
    indFUN = "sum", colFUN = "mean"))
  expect_equal(excluded$NEPP, mean(c(6, 12)))
  expect_equal(excluded$NEPPrecency, mean(c(3, 6)))
})

test_that("ctmaPub validates required inputs", {
  expect_error(ctmaPub(), "getPubObj")
  expect_error(
    ctmaPub(list(Authors = "A", `Publication Frequencies` = list(1))),
    "primaryStudyList"
  )
})
