test_that("maintainpopsizetype maintains int type", {
  tblcheck_int <- tibble::tibble(testcol = c(1, 2, 3))
  expect_equal(is.integer(tblcheck_int$testcol), FALSE)
  tblcheck_int <- maintainpopsizetype(tblcheck_int, "testcol")
  expect_equal(is.integer(tblcheck_int$testcol), TRUE)

  tblcheck_notint <- tibble::tibble(testcol = c(1.1, 2.1, 3.1))
  tblcheck_notint <- maintainpopsizetype(tblcheck_notint, "testcol")
  expect_equal(is.integer(tblcheck_notint$testcol), FALSE)
})
