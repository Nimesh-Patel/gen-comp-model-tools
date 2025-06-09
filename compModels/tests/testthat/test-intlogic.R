test_that("intlogic identifies ints", {
  yes_int <- c(1, 2, 3)
  no_int <- c(.99, .01, 0)
  no_int_becausesum <- c(1, 0, 0)

  expect_true(intlogic(yes_int))
  expect_false(intlogic(no_int))
  expect_false(intlogic(no_int_becausesum))
})
