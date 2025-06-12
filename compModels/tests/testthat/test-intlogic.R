test_that("intlogic identifies ints", {
  yes_int <- c(1, 2, 3)
  no_int <- c(.99, .01, 0)
  no_int_becausesum <- c(1, 0, 0)

  tbl_test1 <- tibble::tibble(popsize = yes_int)
  expect_true(intlogic(tbl_test1))

  tbl_test2 <- tibble::tibble(popsize = no_int)
  expect_false(intlogic(tbl_test2))

  tbl_test3 <- tibble::tibble(popsize = no_int_becausesum)
  expect_false(intlogic(tbl_test3))

  # with metapopulations
  mtbl_test1 <- dplyr::bind_rows(
    tbl_test1 |> dplyr::mutate(metapopulation = 1),
    tbl_test1 |> dplyr::mutate(metapopulation = 2)
  )
  expect_true(intlogic(mtbl_test1))

  mtbl_test2 <- dplyr::bind_rows(
    tbl_test2 |> dplyr::mutate(metapopulation = 1),
    tbl_test2 |> dplyr::mutate(metapopulation = 2)
  )
  expect_false(intlogic(mtbl_test2))

  mtbl_test3 <- dplyr::bind_rows(
    tbl_test3 |> dplyr::mutate(metapopulation = 1),
    tbl_test3 |> dplyr::mutate(metapopulation = 2)
  )
  expect_true(intlogic(mtbl_test3))
})
