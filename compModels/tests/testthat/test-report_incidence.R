test_that("report_incidence aggregates correctly", {
  t_inc <- c(.1, .4, 1.1, 3.1, 3.3)
  cum_inc <- c(.1, .3, .5, .8, 1.2)
  t_report <- seq(7)

  check1 <- report_incidence(t_inc,
    t_report,
    cum_incidence = cum_inc
  )

  check2 <- report_incidence(
    t_inc,
    t_report
  )

  expect_equal(check1, c(.3, .2, 0, .7, 0, 0, 0))

  expect_equal(check2, c(2, 1, 0, 2, 0, 0, 0))
})
