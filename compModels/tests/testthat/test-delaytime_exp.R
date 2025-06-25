test_that("delaytime_exp delays time", {
  numtest <- 10^4
  tvec <- rep(0, numtest)
  meandelay <- 10
  delaytvec <- delaytime_exp(tvec, meandelay)

  # same processes
  expect_equal(
    mean(delaytvec),
    meandelay,
    tolerance = .1
  )

  expect_equal(
    mean(delaytvec),
    sqrt(var(delaytvec)),
    tolerance = .1
  )
})
