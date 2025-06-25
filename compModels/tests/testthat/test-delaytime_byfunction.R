test_that("delaytime_byfunction delays time", {
  numtest <- 10^4
  t0 <- 2
  tvec <- rep(t0, numtest)
  meandelay <- 10
  delaytvec_det <- delaytime_byfunction(tvec, meandelay)
  delaytvec_exp <-
    delaytime_byfunction(tvec, function() {
      stats::rexp(1, rate = 1 / meandelay)
    })
  delaytvec_exp_compare <- delaytime_exp(tvec, meandelay)

  expect_equal(
    min(delaytvec_det),
    meandelay + t0
  )

  expect_equal(
    mean(delaytvec_exp),
    mean(delaytvec_exp_compare),
    tolerance = .1
  )
})
