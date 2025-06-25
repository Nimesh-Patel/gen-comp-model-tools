test_that("delaytime_gamma delays time", {
  numtest <- 10^4
  tvec <- rep(0, numtest)
  meandelay <- 10
  numchains <- 4
  delaytvec <- delaytime_gamma(tvec, meandelay, numchains)
  delaytsame <- tvec

  for (currchain in seq(numchains)) {
    delaytsame <- delaytime_exp(delaytsame, meandelay / numchains)
  }

  expect_equal(
    mean(delaytvec),
    meandelay,
    tolerance = .1
  )

  expect_equal(
    mean(delaytvec),
    mean(delaytsame),
    tolerance = .1
  )

  expect_equal(
    sqrt(var(delaytvec)),
    sqrt(var(delaytsame)),
    tolerance = .1
  )
})
