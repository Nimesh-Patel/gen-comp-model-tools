test_that("generalized_rates handles time-varying rates correctly", {
  # Define rate equations that vary with time
  rate_eqns <- c("a * x + b * y + sin(t)", "d * x - e * y + cos(t)")
  rates_func <- generalized_rates(rate_eqns)
  state <- list(x = 1, y = 2)
  params <- list(a = 0.5, b = -0.3, d = -0.2, e = 0.4)
  times <- c(0, pi / 2, pi)

  expected_results <- list(
    c(
      1 * 0.5 + (-0.3) * 2 + sin(0),
      (-0.2) * 1 - (0.4) * 2 + cos(0)
    ),
    c(
      1 * 0.5 + (-0.3) * 2 + sin(pi / 2),
      (-0.2) * 1 - (0.4) * 2 + cos(pi / 2)
    ),
    c(
      1 * 0.5 + (-0.3) * 2 + sin(pi),
      (-0.2) * 1 - (0.4) * 2 + cos(pi)
    )
  )

  for (i in seq_along(times)) {
    result <- rates_func(state, params, times[i])
    expect_equal(result[[1]], expected_results[[i]][[1]], tolerance = 1e-8)
    expect_equal(result[[2]], expected_results[[i]][[2]], tolerance = 1e-8)
  }
})
