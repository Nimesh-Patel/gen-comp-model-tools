test_that("test that output is a ggplot object", {
  modelout <- run_sir_stochastic(
    0.001, 0.1, 499, 100, 0, 500, 100,
    c("beta*s*i", "gamma*i"),
    matrix(c(-1, 0, +1, -1, 0, +1),
      nrow = 3, byrow = TRUE
    ),
    10
  )
  outplot <- plot_stoch_model(modelout)
  expect_equal(ggplot2::is.ggplot(outplot), TRUE)
})
