test_that("test that output is a ggplot object", {
  base_states <- c("S", "I", "R")
  sir <- define_states(base_states) |>
    add_infection("I", "S", "I", "beta") |>
    add_transition("I", "R", "tau")
  sircompiled <- compilemodel(sir)

  dyn <- wrap_gillespie(
    c("S" = 999, "I" = 1, "R" = 0),
    sircompiled,
    c(beta = 2, tau = 1),
    25,
    1
  )
  outplot <- plot_stoch_model(dyn, compartments = c("S", "I", "R"))
  expect_equal(ggplot2::is.ggplot(outplot), TRUE)
})
