test_that("test that output is a ggplot object", {
  base_states <- c("S", "I", "R")
  sir <- define_states(base_states) |>
    add_infection("I", "S", "I", "beta") |>
    add_transition("I", "R", "tau")
  sircompiled <- compilemodel(sir)

  odefun <- model2desolvefunction(sircompiled)
  x0 <- define_initialstate(sircompiled, c("S" = 999, "I" = 1, "R" = 0)) |>
    output_initialstate()
  parameters <- c(beta = 2, tau = 1)
  times <- seq(0, 25, by = 0.01)
  dyn <- deSolve::ode(y = x0, times = times, func = odefun, parms = parameters)
  outplot <- plot_determin_model(dyn)
  expect_equal(ggplot2::is.ggplot(outplot), TRUE)
})
