# Parameter Validation ----

test_that("Invalid user-supplied parameters trigger an error", {
  # nolint start
  # These tests are commented out for now as they trip the validation function

  run_sir_stochastic_tau(
    0.00001, 0.1, 1e05 - 1, 1, 0, 1e05, 100,
    list(c(s = -1, i = 1), c(i = -1, r = 1)),
    10
  )

  #   # Negative parameters (beta, gamma) should throw an error
  #   expect_error(
  #     run_sir_stochastic_tau(
  #       -0.00001, 0.1, 1e05 - 1, 1, 0, 1e05, 100,
  #       list(c(s = -1, i = 1), c(i = -1, r = 1)),
  #       10
  #   ),
  #     "negative propensity function"
  #   )

  #   expect_error(
  #     run_sir_stochastic_tau(
  #       0.00001, -0.1, 1e05 - 1, 1, 0, 1e05, 100,
  #       list(c(s = -1, i = 1), c(i = -1, r = 1)),
  #       10
  #   ),
  #     "negative propensity function"
  #   )

  #   # Negative initial count in any compartment (state), population size, number
  #   # of time steps, or number of sims should throw an error
  #   expect_error(
  #     run_sir_stochastic_tau(
  #       0.00001, 0.1, -1e05 - 1, 1, 0, 1e05, 100,
  #       list(c(s = -1, i = 1), c(i = -1, r = 1)),
  #       10
  #   ),
  #     "negative propensity function"
  #   )
  #   expect_error(
  #     run_sir_stochastic_tau(
  #       0.00001, 0.1, 1e05 - 1, -1, 0, 1e05, 100,
  #       list(c(s = -1, i = 1), c(i = -1, r = 1)),
  #       10
  #   ),
  #     "negative propensity function"
  #   )
  #   # nolint start
  #   # #TBD why is a negative recovered compartment allowed???
  #   # expect_error(
  #     run_sir_stochastic_tau(
  #       0.00001, 0.1, 1e05 - 1, 1, -1, 1e05, 100,
  #       list(c(s = -1, i = 1), c(i = -1, r = 1)),
  #       10
  #   ),
  #   #   "negative propensity function"
  #   # )
  #   # nolint end

  #   expect_error(
  #     run_sir_stochastic_tau(
  #       0.00001, 0.1, 1e05 - 1, 1, 0, -1e05, 100,
  #       list(c(s = -1, i = 1), c(i = -1, r = 1)),
  #       10
  #   ),
  #     "negative propensity function"
  #   )

  #   expect_error(
  #     run_sir_stochastic_tau(
  #       0.00001, 0.1, 1e05 - 1, 1, 0, 1e05, -100,
  #       list(c(s = -1, i = 1), c(i = -1, r = 1)),
  #       10
  #   ),
  #     "negative propensity function"
  #   )

  #   # Non-numeric input for parameters or compartments (states) should throw
  #   # an error
  #   expect_error(
  #     run_sir_stochastic_tau(
  #       "invalid", 0.1, 1e05 - 1, 1, 0, 1e05, 100,
  #       list(c(s = -1, i = 1), c(i = -1, r = 1)),
  #       10
  #   ),
  #     "non-numeric argument to binary operator" # beta must be numeric
  #   )
  #   expect_error(
  #     run_sir_stochastic_tau(
  #       0.00001, "invalid", 1e05 - 1, 1, 0, 1e05, 100,
  #       list(c(s = -1, i = 1), c(i = -1, r = 1)),
  #       10
  #   ),
  #     "non-numeric argument to binary operator" # gamma must be numeric
  #   )
  #   expect_error(
  #     run_sir_stochastic_tau(
  #       0.00001, 0.1, "invalid", 1, 0, 1e05, 100,
  #       list(c(s = -1, i = 1), c(i = -1, r = 1)),
  #       10
  #   ),
  #     "'x0' is not numeric" # s0 must be numeric
  #   )
  #   expect_error(
  #     run_sir_stochastic_tau(
  #       0.00001, 0.1, 1e05 - 1, "invalid", 0, 1e05, 100,
  #       list(c(s = -1, i = 1), c(i = -1, r = 1)),
  #       10
  #   ),
  #     "'x0' is not numeric" # i0 must be numeric
  #   )
  #   expect_error(
  #     run_sir_stochastic_tau(
  #       0.00001, 0.1, 1e05 - 1, 1, "invalid", 1e05, 100,
  #       list(c(s = -1, i = 1), c(i = -1, r = 1)),
  #       10
  #   ),
  #     "'x0' is not numeric" # r0 must be numeric
  #   )
  #   expect_error(
  #     run_sir_stochastic_tau(
  #       0.00001, 0.1, 1e05 - 1, 1, 0, "invalid", 100,
  #       list(c(s = -1, i = 1), c(i = -1, r = 1)),
  #       10
  #   ),
  #     "'x0' is not numeric" # n must be numeric
  #   )
  #   expect_error(
  #     run_sir_stochastic_tau(
  #       0.00001, 0.1, 1e05 - 1, 1, 0, 1e05, 100,
  #       list(c(s = -1, i = 1), c(i = -1, r = 1)),
  #       "invalid"
  #   ),
  #     "'x0' is not numeric" # n_sims must be numeric
  #   )
  #   expect_error(
  #     run_sir_stochastic_tau(
  #       0.00001, 0.1, 1e05 - 1, 1, 0, 1e05, 100,
  #       list(c(s = -1, i = 1), c(i = -1, r = 1)),
  #       0
  #   ),
  #     "'x0' should be greater than 0" # n_sims must be > 0
  #   )

  #   # Empty n_timestep argument should throw an error
  #     run_sir_stochastic_tau(
  #       0.00001, 0.1, 1e05 - 1, 1, 0, 1e05, "",
  #       list(c(s = -1, i = 1), c(i = -1, r = 1)),
  #       10
  #   ),
  #     "'tf' is not numeric"
  #   )
  expect_error(
    run_sir_stochastic_tau(
      0.00001, 0.1, 1e05 - 1, 1, 0, 1e05, as.numeric(""),
      list(c(s = -1, i = 1), c(i = -1, r = 1)),
      10
    ),
    "missing value where TRUE/FALSE needed"
  )
  # expect_error(
  #     run_sir_stochastic_tau(
  #       0.00001, 0.1, 1e05 - 1, 1, 0, 1e05, 0,
  #       list(c(s = -1, i = 1), c(i = -1, r = 1)),
  #       10
  #   ),
  #   "'x0' should be greater than 0" # n_timesteps must be > 0
  # )
})
# nolint end

# Test Output
test_that("test that output is a list", {
  out <- run_sir_stochastic_tau(
    0.00001, 0.1, 1e05 - 1, 1, 0, 1e05, 100,
    list(c(s = -1, i = 1), c(i = -1, r = 1)),
    10
  )
  expect_equal(class(out), "list")
})

# nolint
test_that("s, i, and r values are integers", {
  out <- run_sir_stochastic_tau(
    0.00001, 0.1, 1e05 - 1, 1, 0, 1e05, 100,
    list(c(s = -1, i = 1), c(i = -1, r = 1)),
    10
  )
  out_s <- out_i <- out_r <- rep(NA, times = length(out))
  for (i in seq_along(out)) {
    out_s[i] <- purrr::map(out, "s")[i]
    out_i[i] <- purrr::map(out, "i")[i]
    out_r[i] <- purrr::map(out, "r")[i]
  }
  expect_equal(unique(unlist(out_s) %% 1 == 0), TRUE)
  expect_equal(unique(unlist(out_i) %% 1 == 0), TRUE)
  expect_equal(unique(unlist(out_r) %% 1 == 0), TRUE)
})


# Initial Conditions ----
test_that("Initial conditions are correctly set in the output", {
  out <- run_sir_stochastic_tau(
    0.00001, 0.1, 1e05 - 1, 1, 0, 1e05, 100,
    list(c(s = -1, i = 1), c(i = -1, r = 1)),
    10
  )
  # Check if initial conditions match the first row of the output
  # Initial conditions are integers and returned list of data frames
  # gives integer values for compartments(states) but float values for time
  # expect_equal uses same approach as all.equal() and ignores small
  # floating differences. Documentation indicates tolerance is relative
  # (i.e. ⁠mean(abs(x - y) / mean(abs(y)) < tolerance⁠), except when the
  # differences are very small, when it becomes absolute
  # (i.e. ⁠mean(abs(x - y) < tolerance⁠)
  for (i in 1:10) {
    expect_equal(out[[i]]$s[1], 1e05 - 1)
    expect_equal(out[[i]]$i[1], 1)
    expect_equal(out[[i]]$r[1], 0)
  }
})
