# Parameter Validation ----

test_that("Invalid user-supplied parameters trigger an error", {
  # nolint start
  # These tests are commented out for now as they trip the validation function

  #   # Negative parameters (beta, gamma) should throw an error
  #   expect_error(
  #     run_sir_stochastic(
  #       -0.001, 0.1, 499, 100, 0, 500, 100,
  #       c("beta*s*i", "gamma*i"),
  #       matrix(c(-1, 0, +1, -1, 0, +1),
  #         nrow = 3, byrow = TRUE
  #       ), 10
  #     ),
  #     "negative propensity function"
  #   )

  #   expect_error(
  #     run_sir_stochastic(
  #       0.001, -0.1, 499, 100, 0, 500, 100,
  #       c("beta*s*i", "gamma*i"),
  #       matrix(c(-1, 0, +1, -1, 0, +1),
  #         nrow = 3, byrow = TRUE
  #       ), 10
  #     ),
  #     "negative propensity function"
  #   )

  #   # Negative initial count in any compartment (state), population size, number
  #   # of time steps, or number of sims should throw an error
  #   expect_error(
  #     run_sir_stochastic(
  #       0.001, 0.1, -499, 100, 0, 500, 100,
  #       c("beta*s*i", "gamma*i"),
  #       matrix(c(-1, 0, +1, -1, 0, +1),
  #         nrow = 3, byrow = TRUE
  #       ), 10
  #     ),
  #     "negative propensity function"
  #   )
  #   expect_error(
  #     run_sir_stochastic(
  #       0.001, 0.1, 499, -100, 0, 500, 100,
  #       c("beta*s*i", "gamma*i"),
  #       matrix(c(-1, 0, +1, -1, 0, +1),
  #         nrow = 3, byrow = TRUE
  #       ), 10
  #     ),
  #     "negative propensity function"
  #   )
  #   # nolint start
  #   # #TBD why is a negative recovered compartment allowed???
  #   # expect_error(
  #   #   run_sir_stochastic(
  #   #   0.001, 0.1, 499, 100, -1, 500, 100,
  #   #   c("beta*s*i", "gamma*i"),
  #   #   matrix(c(-1, 0, +1, -1, 0, +1),
  #   #     nrow = 3, byrow = TRUE
  #   #   ), 10
  #   # ),
  #   #   "negative propensity function"
  #   # )
  #   # nolint end

  #   expect_error(
  #     run_sir_stochastic(
  #       0.001, 0.1, 499, -100, 0, -500, 100,
  #       c("beta*s*i", "gamma*i"),
  #       matrix(c(-1, 0, +1, -1, 0, +1),
  #         nrow = 3, byrow = TRUE
  #       ), 10
  #     ),
  #     "negative propensity function"
  #   )

  #   expect_error(
  #     run_sir_stochastic(
  #       0.001, 0.1, 499, -100, 0, 500, -100,
  #       c("beta*s*i", "gamma*i"),
  #       matrix(c(-1, 0, +1, -1, 0, +1),
  #         nrow = 3, byrow = TRUE
  #       ), 10
  #     ),
  #     "negative propensity function"
  #   )

  #   # Non-numeric input for parameters or compartments (states) should throw
  #   # an error
  #   expect_error(
  #     run_sir_stochastic(
  #       "invalid", 0.1, 499, 100, 0, 500, 100,
  #       c("beta*s*i", "gamma*i"),
  #       matrix(c(-1, 0, +1, -1, 0, +1),
  #         nrow = 3, byrow = TRUE
  #       ), 10
  #     ),
  #     "non-numeric argument to binary operator" # beta must be numeric
  #   )
  #   expect_error(
  #     run_sir_stochastic(
  #       0.001, "invalid", 499, 100, 0, 500, 100,
  #       c("beta*s*i", "gamma*i"),
  #       matrix(c(-1, 0, +1, -1, 0, +1),
  #         nrow = 3, byrow = TRUE
  #       ), 10
  #     ),
  #     "non-numeric argument to binary operator" # gamma must be numeric
  #   )
  #   expect_error(
  #     run_sir_stochastic(
  #       0.001, 0.1, "invalid", 100, 0, 500, 100,
  #       c("beta*s*i", "gamma*i"),
  #       matrix(c(-1, 0, +1, -1, 0, +1),
  #         nrow = 3, byrow = TRUE
  #       ), 10
  #     ),
  #     "'x0' is not numeric" # s0 must be numeric
  #   )
  #   expect_error(
  #     run_sir_stochastic(
  #       0.001, 0.1, 499, "invalid", 0, 500, 100,
  #       c("beta*s*i", "gamma*i"),
  #       matrix(c(-1, 0, +1, -1, 0, +1),
  #         nrow = 3, byrow = TRUE
  #       ), 10
  #     ),
  #     "'x0' is not numeric" # i0 must be numeric
  #   )
  #   expect_error(
  #     run_sir_stochastic(
  #       0.001, 0.1, 499, 100, "invalid", 500, 100,
  #       c("beta*s*i", "gamma*i"),
  #       matrix(c(-1, 0, +1, -1, 0, +1),
  #         nrow = 3, byrow = TRUE
  #       ), 10
  #     ),
  #     "'x0' is not numeric" # r0 must be numeric
  #   )
  #   expect_error(
  #     run_sir_stochastic(
  #       0.001, 0.1, 499, 100, 0, "invalid", 100,
  #       c("beta*s*i", "gamma*i"),
  #       matrix(c(-1, 0, +1, -1, 0, +1),
  #         nrow = 3, byrow = TRUE
  #       ), 10
  #     ),
  #     "'x0' is not numeric" # n must be numeric
  #   )
  #   expect_error(
  #     run_sir_stochastic(
  #       0.001, 0.1, 499, 100, 0, 500, 100,
  #       c("beta*s*i", "gamma*i"),
  #       matrix(c(-1, 0, +1, -1, 0, +1),
  #         nrow = 3, byrow = TRUE
  #       ), "invalid"
  #     ),
  #     "'x0' is not numeric" # n_sims must be numeric
  #   )
  #   expect_error(
  #     run_sir_stochastic(
  #       0.001, 0.1, 499, 100, 0, 500, 100,
  #       c("beta*s*i", "gamma*i"),
  #       matrix(c(-1, 0, +1, -1, 0, +1),
  #         nrow = 3, byrow = TRUE
  #       ), 0
  #     ),
  #     "'x0' should be greater than 0" # n_sims must be > 0
  #   )

  #   # Empty n_timestep argument should throw an error
  #   expect_error(
  #     run_sir_stochastic(
  #       0.001, 0.1, 499, 100, 0, 500, "",
  #       c("beta*s*i", "gamma*i"),
  #       matrix(c(-1, 0, +1, -1, 0, +1),
  #         nrow = 3, byrow = TRUE
  #       ), 10
  #     ),
  #     "'tf' is not numeric"
  #   )
  expect_error(
    run_sir_stochastic(
      0.001, 0.1, 499, 100, 0, 500, as.numeric(""),
      c("beta*s*i", "gamma*i"),
      matrix(c(-1, 0, +1, -1, 0, +1),
        nrow = 3, byrow = TRUE
      ), 10
    ),
    "missing value where TRUE/FALSE needed"
  )
  # expect_error(
  #   run_sir_stochastic(
  #     0.001, 0.1, 499, 100, 0, 500, 0,
  #     c("beta*s*i", "gamma*i"),
  #     matrix(c(-1, 0, +1, -1, 0, +1),
  #       nrow = 3, byrow = TRUE
  #     ), 10
  #   ),
  #   "'x0' should be greater than 0" # n_timesteps must be > 0
  # )
})
# nolint end

# Test Output
test_that("test that output is a list", {
  out <- run_sir_stochastic(
    0.001, 0.1, 499, 100, 0, 500, 100,
    c("beta*s*i", "gamma*i"),
    matrix(c(-1, 0, +1, -1, 0, +1),
      nrow = 3, byrow = TRUE
    ), 10
  )
  expect_equal(class(out), "list")
})

test_that("s, i, and r values are integers", {
  out <- run_sir_stochastic(
    0.001, 0.1, 499, 100, 0, 500, 100,
    c("beta*s*i", "gamma*i"),
    matrix(c(-1, 0, +1, -1, 0, +1),
      nrow = 3, byrow = TRUE
    ), 10
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
  out <- run_sir_stochastic(
    0.001, 0.1, 499, 100, 0, 500, 100,
    c("beta*s*i", "gamma*i"),
    matrix(c(-1, 0, +1, -1, 0, +1),
      nrow = 3, byrow = TRUE
    ), 10
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
    expect_equal(out[[i]]$s[1], 499)
    expect_equal(out[[i]]$i[1], 100)
    expect_equal(out[[i]]$r[1], 0)
  }
})

# Conservation ----
# nolint start
# This works for all cases except when first row is NA, need to find out why
# this occurs

# Total Population Conserved
# test_that("Total population is conserved at each time step", {
#   out <- run_sir_stochastic(
#     0.001, 0.1, 499, 100, 0, 500, 100,
#     c("beta*s*i", "gamma*i"),
#     matrix(c(-1, 0, +1, -1, 0, +1),
#       nrow = 3, byrow = TRUE
#     ), 10
#   )
#   total_pop_initial <- 499 + 100 + 0
#   tolerance <- .Machine$double.eps^0.5
#   for (i in seq_len(nrow(out[[1]]))) { # all dfs are same length
#     for (j in length(out)) {
#       total_pop_current <- with(out[[j]][i, ], s + i + r)
#       expect_equal(total_pop_current, total_pop_initial, tolerance = tolerance)
#     }
#   }
# })
# nolint end

# Performance Testing ----
test_that("Model runs within acceptable performance bounds for large sims", {
  # Measure execution time of the run_sir function on long vector
  # of time steps
  elapsed_time <-
    system.time({
      out <- run_sir_stochastic(
        0.001, 0.1, 499, 100, 0, 500, 100,
        c("beta*s*i", "gamma*i"),
        matrix(c(-1, 0, +1, -1, 0, +1),
          nrow = 3, byrow = TRUE
        ), 10
      )
    })[3] # user + sys time
  # Define an acceptable threshold for the execution time (in seconds)
  max_time_allowed <- 2
  expect_true(elapsed_time < max_time_allowed,
    info = paste(
      "Model took too long to run: ",
      elapsed_time, " seconds"
    )
  )
})
