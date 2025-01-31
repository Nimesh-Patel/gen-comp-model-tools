# Basic Functionality ----
test_that("Test output is a data frame with expected dimensions", {
  out <-
    run_sir(
      init = c(s = 1e05 - 1, i = 1, r = 0),
      time = seq(0.1, 100, by = 0.1),
      parms = c(beta = 0.5, gamma = 0.1)
    )
  expect_equal(class(out), "data.frame")
  # Expect ncol = # compartments + 1 & nrows = length of provided time sequence
  expect_equal(ncol(out), 4)
  expect_equal(nrow(out), length(seq(0.1, 100, by = 0.1)))
})

# Parameter Validation ----
# nolint start
test_that("Invalid user-supplied parameters trigger an error", {
  # Negative parameters (beta, gamma) should throw an error
  expect_error(
    run_sir(
      init = c(s = 1e05 - 1, i = 1, r = 0),
      time = seq(0.1, 100, by = 0.1),
      parms = c(beta = -0.00001, gamma = 0.1)
    ),
    "beta should be non-negative"
  )
  expect_error(
    run_sir(
      init = c(s = 1e05 - 1, i = 1, r = 0),
      time = seq(0.1, 100, by = 0.1),
      parms = c(beta = 0.5, gamma = -0.1)
    ),
    "gamma should be non-negative"
  )

  # Negative initial count in any compartment (state) should throw an error
  expect_error(
    run_sir(
      init = c(s = -1, i = 1, r = 0),
      time = seq(0.1, 100, by = 0.1),
      parms = c(beta = 0.5, gamma = 0.1)
    ),
    "s should be non-negative"
  )
  expect_error(
    run_sir(
      init = c(s = 1e05 - 1, i = -1, r = 0),
      time = seq(0.1, 100, by = 0.1),
      parms = c(beta = 0.5, gamma = 0.1)
    ),
    "i should be non-negative"
  )
  expect_error(
    run_sir(
      init = c(s = 1e05 - 1, i = 1, r = -1),
      time = seq(0.1, 100, by = 0.1),
      parms = c(beta = 0.5, gamma = 0.1)
    ),
    "r should be non-negative"
  )

  # Non-numeric input for parameters or compartments (states) should throw
  # an error
  expect_error(
    run_sir(
      init = c(s = 1e05 - 1, i = 1, r = 0),
      time = seq(0.1, 100, by = 0.1),
      parms = c(beta = "invalid", gamma = 0.1)
    ),
    "beta must be numeric"
  )
  expect_error(
    run_sir(
      init = c(s = 1e05 - 1, i = 1, r = 0),
      time = seq(0.1, 100, by = 0.1),
      parms = c(beta = 0.5, gamma = "invalid")
    ),
    "gamma must be numeric"
  )
  expect_error(
    run_sir(
      init = c(s = "invalid", i = 1, r = 0),
      time = seq(0.1, 100, by = 0.1),
      parms = c(beta = 0.5, gamma = 0.1)
    ),
    "s must be numeric"
  )
  expect_error(
    run_sir(
      init = c(s = 1e05 - 1, i = "invalid", r = 0),
      time = seq(0.1, 100, by = 0.1),
      parms = c(beta = 0.5, gamma = 0.1)
    ),
    "i must be numeric"
  )
  expect_error(
    run_sir(
      init = c(s = 1e05 - 1, i = 1, r = "invalid"),
      time = seq(0.1, 100, by = 0.1),
      parms = c(beta = 0.5, gamma = 0.1)
    ),
    "r must be numeric"
  )

  # Empty time vector should throw an error
  expect_error(
    run_sir(
      init = c(s = 1e05 - 1, i = 1, r = 0),
      time = c(),
      parms = c(beta = 0.5, gamma = 0.1)
    ),
    "time vector must have at least one element"
  )

  # Non-increasing time vector should throw an error
  expect_error(
    run_sir(
      init = c(s = 1e05 - 1, i = 1, r = 0),
      time = c(10:5),
      parms = c(beta = 0.5, gamma = 0.1)
    ),
    "time vector must be strictly increasing"
  )
})

# nolint end

# Initial Conditions ----
test_that("Initial conditions are correctly set in the output", {
  out <- run_sir(
    init = c(s = 1e05 - 1, i = 1, r = 0),
    time = seq(0.1, 100, by = 0.1),
    parms = c(beta = 0.5, gamma = 0.1)
  )
  # Check if initial conditions match the first row of the output
  # Initial conditions are integers and returned data frame gives float
  # expect_equal uses same approach as all.equal() and ignores small
  # floating differences. Documentation indicates tolerance is relative
  # (i.e. ⁠mean(abs(x - y) / mean(abs(y)) < tolerance⁠), except when the
  # differences are very small, when it becomes absolute
  # (i.e. ⁠mean(abs(x - y) < tolerance⁠)
  expect_equal(out$s[1], 1e05 - 1)
  expect_equal(out$i[1], 1)
  expect_equal(out$r[1], 0)
})

# Conservation ----
# Total Population Conserved
test_that("Total population is conserved at each time step", {
  out <- run_sir(
    init = c(s = 1e05 - 1, i = 1, r = 0),
    time = seq(0.1, 100, by = 0.1),
    parms = c(beta = 0.5, gamma = 0.1)
  )
  total_pop_initial <- (1e05 - 1) + 1 + 0
  tolerance <- .Machine$double.eps^0.5
  for (i in seq_len(nrow(out))) {
    total_pop_current <- with(out[i, ], s + i + r)
    expect_equal(total_pop_current, total_pop_initial, tolerance = tolerance)
  }
})

# Deterministic Output ----
test_that("Model output is deterministic", {
  # Run SIR model simulation twice with the same parameters
  out1 <- run_sir(
    init = c(s = 1e05 - 1, i = 1, r = 0),
    time = seq(0.1, 100, by = 0.1),
    parms = c(beta = 0.5, gamma = 0.1)
  )
  out2 <- run_sir(
    init = c(s = 1e05 - 1, i = 1, r = 0),
    time = seq(0.1, 100, by = 0.1),
    parms = c(beta = 0.5, gamma = 0.1)
  )
  # Check if both results are identical using all.equal()
  expect_true(all.equal(out1, out2))
})

# Time Step Size Sensitivity ----
test_that("Model output is stable across different time step granularities", {
  # Run SIR model simulation twice with different granularity time steps
  out1 <- run_sir(
    init = c(s = 1e05 - 1, i = 1, r = 0),
    time = seq(0.1, 100, by = 1),
    parms = c(beta = 0.5, gamma = 0.1)
  )
  out2 <- run_sir(
    init = c(s = 1e05 - 1, i = 1, r = 0),
    time = seq(0.1, 100, by = 0.1),
    parms = c(beta = 0.5, gamma = 0.1)
  )
  # Interpolate fine results at coarse time points for comparison
  interp_fine_s <- approx(out2$time, out2$s, xout = seq(0.1, 100, by = 1))$y
  interp_fine_i <- approx(out2$time, out2$i, xout = seq(0.1, 100, by = 1))$y
  interp_fine_r <- approx(out2$time, out2$r, xout = seq(0.1, 100, by = 1))$y
  # Define a tolerance level for comparing values
  # (this may need to be adjusted based on expected differences)
  tolerance <- max(abs(out2$s - interp_fine_s)) *
    .Machine$double.eps^0.5

  # Compare coarse and interpolated fine results using all.equal()
  # within a tolerance level
  expect_true(all.equal(out1$s, interp_fine_s, tolerance = tolerance))
  expect_true(all.equal(out1$i, interp_fine_i, tolerance = tolerance))
  expect_true(all.equal(out1$r, interp_fine_r, tolerance = tolerance))
})

# Performance Testing ----
test_that("Model runs within acceptable performance bounds for large sims", {
  # Measure execution time of the run_sir function on long vector
  # of time steps
  elapsed_time <-
    system.time({
      out1 <- run_sir(
        init = c(s = 1e05 - 1, i = 1, r = 0),
        time = seq(0.1, 10000, by = 0.1),
        parms = c(beta = 0.5, gamma = 0.1)
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

# Testing Edge Cases and Model Behavior ----
test_that("No epidemic occurs when there are 0 initial infecteds", {
  # Run SIR model simulation with no initial infected individuals
  result_no_infecteds <- run_sir(
    init = c(s = 1e05 - 1, i = 0, r = 0),
    time = seq(0.1, 100, by = 0.1),
    parms = c(beta = 0.5, gamma = 0.1)
  )
  # Check that there are still no infected over time since i was zero
  expect_equal(sum(result_no_infecteds$i), 0)
})

test_that("No changes occur when all individuals start as recovered", {
  # Run SIR model simulation with everyone recovered
  result_all_recovered <- run_sir(
    init = c(s = 0, i = 0, r = 1e05),
    time = seq(0.1, 100, by = 0.1),
    parms = c(beta = 0.5, gamma = 0.1)
  )
  expect_equal(result_all_recovered$s[which.max(seq(0.1, 100, by = 0.1))], 0)
  expect_equal(result_all_recovered$i[which.max(seq(0.1, 100, by = 0.1))], 0)
  expect_equal(result_all_recovered$r[which.max(seq(0.1, 100, by = 0.1))], 1e05)
})

test_that("Model behaves reasonably with extreme parameter values", {
  # Define extreme parameters and initial conditions
  beta_extreme <- 1e-9 # Extremely low transmission rate
  gamma_extreme <- 1 # Extremely high recovery rate, very fast recovery

  # Run SIR model simulation with extreme parameters
  result_extreme <- run_sir(
    init = c(s = 1e05 - 1, i = 1, r = 0),
    time = seq(0.1, 100, by = 0.1),
    parms = c(beta = beta_extreme, gamma = gamma_extreme)
  )

  # Check if outputs remain reasonable (e.g., no negative values)
  expect_true(all(result_extreme$s >= 0))
  expect_true(all(result_extreme$i >= 0))
  expect_true(all(result_extreme$r >= 0))

  # With an extremely low beta, we expect very few new infections,
  # so the susceptible population should not decrease significantly.
  # User defined tolerance_for_new_infections on the far right
  expect_true(max((1e05 - 1) - result_extreme$s) < 1)

  # With an extremely high gamma, we expect infected individuals to recover
  # quickly, so there should be a rapid decrease in infected individuals.
  # User defined tolerance_for_recovery on the far right
  expect_true(min(diff(result_extreme$i)) < -.09)
})

# Testing analytical solutions ----
test_that("Iterative test of final size", {
  # Run SIR model simulation with no initial recovered individuals
  n <- 1e05
  i <- 1
  s <- n - i
  r <- n - s - i
  beta <- 0.5
  gamma <- 0.1
  time_seq <- seq(0.1, 150, by = 0.1)

  out <- run_sir(
    init = c(s = s, i = i, r = r),
    time = time_seq,
    parms = c(beta = beta, gamma = gamma)
  )
  # Get final size (count) from simulation
  print(tail(out))
  fs_sim <- tail(out$r, n = 1)

  # Iterative process for final size
  r_0 <- beta / gamma # same as parms above
  homog_sir_final_size_pred <- function(r_0,
                                        eps = 10^(-8),
                                        maxit = 1000000) {
    cat("Running with R0 =", r_0, "\n")
    tol <- eps^2 # for determining when to stop calculations
    a_old <- eps
    counter <- 0

    while (counter < maxit) {
      a_new <- 1 - exp(-r_0 * a_old)
      if (abs(a_old - a_new) < tol) { # stop if A doesn't change much
        break
      }
      a_old <- a_new
      counter <- counter + 1
    }

    cat("the convergence required", counter, "steps\n")
    return(a_new)
  }

  print(paste("provided R0:", r_0))
  r_0_values <- c(
    0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1,
    1.01, 1.1, 1.25, 1.5, 1.75, 2.0, r_0
  )
  print(paste("R0 values:", r_0_values))

  plot(NULL,
    xlim = range(r_0_values), ylim = c(0, 1),
    xlab = "r0", ylab = "A", main = "Plot of A vs r0"
  )

  for (r_0 in r_0_values) {
    a <- homog_sir_final_size_pred(r_0)
    cat("found A=", a, "\n\n")
    points(r_0, a, col = "blue", pch = 19)
  }
  points(r_0, homog_sir_final_size_pred(r_0), col = "red", pch = 3) # iterative
  points(r_0, fs_sim / n, col = "green", pch = 4) # simulation

  # Compare (with conversion from probability to count for fs_analytic)
  cat("Final size (simulation)", fs_sim / 1e05, "\n")
  cat("Analytical expected final size:", a, "\n")

  # Check if equal (within a certain tolerance)
  expect_equal(fs_sim / 1e05, a, tolerance = 0.00001)
})
