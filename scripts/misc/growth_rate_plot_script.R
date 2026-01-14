library(ggplot2)
library(dplyr)
library(tidyr)
library(rlang)

# Loading compModels functions

validate_sir_input <- function(init, time, parms, intervention = FALSE) {
  if (!is.numeric(parms[["beta"]])) stop("beta must be numeric")
  if (!is.numeric(parms[["gamma"]])) stop("gamma must be numeric")
  if (!is.numeric(init["s"])) stop("s must be numeric")
  if (!is.numeric(init["i"])) stop("i must be numeric")
  if (!is.numeric(init["r"])) stop("r must be numeric")
  if (parms[["beta"]] < 0) stop("beta should be non-negative")
  if (parms[["gamma"]] < 0) stop("gamma should be non-negative")
  if (init["s"] < 0) stop("s should be non-negative")
  if (init["i"] < 0) stop("i should be non-negative")
  if (init["r"] < 0) stop("r should be non-negative")
  if (length(time) == 0) {
    stop("time vector must have at least one element")
  }
  if (!all(diff(time) > 0)) {
    stop("time vector must be strictly increasing")
  }
  if (any(c(init["s"], init["i"], init["r"]) > .Machine$integer.max)) {
    stop("Initial conditions exceed the maximum integer size.")
  }
  print("Input checks passed")
}


sir_deriv <- function(time, state, parms) {
  s <- state[1]
  i <- state[2]
  r <- state[3]
  n <- s + i + r
  ds <- (-parms["beta"] * s * i) / n
  di <- ((parms["beta"] * s * i) / n) - parms["gamma"] * i
  dr <- parms["gamma"] * i
  list(c(s = ds, i = di, r = dr))
}


run_sir <- function(init, time, parms) {
  validate_sir_input(init, time, parms)

  sir_out <- deSolve::ode(
    y = init, times = time,
    func = sir_deriv,
    parms = parms
  )
  as.data.frame(sir_out)
}


# Example usage with deterministic model:
out <- run_sir(
  init = c(s = 1e05 - 1, i = 1, r = 0),
  time = seq(0.1, 100, by = 0.1),
  parms = c(beta = 0.5, gamma = 0.1)
)

apply_growth_rate_calc <- function(output, time_var = NULL,
                                   i_var = "i",
                                   beta = NULL, gamma = NULL,
                                   lower_time = NULL,
                                   upper_time = NULL) {
  output$inc <- output[[i_var]]
  output$log_inc <- log(output[[i_var]])
  max_log_inc <- max(output$log_inc)
  names(output)[names(output) == time_var] <- "time"


  formula_string <- paste("log_inc ~ time")
  # This line below needed for stochastic applications
  output <- output %>%
    filter(!is.na(log_inc) & !is.nan(log_inc) & is.finite(log_inc))
  output$time <- as.numeric(output$time)
  output2 <- output |> filter(between(time, lower_time, upper_time))

  fit <- lm(as.formula(formula_string), data = output2)
  r <- fit$coef[["time"]]
  print(paste("Growth rate is:", r))

  i_0 <- output$log_inc[1]
  print(paste("log(i_0) is:", i_0))
  stoch_t <- output$time
  stoch_t_f <- stoch_t[stoch_t >= lower_time & stoch_t <= upper_time]

  theorygrowth <-
    function(beta, gamma, t, i0) {
      outdyn <- exp(i0) * exp((beta - gamma) * t)
    }
  predicted_inc <- theorygrowth(beta, gamma, stoch_t_f, i_0)

  r_growth <-
    function(r, t, i0) {
      outdyn <- exp(i0) * exp(r * t)
    }
  second_r_line <- r_growth(r, stoch_t_f, i_0)

  output3 <- output |>
    select(time, log_inc) |>
    filter(between(time, lower_time, upper_time))

  exp_growth_df <- data.frame(
    time = stoch_t_f,
    pred = predicted_inc,
    rline = log(second_r_line), # theory using calculated r
    log_pred_inc = log(predicted_inc), # theory
    log_inc = output3$log_inc # observations
  )
  str(exp_growth_df)

  exp_growth_df2 <- exp_growth_df |>
    pivot_longer(!time,
      names_to = "group",
      values_to = "value"
    )
  str(exp_growth_df2)


  p <- ggplot(
    data = subset(exp_growth_df2, !group %in% c("log_inc", "pred")),
    aes(x = time, y = value, color = group)
  ) +
    geom_point(data = output, aes(x = time, y = log_inc, color = "log_inc")) +
    geom_line() +
    labs(
      title = "Log-Transformed Incidence Over Time with
                  Exponential Growth Curve",
      y = "Log(Incidence)",
      color = "Legend"
    ) +
    scale_color_manual(
      labels = c(
        "Simulated Observations",
        "Theoretical Expected Growth",
        "Expected Growth (Calculated r)"
      ),
      values = c(
        log_inc = "grey",
        log_pred_inc = "red",
        rline = "blue"
      )
    ) +
    theme_bw() +
    theme(legend.position = "top") +
    scale_x_continuous(limits = c(
      min(output$time),
      max(output$time)
    ))
}

calculate_initial_growth_rate <- function(output, time_var = NULL,
                                          i_var = "i",
                                          beta = NULL, gamma = NULL,
                                          lower_time = NULL,
                                          upper_time = NULL) {
  if (is.data.frame(output)) {
    if (!i_var %in% colnames(output)) {
      stop(paste("Column", i_var, "not found in the output data frame."))
    }
    if (!time_var %in% colnames(output)) {
      stop(paste("Column", time_var, "not found in the output data frame."))
    }
    print("output is a single data frame")
    return(apply_growth_rate_calc(
      output, time_var, i_var, beta, gamma,
      lower_time, upper_time
    ))
  } else if (is.list(output)) {
    print("output is a list of data frames")
    lapply(output, function(df) {
      apply_growth_rate_calc(
        df, time_var, i_var, beta, gamma,
        lower_time, upper_time
      )
    })
  } else {
    stop("Invalid input: output must be a data frame or a list of data frames.")
  }
}


# Deterministic call
calculate_initial_growth_rate(out,
  time_var = "time",
  beta = 0.5, gamma = 0.1,
  lower_time = 10,
  upper_time = 25
)




# Example with stochastic model
run_gen_stochastic <- function(parms_vec, propensity_fns, init_vals,
                               n_timesteps, change_matrix, n_sims) {
  validate_gen_stoch_input(parms_vec, init_vals, n_timesteps,
    method = "GillespieSSA",
    propensity_fns = propensity_fns,
    change_matrix = change_matrix,
    transitions = NULL,
    n_sims,
    intervention_start_time = NULL,
    intervention_end_time = NULL,
    modifier = NULL
  )

  nu <- change_matrix # State Change Matrix
  t <- n_timesteps
  parms_vec <- c(parms_vec, n = sum(init_vals))

  sims <- lapply(1:n_sims, function(i) {
    GillespieSSA::ssa(
      x0 = init_vals,
      a = propensity_fns,
      nu = nu,
      parms = parms_vec,
      tf = t,
      method = GillespieSSA::ssa.d(),
      simName = "General Stochastic Model"
    )$data
  })

  lapply(sims, function(x) data.frame(x)) # output list object of data frames
}

validate_gen_stoch_input <- function(parms_vec, init_vals,
                                     n_timesteps,
                                     method = NULL,
                                     propensity_fns = NULL,
                                     change_matrix = NULL,
                                     transitions = NULL,
                                     n_sims,
                                     intervention_start_time = NULL,
                                     intervention_end_time = NULL,
                                     modifier = NULL) {
  check_parms_vec(parms_vec)
  check_init_vals(init_vals)
  check_time_step(n_timesteps)
  if (method == "GillespieSSA") {
    check_gillespie_args(propensity_fns, change_matrix)
  } else if (method == "adaptivetau") {
    check_adaptivetau_args(transitions)
  }
  check_nsims(n_sims)
  check_intervention_times(intervention_start_time, intervention_end_time)
  check_modifier2(modifier)

  print("All input checks passed")
}

check_parms_vec <- function(parms_vec) {
  if (!is.numeric(parms_vec)) {
    stop("Parameter values must be numeric.")
  }

  if (any(parms_vec <= 0)) {
    stop("All parameter values must be positive.")
  }

  print("parameter values checks passed")
}

check_init_vals <- function(init_vals) {
  if (!all(is.numeric(init_vals))) {
    stop("All initial values must be numeric.")
  }

  if (any(init_vals < 0)) {
    stop("All initial values must be non-negative.")
  }

  if (any(init_vals > .Machine$integer.max)) {
    stop("Initial conditions exceed the maximum integer size.")
  }

  print("initial values checks passed")
}

check_time_step <- function(n_timesteps) {
  if (!is.numeric(
    n_timesteps
  ) ||
    length(n_timesteps) != 1 ||
    n_timesteps <= 0) {
    stop("The number of time steps must be a single positive numeric value.")
  }

  print("times values checks passed")
}

check_gillespie_args <- function(propensity_fns, change_matrix) {
  if (!is.character(propensity_fns)) stop("propensity_fns must be character")
  if (!is.matrix(change_matrix)) stop("change_matrix must be a matrix")

  print("Gillespie argument value checks passed")
}

check_adaptivetau_args <- function(transitions) {
  if (!is.list(transitions)) stop("transitions must be a list")
  sapply(
    X = transitions,
    FUN = function(x) {
      if (!is.numeric(x)) {
        stop("All elements of transitions list must be numeric")
      }
    }
  )
  sapply(
    X = transitions,
    FUN = function(x) {
      if (!all(!is.null(names(x)))) {
        stop("All elements of transitions list must be named numeric vectors")
      }
    }
  )

  print("adaptivetau argument value checks passed")
}

check_nsims <- function(n_sims) {
  if (!is.numeric(n_sims) || length(n_sims) != 1 || n_sims <= 0) {
    stop("The number of time steps must be a single positive numeric value.")
  }

  print("n_sims values checks passed")
}

check_intervention_times <- function(intervention_start_time,
                                     intervention_end_time) {
  if (!is.null(
    intervention_start_time
  ) &&
    !(
      is.numeric(intervention_start_time) &&
        length(intervention_start_time) == 1
    )) {
    stop("intervention_start_time must be a single numeric value or NULL.")
  }

  if (!is.null(
    intervention_end_time
  ) &&
    !(
      is.numeric(intervention_end_time) &&
        length(intervention_end_time) == 1
    )) {
    stop("intervention_end_time must be a single numeric value or NULL.")
  }

  print("intervention start and end time checks passed")
}

check_modifier1 <- function(modifier) {
  is_null <- is.null(modifier)
  is_single_numeric <- is.numeric(modifier) && length(modifier) == 1
  is_numeric_matrix <- is.matrix(modifier) && is.numeric(modifier)
  is_function <- is.function(modifier)

  return(is_null || is_single_numeric || is_numeric_matrix || is_function)
  print("modifier check 1 passed")
}


check_modifier2 <- function(modifier) {
  if (!check_modifier1(modifier)) {
    stop("Modifier must be either a single numeric value, a matrix of numeric
          values, or a user-defined function.")
  }

  print("modifier checks 1 and 2 passed")
}

init_vals <- c(s = 499, i = 1, r = 0)
change_matrix <- matrix(
  c(
    -1, 0, # S -> I #nolint
    +1, -1, # I -> R #nolint
    0, +1
  ),
  nrow = length(init_vals),
  byrow = TRUE,
  dimnames = list(names(init_vals), NULL)
)

modelout <- run_gen_stochastic(
  parms_vec = c(
    beta = 0.5, gamma = 0.1
  ),
  propensity_fns = c(
    "(beta * s * i)/n",
    "gamma * i"
  ),
  init_vals = init_vals,
  n_timesteps = 100,
  change_matrix = change_matrix,
  n_sims = 10
)
plot_stoch_model(modelout, compartments = c("s", "i", "r"))

# Stochastic Call
calculate_initial_growth_rate(modelout,
  time_var = "t", i_var = "i",
  beta = 0.5, gamma = 0.1,
  lower_time = 0,
  upper_time = 25
)
# If there are no observations in the user specified time range this will
# give an R error

# Trying one df in the list
calculate_initial_growth_rate(modelout[[2]],
  time_var = "t",
  beta = 0.5, gamma = 0.1,
  lower_time = 10,
  upper_time = 25
)
