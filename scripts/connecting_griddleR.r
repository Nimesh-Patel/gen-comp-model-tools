library(devtools)
load_all("compModels")
library(griddleR)
library(tibble)
library(dplyr)
library(ggplot2)


# SIR model without griddleR
base_states <- c("S", "I", "R")
sir <- define_states(base_states) |>
  add_infection("I", "S", "I", "beta") |>
  add_transition("I", "R", "tau")
sircompiled <- compilemodel(sir)
sir_rates <- sircompiled$modeloutstructions$processrates
sir_peter <- sircompiled$modeloutstructions$petermatrix
sir_states <- sircompiled$modeloutstructions$updatedstates

init_vals <- c("S" = 999, "I" = 1, "R" = 0)
parameters <- c(beta = 2, tau = 1)
# Deterministic SIR ODE Model
dyn1 <- wrap_ode(
  c("S" = 999, "I" = 1, "R" = 0),
  seq(0, 25, by = 0.01),
  sircompiled,
  c(beta = 2, tau = 1)
)
plot_determin_model(dyn1)

# Stochastic SIR model
dyn2 <- wrap_adaptivetau(init_vals, sircompiled,
  rate_func = NULL,
  parameters, 25, 10
)
plot_stoch_model(dyn2, compartments = sir_states)


# Bringing in configuration file for use with griddleR
readLines("scripts/sir_params.yaml") |>
  paste(collapse = "\n") |>
  cat()
parameter_sets <- griddleR::read_griddle("scripts/sir_params.yaml")
cat(yaml::as.yaml(parameter_sets))


# Modified run() function from griddleR
run <- function(fun, parameter_sets, ...) {
  parameter_hashes <- purrr::map_chr(parameter_sets, rlang::hash)
  # named list: hash -> parameter set
  parameter_map <- rlang::set_names(parameter_sets, parameter_hashes)

  dots <- list(...)

  simulations <- furrr::future_map(
    seq_along(parameter_sets),
    function(i) {
      # Explicitly call fun with parameter_set and the additional arguments
      args <- c(list(parameter_sets[[i]]), dots)
      list(
        parameter_hash = parameter_hashes[i],
        result = do.call(fun, args)
      )
    },
    # Ensure independent rng for stochastic runs
    # (or could set seed based on parameters)
    .options = furrr::furrr_options(seed = TRUE)
  )

  list(
    parameter_map = parameter_map,
    simulations = simulations
  )
}

# Deterministic SIR ODE model with griddleR
# Note: will need a more flexible way to better handle the initial compartment
# values that is less hard coded.
double_wrap_ode <- function(parameter_set,
                            compiledmodel) {
  y <- c(S = parameter_set$S0, I = parameter_set$I0, R = parameter_set$R0)
  start_time <- parameter_set$start_time
  max_time <- parameter_set$n_timesteps
  times_by <- parameter_set$timesby
  n_sims <- parameter_set$n_sims

  # arguments for wrap_ode are:
  # wrap_ode <- function(y, times, compiledmodel, parms) {
  result <- wrap_ode(
    y = y,
    times = seq(start_time, max_time, by = times_by),
    compiledmodel = compiledmodel,
    parms = parameter_set
  )
  result |>
    as_tibble() |>
    mutate(across(everything(), as.numeric))
}

output <- run(
  fun = double_wrap_ode,
  parameter_sets = parameter_sets,
  compiledmodel = sircompiled
)

str(output)
griddleR::flatten_run(output)
griddleR::flatten_run(output) |>
  dplyr::mutate(across(c(beta, tau), factor)) |>
  ggplot2::ggplot(aes(time, I)) +
  # facet_wrap(vars(beta), labeller = label_both) +
  ggplot2::geom_line(aes(color = beta, linetype = tau)) +
  ggplot2::ggtitle("SIR ODE for varying beta (transmission) & tau (recovery)") +
  theme_classic()


# Stochastic SIR ODE model with griddleR
# Note: will need a more flexible way to better handle the initial compartment
# values and the parameters that is less hard coded.
double_wrap_adaptivetau <- function(parameter_set,
                                    compiledmodel,
                                    rate_func = NULL) {
  init_vals <- c(
    S = parameter_set$S0, I = parameter_set$I0,
    R = parameter_set$R0
  )
  parameters <- c(beta = parameter_set$beta, tau = parameter_set$tau)
  n_timesteps <- parameter_set$n_timesteps
  n_sims <- parameter_set$n_sims

  # arguments for wrap_adaptivetau are:
  # wrap_adaptivetau <- function(init_vals, compiledmodel, rate_func = NULL,
  #                              parameters, n_timesteps, n_sims,
  #                              method = "exact", usestates = "initial")
  result <- wrap_adaptivetau(
    init_vals = init_vals,
    compiledmodel = compiledmodel,
    rate_func = rate_func,
    parameters = parameters,
    n_timesteps = n_timesteps,
    n_sims = n_sims
  )

  result_tibbles <- lapply(seq_along(result), function(i) {
    tibble::as_tibble(result[[i]]) %>%
      dplyr::mutate(across(everything(), as.numeric)) %>%
      dplyr::mutate(simulation = i)
  })

  dplyr::bind_rows(result_tibbles)
}

output2 <- run(
  fun = double_wrap_adaptivetau,
  parameter_sets = parameter_sets,
  compiledmodel = sircompiled
)

str(output2)
griddleR::flatten_run(output2)
griddleR::flatten_run(output2) |>
  dplyr::mutate(across(c(beta, tau), factor)) |>
  ggplot2::ggplot(aes(time, I)) +
  # facet_wrap(vars(beta), labeller = label_both) +
  ggplot2::geom_line(aes(color = beta, linetype = tau), alpha = 0.4) +
  ggplot2::ggtitle("Stochastic SIR for varying transmission & recovery rates") +
  ggplot2::theme_classic() +
  ggplot2::labs(
    x = "Time",
    y = "Infected Individuals",
    color = "Transmission Rate (β)",
    linetype = "Recovery Rate (τ)"
  )


griddleR::flatten_run(output2) |>
  dplyr::mutate(across(c(beta, tau, simulation), factor)) |>
  ggplot2::ggplot(aes(time, I, group = interaction(beta, tau, simulation))) +
  ggplot2::geom_line(aes(color = beta, linetype = tau), alpha = 0.4) +
  ggplot2::ggtitle("Stochastic SIR: Simulations by Parameter Set") +
  ggplot2::facet_grid(tau ~ beta, labeller = label_both) +
  ggplot2::theme_classic() +
  ggplot2::labs(
    x = "Time",
    y = "Infected Individuals",
    color = "Transmission Rate (β)",
    linetype = "Recovery Rate (τ)"
  )
