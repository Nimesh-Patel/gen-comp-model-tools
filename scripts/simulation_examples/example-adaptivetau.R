# Examples of using compModels::wrap_adaptivetau() with various compiled model
# structures and specifically showing us of the two methods available in
# adaptivetau package: ssa.exact() (Gillespie) and ssa.adaptivetau()

library(devtools)
library(Matrix) # nolint object_name_linter
library(dplyr)
library(tibble)
library(tidyr)
library(tidyselect)
library(purrr)
library(stringr)
library(deSolve) # nolint object_name_linter
library(ggplot2)
# load package w/o installing.
load_all("compModels")

# Loading a modified wrap_adaptivetau function to enable run of any model below
# in which base states get modified (eg chains, groups, metapopulations)
# This may get incorporated into compModels in the near future

wrap_adaptivetau <- function(init_vals, compiledmodel, rate_func = NULL,
                             parameters, n_timesteps, n_sims,
                             method = "exact") {
  # rateeqns, transitions, parameters
  model_rates <- as.list(compiledmodel$modeloutstructions$processrates)
  model_peter <- as.matrix(compiledmodel$modeloutstructions$petermatrix)
  parameters <- parameters

  # Ensure init_vals matches model dimensions
  if (length(init_vals) != nrow(model_peter)) {
    x0_tbl <- define_initialstate(compiledmodel, init_vals)
    x0 <- x0_tbl$X0
    names(x0) <- x0_tbl$updatedstate
  } else {
    x0 <- init_vals
    if (is.null(names(x0))) {
      names(x0) <- compiledmodel$modeloutstructions$updatedstates
    }
  }

  t <- n_timesteps

  if (is.null(rate_func)) {
    rate_func <- generalized_rates
  }

  # Prepare rate function with named state vector
  # adaptivetau may strip names from the state vector passed to the rate fn
  # Workaround: use wrap rate fn function to reapply names based on the model
  # so that generalized_rates() can find the variables in the environment
  state_names <- compiledmodel$modeloutstructions$updatedstates
  actual_rate_func <- rate_func(model_rates)
  wrapped_rate_func <- function(state, params, t) {
    names(state) <- state_names
    actual_rate_func(state, params, t)
  }

  if (method == "exact") {
    sims <- lapply(1:n_sims, function(i) {
      adaptivetau::ssa.exact(
        init.values = x0,
        transitions = model_peter,
        rateFunc = wrapped_rate_func,
        params = parameters,
        tf = t
      ) |> data.frame()
    })
  } else if (method == "adaptivetau") {
    sims <- lapply(1:n_sims, function(i) {
      adaptivetau::ssa.adaptivetau(
        init.values = x0,
        transitions = model_peter,
        rateFunc = wrapped_rate_func,
        params = parameters,
        tf = t
      ) |> data.frame()
    })
  } else {
    stop("Method not recognized. Please use 'exact' or 'adaptivetau'")
  }
}


# SIR instructions
base_states <- c("S", "I", "R")
sir <- define_states(base_states) |>
  add_infection("I", "S", "I", "beta") |>
  add_transition("I", "R", "tau")
sircompiled <- compilemodel(sir)
sir_rates <- sircompiled$modeloutstructions$processrates
sir_peter <- sircompiled$modeloutstructions$petermatrix
sir_states <- sircompiled$modeloutstructions$updatedstates

# SIR
init_vals <- c("S" = 999, "I" = 1, "R" = 0)
parameters <- c(beta = 2, tau = 1)

# This function is needed to reapply state names within the rate function for
# the use of adaptivetau package functions
named_rate_func <- function(rates, state_names) {
  actual_rate_func <- generalized_rates(as.list(rates))
  function(state, params, t) {
    names(state) <- state_names
    actual_rate_func(state, params, t)
  }
}

# adaptivetau::ssa.exact - Gillespie
dyn <- adaptivetau::ssa.exact(
  init.values = init_vals,
  transitions = as.matrix(sir_peter),
  rateFunc = named_rate_func(sir_rates, sir_states),
  params = parameters,
  tf = 100
)

# adaptivetau::adaptivetau - adaptivetau leaping
dynt <- adaptivetau::ssa.adaptivetau(
  init.values = init_vals,
  transitions = as.matrix(sir_peter),
  rateFunc = named_rate_func(sir_rates, sir_states),
  params = parameters,
  tf = 100
)

# ssa.exact
dynplt <- tibble::as_tibble(dyn) |>
  dplyr::rename(t = time) |>
  tidyr::pivot_longer(!t) |>
  dplyr::mutate(density = as.numeric(value), t = as.numeric(t))

ggplot2::ggplot(dynplt, aes(x = t, y = density, color = name)) +
  ggplot2::geom_line(linewidth = 2) +
  ggplot2::theme_classic(base_size = 16) +
  ggplot2::labs(
    title = "SIR Model Simulations \n SSA Exact",
    x = "Time",
    y = "Number of Individuals"
  )

# adaptivetau
dyntplt <- tibble::as_tibble(dynt) |>
  dplyr::rename(t = time) |>
  tidyr::pivot_longer(!t) |>
  dplyr::mutate(density = as.numeric(value), t = as.numeric(t))

ggplot2::ggplot(dyntplt, aes(x = t, y = density, color = name)) +
  ggplot2::geom_line(linewidth = 2) +
  ggplot2::theme_classic(base_size = 16) +
  ggplot2::labs(
    title = "SIR Model Simulations \n SSA Adaptivetau",
    x = "Time",
    y = "Number of Individuals"
  )

# Running with wrap_adaptivetau for multiple sims for Gillespie and tau-leaping
# algorithms, respectively
dyn2 <- wrap_adaptivetau(init_vals, sircompiled,
  rate_func = NULL,
  parameters, 25, 10
)
plot_stoch_model(dyn2, compartments = sir_states)
dyn2t <- wrap_adaptivetau(init_vals, sircompiled,
  rate_func = NULL,
  parameters, 25, 10, "adaptivetau"
)
plot_stoch_model(dyn2, compartments = sir_states)

# SEIR
base_states <- c("S", "I", "R", "E")
seir <- define_states(base_states) |>
  add_infection("I", "S", "E", "beta") |>
  add_transition("I", "R", "tau") |>
  add_transition("E", "I", "taue")
seircompiled <- compilemodel(seir)
seir_rates <- seircompiled$modeloutstructions$processrates
seir_peter <- seircompiled$modeloutstructions$petermatrix
seir_states <- seircompiled$modeloutstructions$updatedstates

init_vals <- c("S" = 999, "I" = 1, "R" = 0, "E" = 0)
parameters <- c(beta = 2, tau = 1, taue = .5)

# Exact method
dyn3 <- wrap_adaptivetau(init_vals, seircompiled,
  rate_func = NULL,
  parameters, 25, 10
)
plot_stoch_model(dyn3, compartments = seir_states)
# Tau method
dyn3t <- wrap_adaptivetau(init_vals, seircompiled,
  rate_func = NULL,
  parameters, 25, 10, "adaptivetau"
)
plot_stoch_model(dyn3t, compartments = seir_states)


# SEI2R
base_states <- c("S", "I", "R", "E")
sei2r <- define_states(base_states) |>
  add_infection("I", "S", "E", "beta") |>
  add_transition("I", "R", "tau", chainlength = 2) |>
  add_transition("E", "I", "taue")
sei2rcompiled <- compilemodel(sei2r)
sei2r_rates <- sei2rcompiled$modeloutstructions$processrates
sei2r_peter <- sei2rcompiled$modeloutstructions$petermatrix
sei2r_states <- sei2rcompiled$modeloutstructions$updatedstates

init_vals <- c("S" = 999, "I" = 1, "R" = 0, "E" = 0)
parameters <- c(beta = 2, tau = 1, taue = .5)

dyn4 <- wrap_adaptivetau(init_vals, sei2rcompiled,
  rate_func = NULL,
  parameters, 25, 10
)
plot_stoch_model(dyn4, compartments = sei2r_states)
dyn4t <- wrap_adaptivetau(init_vals, sei2rcompiled,
  rate_func = NULL,
  parameters, 25, 10, "adaptivetau"
)
plot_stoch_model(dyn4t, compartments = sei2r_states)


# SEI2{RH}
base_states <- c("S", "E", "I", "R", "H")
sei2rh <- define_states(base_states) |>
  add_infection("I", "S", "E", "beta") |>
  add_transition("I", c("R", "H"), "tau",
    chainlength = 2, forkprobability = c("1-probH", "probH")
  ) |>
  add_transition("E", "I", "taue")
sei2rhcompiled <- compilemodel(sei2rh)
sei2rh_rates <- sei2rhcompiled$modeloutstructions$processrates
sei2rh_peter <- sei2rhcompiled$modeloutstructions$petermatrix
sei2rh_states <- sei2rhcompiled$modeloutstructions$updatedstates

init_vals <- c("S" = 999, "E" = 0, "I" = 1, "R" = 0, "H" = 0)
parameters <- c(beta = 2, tau = 1, taue = .5, probH = .2)

dyn5 <- wrap_adaptivetau(init_vals, sei2rhcompiled,
  rate_func = NULL,
  parameters, 25, 10
)
plot_stoch_model(dyn5, compartments = sei2rh_states)
dyn5t <- wrap_adaptivetau(init_vals, sei2rhcompiled,
  rate_func = NULL,
  parameters, 25, 10, "adaptivetau"
)
plot_stoch_model(dyn5t, compartments = sei2rh_states)


# SEI2R and Treatment
base_states <- c("S", "E", "I", "R")
sei2r_treat <- define_states(base_states) |>
  add_infection("I", "S", "E", "beta") |>
  add_transition("I", "R", "tau") |>
  add_transition("E", "I", "taue") |>
  add_transition("I", "R", "taut")
sei2r_treatcompiled <- compilemodel(sei2r_treat)
sei2r_treat_rates <- sei2r_treatcompiled$modeloutstructions$processrates
sei2r_treat_peter <- sei2r_treatcompiled$modeloutstructions$petermatrix
sei2r_treat_states <- sei2r_treatcompiled$modeloutstructions$updatedstates

init_vals <- c("S" = 999, "I" = 1, "R" = 0, "E" = 0)
parameters <- c(beta = 2, tau = 1, taue = .5, taut = 2)

dyn6 <- wrap_adaptivetau(init_vals, sei2r_treatcompiled,
  rate_func = NULL,
  parameters, 25, 10
) # try end time = 50 to see more
plot_stoch_model(dyn6, compartments = sei2r_treat_states)
dyn6t <- wrap_adaptivetau(init_vals, sei2r_treatcompiled,
  rate_func = NULL,
  parameters, 25, 10, "adaptivetau"
)
plot_stoch_model(dyn6t, compartments = sei2r_treat_states)


# Metapopulations
sirmeta <- define_states(c("S", "I", "R")) |>
  add_infection("I", "S", "I", "beta") |>
  add_transition("I", "R", "tau") |>
  define_metapopulations(
    metapopulation = c("R0", "R0times2"),
    scaletransitions = c(1, 2)
  ) |>
  add_travel("mu")
outmeta <- compilemodel(sirmeta)
outmeta_rates <- outmeta$modeloutstructions$processrates
outmeta_peter <- outmeta$modeloutstructions$petermatrix
outmeta_states <- outmeta$modeloutstructions$updatedstates

init_vals <- c("S" = 999, "I" = 1, "R" = 0)
parameters <- c(beta = 2, tau = 1, mu = .1)

dyn7 <- wrap_adaptivetau(init_vals, outmeta,
  rate_func = NULL,
  parameters, 25, 10
)
plot_stoch_model(dyn7, compartments = outmeta_states)
dyn7t <- wrap_adaptivetau(init_vals, outmeta,
  rate_func = NULL,
  parameters, 25, 10, "adaptivetau"
)
plot_stoch_model(dyn7t, compartments = outmeta_states)


# Incorporating Groups
# One Group
sir1group <- define_states(c("S", "I", "R")) |>
  add_infection("I", "S", "I", "beta") |>
  add_transition("I", "R", "tau") |>
  add_group(
    groupname = c("Social", "Antisocial"),
    scaleprocessbyname = list(infection = c(2, 1))
  )
sir1groupcomp <- compilemodel(sir1group)
sir1groupcomp_rates <- sir1groupcomp$modeloutstructions$processrates
sir1groupcomp_peter <- sir1groupcomp$modeloutstructions$petermatrix
sir1groupcomp_states <- sir1groupcomp$modeloutstructions$updatedstates

init_vals <- c("S" = 999, "I" = 1, "R" = 0)
parameters <- c(beta = 1.5, tau = 1)

dyn8 <- wrap_adaptivetau(init_vals, sir1groupcomp,
  rate_func = NULL,
  parameters, 25, 10
)
plot_stoch_model(dyn8, compartments = sir1groupcomp_states)
dyn8t <- wrap_adaptivetau(init_vals, sir1groupcomp,
  rate_func = NULL,
  parameters, 25, 10, "adaptivetau"
)
plot_stoch_model(dyn8t, compartments = sir1groupcomp_states)


# Two Groups
sir2group <- define_states(c("S", "I", "R")) |>
  add_infection("I", "S", "I", "beta") |>
  add_transition("I", "R", "tau") |>
  add_group(
    groupname = c("Young", "Old"),
    grouptype = "Age",
    scaleprocessbyname = list(infection = c(2, 1))
  ) |>
  add_group(
    groupname = c("Patient", "HCW"),
    grouptype = "Hospital",
    scaleprocessbyname = list(infection = c(3, 4))
  ) |>
  compilemodel()
sir2group_rates <- sir2group$modeloutstructions$processrates
sir2group_peter <- sir2group$modeloutstructions$petermatrix
sir2group_states <- sir2group$modeloutstructions$updatedstates

init_vals <- c("S" = 999, "I" = 1, "R" = 0)
parameters <- c(beta = 2, tau = 1)

dyn9 <- wrap_adaptivetau(init_vals, sir2group,
  rate_func = NULL,
  parameters, 25, 10
)
plot_stoch_model(dyn9,
  compartments = sir2group_states,
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(12)
)
dyn9t <- wrap_adaptivetau(init_vals, sir2group,
  rate_func = NULL,
  parameters, 25, 10, "adaptivetau"
)
plot_stoch_model(dyn9t,
  compartments = sir2group_states,
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(12)
)


# Full model
sirfull <- define_states(c("S", "I", "R")) |>
  define_metapopulations(metapopulation = c("UK", "USA")) |>
  add_travel("mu") |>
  add_group(c("Young", "Old"), grouptype = "Age") |>
  add_group(c("Patient", "HCW"), grouptype = "Hospital") |>
  add_infection("I", "S", "I", "beta") |>
  add_transition("I", "R", "tau", chainlength = 3)

outlistfull <- compilemodel(sirfull)
outlistfull_rates <- outlistfull$modeloutstructions$processrates
outlistfull_peter <- outlistfull$modeloutstructions$petermatrix
outlistfull_states <- outlistfull$modeloutstructions$updatedstates

init_vals <- c("S" = 999, "I" = 1, "R" = 0)
parameters <- c(beta = 2, tau = 1, mu = 0.1) # arbitrarily selected mu

# This model with 10 sims could take some time to run, 1 sim took 3-4 min
# Recommend popping out plots to view in html browser, due to large legend size
dyn10 <- wrap_adaptivetau(init_vals, outlistfull,
  rate_func = NULL,
  parameters, 25, 10
)
plot_stoch_model(dyn10,
  compartments = outlistfull_states,
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(40)
)
# This model with 10 sims could run ~ 3-4 min
dyn10t <- wrap_adaptivetau(init_vals, outlistfull,
  rate_func = NULL,
  parameters, 25, 10, "adaptivetau"
)
plot_stoch_model(dyn10t,
  compartments = outlistfull_states,
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(40)
)
