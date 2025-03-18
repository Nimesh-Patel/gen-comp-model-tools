# Examples of adaptivetau:ssa.exact() and adaptivetau::ssa.adaptivetau() with
# the SIRModelBuilder Approach

library(devtools)
library(Matrix)
library(dplyr)
library(tibble)
library(tidyr)
library(tidyselect)
library(purrr)
library(stringr)
library(deSolve)
library(ggplot2)
# load package w/o installing.
load_all("compModels")


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

# adaptivetau::ssa.exact - Gillespie
dyn <- adaptivetau::ssa.exact(
  init.values = init_vals,
  transitions = as.matrix(sir_peter),
  rateFunc = generalized_rates(as.list(sir_rates)),
  params = parameters,
  tf = 100
)

# adaptivetau::adaptivetau - adaptivetau leaping
dynt <- adaptivetau::ssa.adaptivetau(
  init.values = init_vals,
  transitions = as.matrix(sir_peter),
  rateFunc = generalized_rates(as.list(sir_rates)),
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


# Running with wrap_adaptivetau for multiple sims
dyn2 <- wrap_adaptivetau(x0, sircompiled, rate_func = NULL, parameters, 25, 10)
plot_stoch_model(dyn2, compartments = sir_states)


# SEIR
base_states <- c("S", "E", "I", "R")
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

dyn3 <- wrap_adaptivetau(init_vals, seircompiled,
  rate_func = NULL,
  parameters, 25, 10
)
plot_stoch_model(dyn3, compartments = seir_states)

dyn3t <- wrap_adaptivetau(init_vals, seircompiled,
  rate_func = NULL,
  parameters, 25, 10, "adaptivetau"
)
plot_stoch_model(dyn3t, compartments = seir_states)


# SEI2R
base_states <- c("S", "E", "I", "R")
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

init_vals <- c("S" = 999, "I" = 1, "R" = 0, "E" = 0, "H" = 0)
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
    metapop_names = c("R0", "R0times2"),
    interactionscale = c(1, 2)
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
  add_group(c("Social", "Antisocial"), interactionscale = c(2, 1))
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
  add_group(c("Young", "Old"),
    grouptype = "Age",
    interactionscale = c(2, 1)
  ) |>
  add_group(c("Patient", "HCW"),
    grouptype = "Hospital",
    interactionscale = c(3, 4)
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
  define_metapopulations(metapop_names = c("UK", "USA")) |>
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
