# Draft model to replicate the modeling work done during the Chicago Measles
# 2024 response. This is a work in progress between Catherine Herzog & Bradford
# Taylor.

# https://www.sciencedirect.com/science/article/pii/S2468266723001305?via%3Dihub#sec1 # nolint
# https://www.cdc.gov/mmwr/volumes/73/wr/mm7319a2.htm
# https://github.com/CDCgov/measles-model-chicago-2024/tree/main

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
load_all("compModels")

# Measles Models
# The approach here is to build up each piece of functionality in the
# OAW Measles Model then add as many of the modifications from the Chicago
# Measles response as possible.

# First we set up a basic SEIR model with a chain of length 2 in the
# E compartment.
# SE2IR
base_states <- c("S", "E", "I", "R")
se2ir <- define_states(base_states) |>
  add_infection("I", "S", "E", "beta") |>
  add_transition("I", "R", "tau") |>
  add_transition("E", "I", "taue", chainlength = 2)
se2ircompiled <- compilemodel(se2ir)
se2ir_rates <- se2ircompiled$modeloutstructions$processrates
se2ir_peter <- se2ircompiled$modeloutstructions$petermatrix
se2ir_states <- se2ircompiled$modeloutstructions$updatedstates

# Examples with each wrapper, going forward will stick with adaptive tau-leaping
# Both the OAW study and Chicago Measles 2024 response used adaptive tau-leaping
dyn1_g <- wrap_gillespie(
  c("S" = 999, "I" = 1, "R" = 0, "E" = 0),
  se2ircompiled,
  c(beta = 2, tau = 1, taue = .5),
  25,
  10
)
plot_stoch_model(dyn1_g, compartments = se2ir_states)

dyn1_t1 <- wrap_adaptivetau(
  c("S" = 999, "I" = 1, "R" = 0, "E" = 0),
  se2ircompiled,
  rate_func = NULL, # this defaults to compModels::generalized_rates()
  c(beta = 2, tau = 1, taue = .5),
  25,
  10
)
plot_stoch_model(dyn1_t1, compartments = se2ir_states)

dyn1_t2 <- wrap_adaptivetau(
  c("S" = 999, "I" = 1, "R" = 0, "E" = 0),
  se2ircompiled,
  rate_func = NULL, # this defaults to compModels::generalized_rates()
  c(beta = 2, tau = 1, taue = .5),
  25,
  10,
  "adaptivetau"
)
plot_stoch_model(dyn1_t2, compartments = se2ir_states)


# Now we add an intervention: vaccination
# We follow the equations in the supplement (pg 8) of Masters et al 2023 on OAW
# 𝑑𝑆𝑖(𝑡)/𝑑𝑡 = −𝛽𝐼(𝑡)𝑆𝑖(𝑡) − 𝜃(𝑡)𝑆𝑖(𝑡)
# 𝑑𝑆𝑣𝑖(𝑡)/𝑑𝑡 = 𝑝 ∗ 𝜃(𝑡)𝑆𝑣𝑖(𝑡) −𝛽𝐼(𝑡)𝑆𝑣𝑖(𝑡)
# 𝑑𝐸𝑖(𝑡)/𝑑𝑡 = 𝛽𝐼(𝑡)𝑆𝑖(𝑡) + 𝛽𝐼(𝑡)𝑆𝑣𝑖(𝑡)− 𝜎𝐸𝑖(𝑡)
# 𝑑𝐼𝑖(𝑡)/𝑑𝑡 = 𝜎𝐸𝑖 (𝑡) − 𝛾𝐼𝑖(𝑡)
# 𝑑𝑅𝑖(𝑡)/𝑑𝑡 = 𝛾𝐼𝑖(𝑡)+ (1 − 𝑝) ∗ 𝜃(𝑡)𝑆𝑖(𝑡)
# where 𝐼(𝑡) = ∑𝐼𝑖(𝑡) and the force of infection, 𝜆(𝑡) = 𝛽𝐼(𝑡)

# First we add vaccination. V represents compartment Sv, which contained
# once vaccinated individuals who remained susceptible (primary vaccine failure)
# and who could contribute to transmission but who were not vaccinated a second
# time during the initial OAW mass vaccination campaigns. Vaccine efficacy is
# represented with a proportion p.

# SE2IRV (no fork)
base_states <- c("S", "E", "I", "R", "V")
se2irv <- define_states(base_states) |>
  add_infection("I", "S", "E", "beta") |>
  add_transition("I", "R", "tau") |>
  add_transition("E", "I", "taue", chainlength = 2) |>
  add_transition("S", "V", "1/(p*theta)") |>
  # the above and below show how you indicate a rate, not a time
  add_transition("S", "R", "1/((1-p)*theta)") |>
  add_infection("I", "V", "E", "beta")

se2irvcompiled <- compilemodel(se2irv)
se2irv_rates <- se2irvcompiled$modeloutstructions$processrates
se2irv_peter <- se2irvcompiled$modeloutstructions$petermatrix
se2irv_states <- se2irvcompiled$modeloutstructions$updatedstates

dyn2 <- wrap_adaptivetau(
  c("S" = 999, "I" = 1, "R" = 0, "E" = 0, "V" = 0),
  se2irvcompiled,
  rate_func = NULL, # this defaults to compModels::generalized_rates()
  c(beta = 2, tau = 1, taue = .5, theta = .01, p = 0.4),
  25,
  10,
  "adaptivetau"
)
plot_stoch_model(dyn2, compartments = se2irv_states)


# Here we make use of the forking functionality to condense two add_transition()
# into one. This is another way to represent the model.

# SE2IRV (with fork)
se2irvf <- define_states(base_states) |>
  add_infection("I", "S", "E", "beta") |>
  add_transition("I", "R", "tau") |>
  add_transition("E", "I", "taue", chainlength = 2) |>
  add_transition("S", c("V", "R"), "1/(theta)", # indicating a rate
    forkprobability = c("p", "1-p")
  ) |>
  add_infection("I", "V", "E", "beta")

se2irvfcompiled <- compilemodel(se2irvf)
se2irvf_rates <- se2irvfcompiled$modeloutstructions$processrates
se2irvf_peter <- se2irvfcompiled$modeloutstructions$petermatrix
se2irvf_states <- se2irvfcompiled$modeloutstructions$updatedstates

dyn3 <- wrap_adaptivetau(
  c("S" = 999, "I" = 1, "R" = 0, "E" = 0, "V" = 0),
  se2irvfcompiled,
  rate_func = NULL, # this defaults to compModels::generalized_rates()
  c(beta = 2, tau = 1, taue = .5, theta = .01, p = 0.4),
  25,
  10,
  "adaptivetau"
)
plot_stoch_model(dyn3, compartments = se2irvf_states)

# Next we try to handle groups. In the original OAW model there were 5 groups
# based on eligibility for the vaccine: < 6 months, 6-11 months, 1-11 years,
# >= 12 years (not pregnant) and >= 12 years (pregnant). Those < 6 months and
# >=12 years (pregnant) did not receive MMR vaccination. Homogenous mixing was
# assumed given the shelter conditions.

# WORK IN PROGRESS
# Currently we can add groups easily, and as notes show, we are working on
# functionality to remove groups from specific functions (eg from being
# vaccinated in the forked add_transition call)

# SE2IRV with groups
se2irvf_grp <- define_states(base_states) |>
  add_infection("I", "S", "E", "beta") |>
  add_transition("I", "R", "tau") |>
  add_transition("E", "I", "taue", chainlength = 2) |>
  add_transition("S", c("V", "R"), "1/(theta)",
    forkprobability = c("p", "1-p") # indicating a rate
  ) |>
  add_infection("I", "V", "E", "beta") |>
  add_group(as.character(1:5))

se2irvfgcompiled <- compilemodel(se2irvf_grp)
se2irvfg_rates <- se2irvfgcompiled$modeloutstructions$processrates
se2irvfg_peter <- se2irvfgcompiled$modeloutstructions$petermatrix
se2irvfg_states <- se2irvfgcompiled$modeloutstructions$updatedstates

dyn4 <- wrap_adaptivetau(
  c("S" = 999, "I" = 1, "R" = 0, "E" = 0, "V" = 0),
  se2irvfgcompiled,
  rate_func = NULL, # this defaults to compModels::generalized_rates()
  c(beta = 2, tau = 1, taue = .5, theta = .01, p = 0.4),
  25,
  10,
  "adaptivetau"
)
plot_stoch_model(dyn4,
  compartments = se2irvfg_states,
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(30)
)


# Intervention timing
# Mass vaccinations were provided in specific time blocks so there is a need
# to run the model in blocks of time, stop on mass vaccination dates and adjust
# the state compartments and rates accordinging, then start again until the
# next intervention time. The results of simulations in each time block will
# need to be stitched together. Here we run 10 time steps, vaccinate, then run
# an additional 15 timesteps.

base_states <- c("S", "E", "I", "R", "Sv")
se2irvf_grp <- define_states(base_states) |>
  add_infection("I", "S", "E", "beta") |>
  add_transition("I", "R", "tau") |>
  add_transition("E", "I", "taue", chainlength = 2) |>
  add_transition("S", c("Sv", "R"), "1/(theta)",
    forkprobability = c("p", "1-p") # indicating a rate
  ) |>
  add_infection("I", "Sv", "E", "beta") |>
  add_group(as.character(1:5))
se2irvfgcompiled <- compilemodel(se2irvf_grp)
se2irvfg_rates <- se2irvfgcompiled$modeloutstructions$processrates
se2irvfg_peter <- se2irvfgcompiled$modeloutstructions$petermatrix
se2irvfg_states <- se2irvfgcompiled$modeloutstructions$updatedstates

pre_vacc <- wrap_adaptivetau(
  c("S" = 999, "I" = 1, "R" = 0, "E" = 0, "Sv" = 0),
  se2irvfgcompiled,
  rate_func = NULL, # this defaults to compModels::generalized_rates()
  c(beta = 2, tau = 1, taue = .5, theta = .01, p = 0.4),
  10,
  1,
  "adaptivetau"
) # takes ~ 1 min for 25 timesteps
last_row <- tail(pre_vacc[[1]], n = 1) |>
  select(-1) |>
  unlist()
last_row

# Vaccinating 200 of each remaining susceptible class (takes them to R, assumes
# no further primary vaccine failure)
current_states <- modify_states(last_row,
  adjust_prefix = c("S_", "R_"),
  adjust_value = c(-200, 200)
)

post_vacc <- wrap_adaptivetau(
  current_states,
  se2irvfgcompiled,
  rate_func = NULL, # this defaults to compModels::generalized_rates()
  c(beta = 2, tau = 1, taue = .5, theta = .01, p = 0.4),
  15,
  1,
  "adaptivetau",
  "current"
)

post_vacc_mod <- lapply(post_vacc, function(df) {
  df$time <- df$time + 10
  df <- df[-1, ]
  return(df)
})
dyn5 <- Map(rbind, pre_vacc, post_vacc_mod)

plot_stoch_model(dyn5,
  compartments = se2irvfg_states,
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(30)
)



# Functionality working on development branch
# Closer approximation to measles model
base_states <- c("S", "E", "I", "R", "Sv")
measlesmodel2 <- define_states(base_states) |>
  add_infection("I", "S", "E", "beta") |>
  add_transition("I", "R", "tau") |>
  add_transition("E", "I", "taue", chainlength = 2) |>
  add_transition("S", c("Sv", "R"), "1/(theta)",
    forkprobability = c("p", "1-p"),
    groupname = 3:4, processname = "vax"
  ) |>
  add_transition("S", c("Sv", "R"), "1/(theta)",
    forkprobability = c("p_infant", "1-p_infant"),
    groupname = 2, processname = "vax"
  ) |>
  add_infection("I", "Sv", "E", "beta") |>
  add_group(1:5, grouptype = "Age")
measlescompiled2 <- compilemodel(measlesmodel2)

measles_rates2 <- measlescompiled2$modeloutstructions$processrates
measles_peter2 <- measlescompiled2$modeloutstructions$petermatrix
measles_states2 <- measlescompiled2$modeloutstructions$updatedstates
measles_rates2
measles_peter2
measles_states2

# Population information
# nolint start
initial_values <- c(
  "S_group1" = 8 * .3333, "E_group1_chain1" = 0, "E_group1_chain2" = 0, "I_group1" = 0, "R_group1" = 1 - 8 * .3333, "Sv_group1" = 0,
  "S_group2" = 13 * .667, "E_group2_chain1" = 0, "E_group2_chain2" = 0, "I_group2" = 0, "R_group2" = 1 - 13 * .667, "Sv_group2" = 0,
  "S_group3" = 704 * .172, "E_group3_chain1" = 0, "E_group3_chain2" = 0, "I_group3" = 0, "R_group3" = 1 - 704 * .172, "Sv_group3" = 0,
  "S_group4" = 1139 * .1, "E_group4_chain1" = 0, "E_group4_chain2" = 0, "I_group4" = 0, "R_group4" = 1 - 1139 * .1, "Sv_group4" = 0,
  "S_group5" = 13 * .1, "E_group5_chain1" = 0, "E_group5_chain2" = 0, "I_group5" = 0, "R_group5" = 1 - 13 * .1, "Sv_group5" = 0
)
# nolint end

pre_vacc2 <- wrap_adaptivetau(
  initial_values,
  measlescompiled2,
  rate_func = NULL,
  c(
    beta = 5, tau = 5, taue = 8, theta = .01,
    p = 1 - 0.925, p_infant = 1 - .84
  ),
  15, # from date of first import to intervention day
  1,
  "adaptivetau"
)
last_row2 <- tail(pre_vacc2[[1]], n = 1) |>
  select(-1) |>
  unlist()
last_row2

current_states2 <- modify_states(last_row2,
  adjust_names = c(
    "S_group2", "S_group3", "S_group4",
    "R_group2", "R_group3", "R_group4"
  ),
  adjust_value = c(-277, -276, -329, 277, 276, 329) # actual vaccinations
)

post_vacc2 <- wrap_adaptivetau(
  current_states2,
  measlescompiled2,
  rate_func = NULL, # this defaults to compModels::generalized_rates()
  c(
    beta = 5, tau = 5, taue = 8, theta = .01, p = 1 - 0.925,
    p_infant = 1 - .84
  ),
  328, # 328 days remaining after first intervention day
  1,
  "adaptivetau",
  "current"
)

post_vacc_mod2 <- lapply(post_vacc2, function(df) {
  df$time <- df$time + 15
  df <- df[-1, ]
  return(df)
})
dyn_int2 <- Map(rbind, pre_vacc2, post_vacc_mod2)

plot_stoch_model(dyn_int2,
  compartments = measles_states2,
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(30)
)
