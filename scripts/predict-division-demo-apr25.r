# Measles model draft replicate the model development and simulation work done
# during the Chicago Measles 2024 response.
# Catherine Herzog & Bradford Taylor.

# Background
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

# Measles Models Details:
# Both the OAW study and Chicago Measles 2024 response used adaptive tau-leaping
# Intervention:  vaccination
# We follow the equations in the supplement (pg 8) of Masters et al 2023 on OAW
# compartment Sv contains once vaccinated individuals who remained susceptible
# (primary vaccine failure) and who could contribute to transmission but who
# were not vaccinated a second time during the initial OAW mass vaccination
# campaigns.
# 𝑑𝑆𝑖(𝑡)/𝑑𝑡 = −𝛽𝐼(𝑡)𝑆𝑖(𝑡) − 𝜃(𝑡)𝑆𝑖(𝑡)
# 𝑑𝑆𝑣𝑖(𝑡)/𝑑𝑡 = 𝑝 ∗ 𝜃(𝑡)𝑆𝑣𝑖(𝑡) −𝛽𝐼(𝑡)𝑆𝑣𝑖(𝑡)
# 𝑑𝐸𝑖(𝑡)/𝑑𝑡 = 𝛽𝐼(𝑡)𝑆𝑖(𝑡) + 𝛽𝐼(𝑡)𝑆𝑣𝑖(𝑡)− 𝜎𝐸𝑖(𝑡)
# 𝑑𝐼𝑖(𝑡)/𝑑𝑡 = 𝜎𝐸𝑖 (𝑡) − 𝛾𝐼𝑖(𝑡)
# 𝑑𝑅𝑖(𝑡)/𝑑𝑡 = 𝛾𝐼𝑖(𝑡)+ (1 − 𝑝) ∗ 𝜃(𝑡)𝑆𝑖(𝑡)
# where 𝐼(𝑡) = ∑𝐼𝑖(𝑡) and the force of infection, 𝜆(𝑡) = 𝛽𝐼(𝑡)

# Original OAW model there were 5 groups based on eligibility for the vaccine:
# < 6 months, 6-11 months, 1-11 years, >= 12 years (not pregnant) and
# >= 12 years (pregnant). Those < 6 months and >=12 years (pregnant) did not
# receive MMR vaccination. Homogenous mixing was assumed given the shelter
# conditions.

# Measles Model: SE2IRSv with groups
base_states <- c("S", "E", "I", "R", "Sv")
measlesmodel <- define_states(base_states) |>
  add_infection("I", "S", "E", "beta") |>
  add_transition("I", "R", "tau") |>
  add_transition("E", "I", "taue", chainlength = 2) |>
  add_transition("S", c("Sv", "R"), "1/(theta)",
    forkprobability = c("p", "1-p")
  ) |>
  add_infection("I", "Sv", "E", "beta") |>
  add_group(as.character(1:5))

measlescompiled <- compilemodel(measlesmodel)

measles_rates <- measlescompiled$modeloutstructions$processrates
measles_peter <- measlescompiled$modeloutstructions$petermatrix
measles_states <- measlescompiled$modeloutstructions$updatedstates
measles_rates
measles_peter
measles_states

# For reference this is the model we would have to write down manually!!
# We can generate output with minimal input
# Compare to seirMeasles/R/rates.R and seirMeasles/R/transitions.R
modeleqn <- modeloutput2odeinput(measlescompiled)
length(modeleqn) # 30 differential equations
modeleqn

# Simulation and Plotting
dyn <- wrap_adaptivetau(
  c("S" = 999, "I" = 1, "R" = 0, "E" = 0, "Sv" = 0),
  measlescompiled,
  rate_func = NULL, # defaults to compModels::generalized_rates()
  c(beta = 2, tau = 1, taue = .5, theta = .01, p = 0.4),
  25,
  1,
  "adaptivetau"
)
plot_stoch_model(dyn,
  compartments = measles_states,
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(30)
)



# Simulation with importation on first time step & intervention day 15
pre_vacc <- wrap_adaptivetau(
  c("S" = 999, "I" = 1, "R" = 0, "E" = 0, "Sv" = 0),
  measlescompiled,
  rate_func = NULL,
  c(beta = 2, tau = 1, taue = .5, theta = .01, p = 0.4),
  15,
  1,
  "adaptivetau"
)
last_row <- tail(pre_vacc[[1]], n = 1) |>
  select(-1) |>
  unlist()
last_row

# Vaccinating the three eligible susceptible classes (takes them to R, assumes
# no further primary vaccine failure)
current_states <- modify_states(last_row,
  adjust_names = c(
    "S_group2", "S_group3", "S_group4",
    "R_group2", "R_group3", "R_group4"
  ),
  adjust_value = c(-200, -200, -200, 200, 200, 200)
)

post_vacc <- wrap_adaptivetau(
  current_states,
  measlescompiled,
  rate_func = NULL, # this defaults to compModels::generalized_rates()
  c(beta = 2, tau = 1, taue = .5, theta = .01, p = 0.4),
  25,
  1,
  "adaptivetau",
  "current"
)

post_vacc_mod <- lapply(post_vacc, function(df) {
  df$time <- df$time + 15
  df <- df[-1, ]
  return(df)
})
dyn_int <- Map(rbind, pre_vacc, post_vacc_mod)

plot_stoch_model(dyn_int,
  compartments = measles_states,
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(30)
)
