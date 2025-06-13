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

# Replicating Chicago Measles Response 2024: same populations & parameters
# Overview:

# Simulation start date 2024-02-01, sim length 365 days, 10,000 sims/scenario
# Parameters
# sigma (latent period) = 8 days (mean)
# infectious period = 5 days (mean),
# r_0 = 25  (recall: beta/gamma = r0; also: beta*gamma = r0 = 25)

# Introduce 1 infected into age group 3 on 2024-02-22 (22 timesteps)
# https://github.com/CDCgov/measles-model-chicago-2024/blob/c857905313192741ebff2c61d9a46ca57b808a42/seirMeasles/R/simulate.R#L26 #nolint
# Interventions:
# 1. Vaccination
#    277, 276, & 329 doses on 2024-03-08, 2024-03-09, 2024-03-10
#    (15, 16, 17 days later)
#    Vaccine efficacy: theta_infant (group 2) and theta_adult (groups 3 & 4)
#    Existing documented_vax_coverage: 0.49
# 2. Active case finding: started 2024-03-08 with efficacy 0.25

# Not included here:
# Tracking compartments: cases, pre-reported cases, reported cases
# Case ascertainment_delay: 2.5 days (mean)
# https://github.com/CDCgov/measles-model-chicago-2024/blob/c857905313192741ebff2c61d9a46ca57b808a42/seirMeasles/R/rates.R#L65C7-L65C55 #nolint
# https://github.com/CDCgov/measles-model-chicago-2024/blob/c857905313192741ebff2c61d9a46ca57b808a42/seirMeasles/R/simulate.R#L345 #nolint


base_states <- c("S", "E", "I", "R", "Sv")
measlesmodel <- define_states(base_states) |>
  add_infection("I", "S", "E", "beta") |>
  add_transition("I", "R", meantime = "gamma*(1-int_eff)") |>
  add_transition("E", "I", meantime = "sigma", chainlength = 2) |>
  add_infection("I", "Sv", "E", "beta") |>
  add_group(groupnames = as.character(1:5))

measlescompiled <- compilemodel(measlesmodel)
measles_rates <- measlescompiled$modeloutstructions$processrates
measles_peter <- measlescompiled$modeloutstructions$petermatrix
measles_states <- measlescompiled$modeloutstructions$updatedstates


# Population sizes
# https://github.com/CDCgov/measles-model-chicago-2024/blob/c857905313192741ebff2c61d9a46ca57b808a42/seirMeasles/inst/params.yaml#L6 #nolint
pop_sizes <- setNames(
  c(8, 13, 704, 1139, 13), # 1877 total people
  measles_states[grep("S_", measles_states)]
)
# Susceptible population sizes
# https://github.com/CDCgov/measles-model-chicago-2024/blob/c857905313192741ebff2c61d9a46ca57b808a42/seirMeasles/inst/params.yaml#L13 #nolint
# https://github.com/CDCgov/measles-model-chicago-2024/blob/c857905313192741ebff2c61d9a46ca57b808a42/seirMeasles/R/simulate.R#L171 #nolint
pop_sizes_s <- floor(pop_sizes * c(0.3333, 0.667, 0.172, 0.1, 0.1))
pop_sizes_r <- setNames(
  (pop_sizes - pop_sizes_s),
  measles_states[grep("R_", measles_states)]
)
initial_pop <- c(pop_sizes_s, pop_sizes_r)
tblpopsize <- define_popsize(measlescompiled, inputpops = initial_pop)
vecpopsize <- setNames(tblpopsize$popsize, tblpopsize$updatedstate)

# First simulation: start to infectious import
# int_eff = 0 (case finding intervention has not yet started)
pre_import <- wrap_adaptivetau(
  vecpopsize,
  measlescompiled,
  rate_func = NULL,
  c(beta = 5, gamma = 5, int_eff = 0, sigma = 8),
  22,
  1,
  "adaptivetau"
)

lr_preimport <- tail(pre_import[[1]], n = 1) |>
  dplyr::select(-1) |> # remove time column
  unlist()
lr_preimport

# Importing 1 infection (note: this makes total population size + 1)
inf_import <- 1 # 1 infected in group 3
tblpopsize <- define_popsize(measlescompiled, inputpops = lr_preimport) |>
  addpopsize_byfeature(inf_import, basestates = "I", groupnames = "3")
vecpopsize <- setNames(tblpopsize$popsize, tblpopsize$updatedstate)
# vecpopsize[grep("I_", names(vecpopsize))] # checking infection #nolint

# Second simulation: infectious import to first vaccination day
# int_eff = 0 (case finding intervention has not yet started)
import_prevacc <- wrap_adaptivetau(
  vecpopsize,
  measlescompiled,
  rate_func = NULL,
  c(beta = 5, gamma = 5, int_eff = 0, sigma = 8),
  15 + 7, # adding on immunity onset delay
  1,
  "adaptivetau"
)

# nolint start
# Visual and other checks
# plot_stoch_model(import_prevacc,
#   compartments = measles_states[grep("I_", measles_states)],
#   colors =
#     colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(5)
# )

# Finding sims where last row has non-zero E values
# check_non_zero_E <- function(df) {
#   last_row <- tail(df, 1)
#   E_columns <- grep("^E_", names(last_row), value = TRUE)
#   any(last_row[E_columns] != 0)
# }
# filtered_df_list <- lapply(import_prevacc, function(df) {
#   if (check_non_zero_E(df)) {
#     return(df)
#   } else {
#     return(NULL)
#   }
# })
# filtered_df_list <- filtered_df_list[!sapply(filtered_df_list, is.null)]
# tail(filtered_df_list[[1]], 1)
# nolint end

lr_postimport_prevacc1 <- tail(import_prevacc[[1]], n = 1) |>
  dplyr::select(-1) |>
  unlist()
lr_postimport_prevacc1


# Vaccination Day 1
# Case Ascertainment intervention starts (int_eff = 0.25)
# https://github.com/CDCgov/measles-model-chicago-2024/blob/c857905313192741ebff2c61d9a46ca57b808a42/seirMeasles/R/simulate.R#L82 #nolint
doses_day1 <- 277
theta_infant <- 0.84
theta_adult <- 0.925
prop_s <- sum(pop_sizes_s[2:4]) / sum(pop_sizes[2:4])
vacc_coverage <- 0.49
prop_vax_to_se <- prop_s / (1.0 - vacc_coverage)

# Get S and E populations and distribute vaccines to specific compartments
se <- lr_postimport_prevacc1[grep("^(S_|E_)", names(lr_postimport_prevacc1))]
se_doses_total <- doses_day1 * prop_vax_to_se
se_doses <- round(pmin(se_doses_total * se / sum(se), se))


tblpopsize <- define_popsize(measlescompiled,
  inputpops = lr_postimport_prevacc1
) |>
  # assume all E vaccinations are failures so not coding a move to R
  movepopsize_byname(round(se_doses["S_dummygroup1X2"] * theta_infant),
    namebefore = "S_dummygroup1X2",
    nameafter = "R_dummygroup1X2"
  ) |>
  movepopsize_byname(round(se_doses["S_dummygroup1X3"] * theta_adult),
    namebefore = "S_dummygroup1X3",
    nameafter = "R_dummygroup1X3"
  ) |>
  movepopsize_byname(round(se_doses["S_dummygroup1X4"] * theta_adult),
    namebefore = "S_dummygroup1X4",
    nameafter = "R_dummygroup1X4"
  ) |>
  movepopsize_byname(round(se_doses["S_dummygroup1X2"] * (1 - theta_infant)),
    namebefore = "S_dummygroup1X2",
    nameafter = "Sv_dummygroup1X2"
  ) |>
  movepopsize_byname(round(se_doses["S_dummygroup1X3"] * (1 - theta_adult)),
    namebefore = "S_dummygroup1X3",
    ameafter = "Sv_dummygroup1X3"
  ) |>
  movepopsize_byname(round(se_doses["S_dummygroup1X4"] * (1 - theta_adult)),
    namebefore = "S_dummygroup1X4",
    nameafter = "Sv_dummygroup1X4"
  )

vecpopsize <- setNames(tblpopsize$popsize, tblpopsize$updatedstate)

# Third simulation: first vaccination day to second vaccination day
postimport_vacc1 <- wrap_adaptivetau(
  vecpopsize,
  measlescompiled,
  rate_func = NULL,
  c(beta = 5, gamma = 5, int_eff = 0.25, sigma = 8),
  1, # immunity onset delay already added in previous sim
  1,
  "adaptivetau"
)

lr_postimport_postvacc1 <- tail(postimport_vacc1[[1]], n = 1) |>
  dplyr::select(-1) |>
  unlist()
lr_postimport_postvacc1

# Vaccination Day 2
# Case Ascertainment intervention continues (int_eff = 0.25)
doses_day2 <- 276
se <- lr_postimport_postvacc1[grep(
  "^(S_|E_)",
  names(lr_postimport_postvacc1)
)]
se_doses_total <- doses_day2 * prop_vax_to_se
se_doses <- round(pmin(se_doses_total * se / sum(se), se))

tblpopsize <- define_popsize(measlescompiled,
  inputpops = lr_postimport_postvacc1
) |>
  # assume all E vaccinations are failures so not coding a move to R
  movepopsize_byname(round(se_doses["S_dummygroup1X2"] * theta_infant),
    namebefore = "S_dummygroup1X2",
    nameafter = "R_dummygroup1X2"
  ) |>
  movepopsize_byname(round(se_doses["S_dummygroup1X3"] * theta_adult),
    namebefore = "S_dummygroup1X3",
    nameafter = "R_dummygroup1X3"
  ) |>
  movepopsize_byname(round(se_doses["S_dummygroup1X4"] * theta_adult),
    namebefore = "S_dummygroup1X4",
    nameafter = "R_dummygroup1X4"
  ) |>
  movepopsize_byname(round(se_doses["S_dummygroup1X2"] * (1 - theta_infant)),
    namebefore = "S_dummygroup1X2",
    nameafter = "Sv_dummygroup1X2"
  ) |>
  movepopsize_byname(round(se_doses["S_dummygroup1X3"] * (1 - theta_adult)),
    namebefore = "S_dummygroup1X3",
    ameafter = "Sv_dummygroup1X3"
  ) |>
  movepopsize_byname(round(se_doses["S_dummygroup1X4"] * (1 - theta_adult)),
    namebefore = "S_dummygroup1X4",
    nameafter = "Sv_dummygroup1X4"
  )
# lr_postimport_postvacc1[grep("^(S_|Sv_|R)", names(lr_postimport_postvacc1))] #nolint
vecpopsize <- setNames(tblpopsize$popsize, tblpopsize$updatedstate)

# Fourth simulation: second vaccination day to third vaccination day
postimport_vacc2 <- wrap_adaptivetau(
  vecpopsize,
  measlescompiled,
  rate_func = NULL,
  c(beta = 5, gamma = 5, int_eff = 0.25, sigma = 8),
  1, # immunity onset delay already added in previous sim
  1,
  "adaptivetau"
)

lr_postimport_postvacc2 <- tail(postimport_vacc2[[1]], n = 1) |>
  dplyr::select(-1) |>
  unlist()
lr_postimport_postvacc2

# Vaccination Day 3
# Case Ascertainment intervention continues (int_eff = 0.25)
doses_day3 <- 329
se <- lr_postimport_postvacc2[grep(
  "^(S_|E_)",
  names(lr_postimport_postvacc2)
)]
se_doses_total <- doses_day3 * prop_vax_to_se
se_doses <- round(pmin(se_doses_total * se / sum(se), se))

tblpopsize <- define_popsize(measlescompiled,
  inputpops = lr_postimport_postvacc2
) |>
  # assume all E vaccinations are failures so not coding a move to R
  movepopsize_byname(round(se_doses["S_dummygroup1X2"] * theta_infant),
    namebefore = "S_dummygroup1X2",
    nameafter = "R_dummygroup1X2"
  ) |>
  movepopsize_byname(round(se_doses["S_dummygroup1X3"] * theta_adult),
    namebefore = "S_dummygroup1X3",
    nameafter = "R_dummygroup1X3"
  ) |>
  movepopsize_byname(round(se_doses["S_dummygroup1X4"] * theta_adult),
    namebefore = "S_dummygroup1X4",
    nameafter = "R_dummygroup1X4"
  ) |>
  movepopsize_byname(round(se_doses["S_dummygroup1X2"] * (1 - theta_infant)),
    namebefore = "S_dummygroup1X2",
    nameafter = "Sv_dummygroup1X2"
  ) |>
  movepopsize_byname(round(se_doses["S_dummygroup1X3"] * (1 - theta_adult)),
    namebefore = "S_dummygroup1X3",
    ameafter = "Sv_dummygroup1X3"
  ) |>
  movepopsize_byname(round(se_doses["S_dummygroup1X4"] * (1 - theta_adult)),
    namebefore = "S_dummygroup1X4",
    nameafter = "Sv_dummygroup1X4"
  )

vecpopsize <- setNames(tblpopsize$popsize, tblpopsize$updatedstate)

# Fifth simulation: third vaccination day to end of simulation (365 days)
postimport_vacc3 <- wrap_adaptivetau(
  vecpopsize,
  measlescompiled,
  rate_func = NULL,
  c(beta = 5, gamma = 5, int_eff = 0.25, sigma = 8),
  320, # remaining days, after immunity onset delay taken into account
  1,
  "adaptivetau"
)

# Stitching together the simulations
# Edit start times of the simulated output files ahead of appending (rbind)
import_prevacc_mod <- lapply(import_prevacc, function(df) {
  df$time <- df$time + 22
  df <- df[-1, ]
  return(df)
})

postimport_vacc1_mod <- lapply(postimport_vacc1, function(df) {
  df$time <- df$time + (22 + 15 + 7)
  df <- df[-1, ]
  return(df)
})

postimport_vacc2_mod <- lapply(postimport_vacc2, function(df) {
  df$time <- df$time + (22 + 15 + 7 + 1)
  df <- df[-1, ]
  return(df)
})

postimport_vacc3_mod <- lapply(postimport_vacc3, function(df) {
  df$time <- df$time + (22 + 15 + 7 + 1 + 1)
  df <- df[-1, ]
  return(df)
})

all_dyn <- Map(
  rbind, pre_import, import_prevacc_mod, postimport_vacc1_mod,
  postimport_vacc2_mod, postimport_vacc3_mod
)
all_dyn


# Plot compartments
plot_stoch_model(all_dyn,
  compartments = measles_states[grep("^(S_)", measles_states)],
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(5)
)

plot_stoch_model(all_dyn,
  compartments = measles_states[grep("^(Sv_)", measles_states)],
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(5)
)

plot_stoch_model(all_dyn,
  compartments = measles_states[grep("^(E_)", measles_states)],
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(10)
)

plot_stoch_model(all_dyn,
  compartments = measles_states[grep("^(I_)", measles_states)],
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(5)
)

plot_stoch_model(all_dyn,
  compartments = measles_states[grep("^(R_)", measles_states)],
  colors =
    colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(5)
)
