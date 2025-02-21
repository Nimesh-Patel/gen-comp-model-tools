# Examples of GillespieSSA:ssa() with SIRModelBuilder Approach

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
load_all("SIRmodelbuilder")


# SIR instructions
base_states <- c("S", "I", "R")
sir <- define_states(base_states) |>
  add_infection("I", "S", "I", "beta") |>
  add_transition("I", "R", "tau")
sircompiled <- compilemodel(sir)
sir_rates <- sircompiled$modeloutstructions$processrates
sir_peter <- sircompiled$modeloutstructions$petermatrix
sir_states <- sircompiled$modeloutstructions$updatedstates

# SIR Gillespie
x0 <- define_initialstate(
  sircompiled,
  c("S" = 999, "I" = 1, "R" = 0)
) |>
  output_initialstate()

parameters <- c(beta = 2, tau = 1, taue = .5)

sir_peter2 <- matrix(as.numeric(sir_peter), ncol = 2, nrow = 3)
dyn <- GillespieSSA::ssa(
  x0 = x0,
  a = sir_rates,
  nu = sir_peter2,
  parms = parameters,
  tf = 25,
  method = GillespieSSA::ssa.d(),
  simName = "General Stochastic Model"
)$data

dynplt <- tibble::as_tibble(dyn) |>
  tidyr::pivot_longer(!t) |>
  dplyr::mutate(density = as.numeric(value), t = as.numeric(t))

ggplot2::ggplot(dynplt, aes(x = t, y = density, color = name)) +
  ggplot2::geom_line(linewidth = 2) +
  ggplot2::theme_classic(base_size = 16)


# Create Gillespie runner for multiple sims
run_gillespie <- function(init_vals, propensity_fns,
                          change_matrix, parameters,
                          n_timesteps, n_sims) {
  nu <- change_matrix
  t <- n_timesteps

  sims <- lapply(1:n_sims, function(i) {
    current_parameters <- parameters

    sim_data <- GillespieSSA::ssa(
      x0 = init_vals,
      a = propensity_fns,
      nu = nu,
      parms = current_parameters,
      tf = t,
      method = GillespieSSA::ssa.d(),
      simName = "General Stochastic Model With Intervention"
    )$data
  })

  lapply(sims, function(x) data.frame(x))
}


dyn2 <- run_gillespie(x0, sir_rates, sir_peter2, parameters, 25, 10)

# Brining over plot functionality from compModels
plot_stoch_model <- function(output,
                             compartments = c("I"),
                             time_var = NULL,
                             show_intervention_period = FALSE,
                             intervention_period = NULL,
                             colors = RColorBrewer::brewer.pal(8, "Dark2")) {
  possible_time_vars <- c("t", "time")

  if (is.null(time_var)) {
    found_time_var <- FALSE
    for (possible_time_var in possible_time_vars) {
      if (possible_time_var %in% colnames(output[[1]])) {
        time_var <- possible_time_var
        found_time_var <- TRUE
        print(paste0("Found time variable:", possible_time_var))
        break
      }
    }
    if (!found_time_var) {
      stop("Time variable not found in output data.")
    }
  }

  combined_sims_df <- dplyr::bind_rows(lapply(seq_along(output), function(i) {
    cbind(simulation = i, output[[i]])
  }), .id = "simulation_id")

  unique_compartments <- unique(compartments)
  num_colors <- length(unique_compartments)
  colors <- colors[1:num_colors]
  print(paste0(
    "There are ", num_colors, " unique compartments called ",
    paste(unique_compartments, collapse = ", "),
    ", which are given these colors:",
    paste(colors, collapse = ", ")
  ))

  combined_sims_df <- combined_sims_df |>
    tidyr::pivot_longer(
      cols = tidyselect::all_of(compartments),
      names_to = "compartment",
      values_to = "value"
    )

  p <- ggplot(
    combined_sims_df,
    aes(
      x = !!sym(time_var),
      y = .data$value,
      color = .data$compartment,
      group = interaction(
        .data$simulation_id,
        .data$compartment
      )
    )
  ) +
    geom_line(alpha = 0.5) +
    theme_classic()

  if (show_intervention_period == TRUE && !is.null(intervention_period)) {
    p <- p + geom_vline(
      xintercept = intervention_period[1],
      linetype = "longdash"
    ) + geom_vline(
      xintercept = intervention_period[2],
      linetype = "longdash"
    )
  }

  # Add labels and colors
  p <- p + labs(
    title = "Stochastic SIR Model Simulations",
    x = "Time",
    y = "Number of Individuals"
  ) + scale_color_manual(values = colors, name = "Compartment")

  return(p)
}

plot_stoch_model(dyn2, compartments = c("S", "I", "R"))


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

x0 <- define_initialstate(
  seircompiled,
  c("S" = 999, "I" = 1, "R" = 0, "E" = 0)
) |>
  output_initialstate()
parameters <- c(beta = 2, tau = 1, taue = .5)

seir_peter2 <- matrix(as.numeric(seir_peter), ncol = 3, nrow = 4)
dyn3 <- run_gillespie(x0, seir_rates, seir_peter2, parameters, 25, 10)
plot_stoch_model(dyn3, compartments = c("S", "E", "I", "R"))



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

x0 <- define_initialstate(
  sei2rcompiled,
  c("S" = 999, "I" = 1, "R" = 0, "E" = 0)
) |>
  output_initialstate()

parameters <- c(beta = 2, tau = 1, taue = .5)

sei2r_peter2 <- matrix(as.numeric(sei2r_peter), ncol = 5, nrow = 5)
dyn4 <- run_gillespie(x0, sei2r_rates, sei2r_peter2, parameters, 25, 10)
plot_stoch_model(dyn4, compartments = c("S", "E", "I_chain1", "I_chain2", "R"))


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

x0 <- define_initialstate(
  sei2rhcompiled,
  c("S" = 999, "I" = 1, "R" = 0, "E" = 0, "H" = 0)
) |>
  output_initialstate()
parameters <- c(beta = 2, tau = 1, taue = .5, probH = .2)

sei2rh_peter2 <- matrix(as.numeric(sei2rh_peter), ncol = 6, nrow = 6)
dyn5 <- run_gillespie(x0, sei2rh_rates, sei2rh_peter2, parameters, 25, 10)
plot_stoch_model(dyn5, compartments = c(
  "S", "E", "I_chain1", "I_chain2",
  "R", "H"
))


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

x0 <- define_initialstate(
  sei2r_treatcompiled,
  c("S" = 999, "I" = 1, "R" = 0, "E" = 0)
) |>
  output_initialstate()
parameters <- c(beta = 2, tau = 1, taue = .5, taut = 2)

sei2r_treat_peter2 <- matrix(as.numeric(sei2r_treat_peter), ncol = 4, nrow = 4)
dyn6 <- run_gillespie(
  x0, sei2r_treat_rates, sei2r_treat_peter2,
  parameters, 25, 10
) # try with end time = 50 to see more
plot_stoch_model(dyn6, compartments = c("S", "E", "I", "R"))

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

x0meta <- define_initialstate(outmeta, c("S" = 999, "I" = 1, "R" = 0)) |>
  output_initialstate()
parameters <- c(beta = 2, tau = 1, mu = .1)


outmeta_peter2 <- matrix(as.numeric(outmeta_peter), ncol = 10, nrow = 6)
dyn7 <- run_gillespie(
  x0meta, outmeta_rates, outmeta_peter2,
  parameters, 25, 10
)
plot_stoch_model(dyn7, compartments = c(
  "S_metapopulationR0",
  "S_metapopulationR0times2",
  "I_metapopulationR0",
  "I_metapopulationR0times2",
  "R_metapopulationR0",
  "R_metapopulationR0times2"
))

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

x0group1 <- define_initialstate(
  sir1groupcomp,
  c("S" = 999, "I" = 1, "R" = 0)
) |>
  output_initialstate()
parameters <- c(beta = 1.5, tau = 1)

sir1groupcomp_peter2 <- matrix(as.numeric(sir1groupcomp_peter),
  ncol = 6, nrow = 6
)
dyn8 <- run_gillespie(
  x0group1, sir1groupcomp_rates, sir1groupcomp_peter2,
  parameters, 25, 10
)
plot_stoch_model(dyn8, compartments = c(
  "S_groupSocial",
  "I_groupSocial",
  "R_groupSocial",
  "S_groupAntisocial",
  "I_groupAntisocial",
  "R_groupAntisocial"
))

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

x0group2 <- define_initialstate(sir2group, c("S" = 999, "I" = 1, "R" = 0)) |>
  output_initialstate()
parameters <- c(beta = 2, tau = 1)

sir2group_peter2 <- matrix(as.numeric(sir2group_peter),
  ncol = 20, nrow = 12
)
dyn9 <- run_gillespie(
  x0group2, sir2group_rates, sir2group_peter2,
  parameters, 25, 10
)
plot_stoch_model(dyn9, compartments = c(
  "S_AgeYoung_HospitalPatient",
  "I_AgeYoung_HospitalPatient",
  "R_AgeYoung_HospitalPatient",
  "S_AgeOld_HospitalPatient",
  "I_AgeOld_HospitalPatient",
  "R_AgeOld_HospitalPatient",
  "S_AgeYoung_HospitalHCW",
  "I_AgeYoung_HospitalHCW",
  "R_AgeYoung_HospitalHCW",
  "S_AgeOld_HospitalHCW",
  "I_AgeOld_HospitalHCW",
  "R_AgeOld_HospitalHCW"
))


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

x0full <- define_initialstate(outlistfull, c("S" = 999, "I" = 1, "R" = 0)) |>
  output_initialstate()
parameters <- c(beta = 2, tau = 1, mu = 0.1) # guess mu for travel

outlistfull_peter2 <- matrix(as.numeric(outlistfull_peter),
  ncol = 160, nrow = 40
)
dyn10 <- run_gillespie(
  x0full, outlistfull_rates, outlistfull_peter2,
  parameters, 50, 1
)
# This plot will need some owrk - ran out of color palette!
plot_stoch_model(dyn10,
  compartments = c(
    "S_AgeYoung_HospitalPatient_metapopulationUK",
    "S_AgeYoung_HospitalPatient_metapopulationUSA",
    "I_AgeYoung_HospitalPatient_metapopulationUK_chain1",
    "I_AgeYoung_HospitalPatient_metapopulationUK_chain2",
    "I_AgeYoung_HospitalPatient_metapopulationUK_chain3",
    "I_AgeYoung_HospitalPatient_metapopulationUSA_chain1",
    "I_AgeYoung_HospitalPatient_metapopulationUSA_chain2",
    "I_AgeYoung_HospitalPatient_metapopulationUSA_chain3",
    "R_AgeYoung_HospitalPatient_metapopulationUK",
    "R_AgeYoung_HospitalPatient_metapopulationUSA",
    "S_AgeOld_HospitalPatient_metapopulationUK",
    "S_AgeOld_HospitalPatient_metapopulationUSA",
    "I_AgeOld_HospitalPatient_metapopulationUK_chain1",
    "I_AgeOld_HospitalPatient_metapopulationUK_chain2",
    "I_AgeOld_HospitalPatient_metapopulationUK_chain3",
    "I_AgeOld_HospitalPatient_metapopulationUSA_chain1",
    "I_AgeOld_HospitalPatient_metapopulationUSA_chain2",
    "I_AgeOld_HospitalPatient_metapopulationUSA_chain3",
    "R_AgeOld_HospitalPatient_metapopulationUK",
    "R_AgeOld_HospitalPatient_metapopulationUSA",
    "S_AgeYoung_HospitalHCW_metapopulationUK",
    "S_AgeYoung_HospitalHCW_metapopulationUSA",
    "I_AgeYoung_HospitalHCW_metapopulationUK_chain1",
    "I_AgeYoung_HospitalHCW_metapopulationUK_chain2",
    "I_AgeYoung_HospitalHCW_metapopulationUK_chain3",
    "I_AgeYoung_HospitalHCW_metapopulationUSA_chain1",
    "I_AgeYoung_HospitalHCW_metapopulationUSA_chain2",
    "I_AgeYoung_HospitalHCW_metapopulationUSA_chain3",
    "R_AgeYoung_HospitalHCW_metapopulationUK",
    "R_AgeYoung_HospitalHCW_metapopulationUSA",
    "S_AgeOld_HospitalHCW_metapopulationUK",
    "S_AgeOld_HospitalHCW_metapopulationUSA",
    "I_AgeOld_HospitalHCW_metapopulationUK_chain1",
    "I_AgeOld_HospitalHCW_metapopulationUK_chain2",
    "I_AgeOld_HospitalHCW_metapopulationUK_chain3",
    "I_AgeOld_HospitalHCW_metapopulationUSA_chain1",
    "I_AgeOld_HospitalHCW_metapopulationUSA_chain2",
    "I_AgeOld_HospitalHCW_metapopulationUSA_chain3",
    "R_AgeOld_HospitalHCW_metapopulationUK",
    "R_AgeOld_HospitalHCW_metapopulationUSA"
  )
)
