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
# SIRmodelbuilder functions are now in compModels/R
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

# SIR Gillespie
dyn <- wrap_gillespie(
  c("S" = 999, "I" = 1, "R" = 0),
  sircompiled,
  c(beta = 2, tau = 1),
  25,
  1
)
plot_stoch_model(dyn, compartments = c("S", "I", "R"))

dyn2 <- wrap_gillespie(
  c("S" = 999, "I" = 1, "R" = 0),
  sircompiled,
  c(beta = 2, tau = 1),
  25,
  10
)
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

dyn3 <- wrap_gillespie(
  c("S" = 999, "I" = 1, "R" = 0, "E" = 0),
  seircompiled,
  c(beta = 2, tau = 1, taue = .5),
  25,
  10
)
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

dyn4 <- wrap_gillespie(
  c("S" = 999, "I" = 1, "R" = 0, "E" = 0),
  sei2rcompiled,
  c(beta = 2, tau = 1, taue = .5),
  25,
  10
)
plot_stoch_model(dyn4, compartments = c("S", "E", "I_chain1", "I_chain2", "R"))

# SE-Pro-I2R with prodromal period
base_states <- c("S", "E", "Pro", "I", "R")
sepi2r <- define_states(base_states) |>
  add_infection("I", "S", "E", "beta") |>
  add_infection("Pro", "S", "E", "beta_pro") |>
  add_transition("I", "R", "tau", chainlength = 2) |>
  add_transition("E", "Pro", "taue") |>
  add_transition("Pro", "I", "taup")
sepi2rcompiled <- compilemodel(sepi2r)
sepi2r_rates <- sepi2rcompiled$modeloutstructions$processrates
sepi2r_peter <- sepi2rcompiled$modeloutstructions$petermatrix
sepi2r_states <- sepi2rcompiled$modeloutstructions$updatedstates

dyn4p <- wrap_gillespie(
  c("S" = 999, "I" = 1, "R" = 0, "E" = 0),
  sepi2rcompiled,
  c(beta = 2, beta_pro = 2.5, tau = 1, taue = .5, taup = 0.6),
  25,
  10
)
plot_stoch_model(dyn4p, compartments = c(
  "S", "E", "Pro",
  "I_chain1", "I_chain2", "R"
))


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

dyn5 <- wrap_gillespie(
  c("S" = 999, "I" = 1, "R" = 0, "E" = 0, "H" = 0),
  sei2rhcompiled,
  c(beta = 2, tau = 1, taue = .5, probH = .2),
  25,
  10
)
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

dyn6 <- wrap_gillespie(
  c("S" = 999, "I" = 1, "R" = 0, "E" = 0),
  sei2r_treatcompiled,
  c(beta = 2, tau = 1, taue = .5, taut = 2),
  25,
  10
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

dyn7 <- wrap_gillespie(
  c("S" = 999, "I" = 1, "R" = 0),
  outmeta,
  c(beta = 2, tau = 1, mu = .1),
  25,
  10
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

dyn8 <- wrap_gillespie(
  c("S" = 999, "I" = 1, "R" = 0),
  sir1groupcomp,
  c(beta = 1.5, tau = 1),
  25,
  10
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

dyn9 <- wrap_gillespie(
  c("S" = 999, "I" = 1, "R" = 0),
  sir2group,
  c(beta = 2, tau = 1),
  25,
  10
)

plot_stoch_model(dyn9,
  compartments = c(
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
  ),
  colors = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(12)
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

dyn10 <- wrap_gillespie(
  c("S" = 999, "I" = 1, "R" = 0),
  outlistfull,
  c(beta = 2, tau = 1, mu = 0.1), # guess mu for travel,
  25,
  10
)

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
  ),
  colors = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(40)
)
