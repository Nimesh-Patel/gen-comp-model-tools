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
load_all("compModels")

# Loading a modified wrap_adaptivetau function to enable run of any model below
# in which base states get modified (eg chains, groups, metapopulations)
# This may get incorporated into compModels in the near future
wrap_gillespie <- function(init_vals, compiledmodel, parameters, n_timesteps,
                           n_sims) {
  x0 <- init_vals
  model_rates <- compiledmodel$modeloutstructions$processrates # propensity
  model_peter <- compiledmodel$modeloutstructions$petermatrix # change matrix
  parameters <- parameters

  # Ensure x0 matches model dimensions (nu rows)
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

  nu <- matrix(as.numeric(model_peter),
    ncol = dim(model_peter)[2],
    nrow = dim(model_peter)[1]
  )
  # rownames(nu) <- names(x0) #nolint does not appear needed

  t <- n_timesteps

  sims <- lapply(1:n_sims, function(i) {
    current_parameters <- parameters

    sim_data <- GillespieSSA::ssa(
      x0 = x0,
      a = model_rates,
      nu = nu,
      parms = current_parameters,
      tf = t,
      method = GillespieSSA::ssa.d()
    )$data |> data.frame()
  })
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
  c("S" = 999, "E" = 1, "I" = 1, "R" = 0),
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
  c("S" = 999, "E" = 1, "I" = 1, "R" = 0),
  sei2rcompiled,
  c(beta = 2, tau = 1, taue = .5),
  25,
  10
)
plot_stoch_model(dyn4, compartments = c(
  "S", "E",
  "I_chainid1X1", "I_chainid1X2",
  "R"
))

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
  c("S" = 999, "E" = 0, "I" = 1, "R" = 0),
  sepi2rcompiled,
  c(beta = 2, beta_pro = 2.5, tau = 1, taue = .5, taup = 0.6),
  25,
  10
)
plot_stoch_model(dyn4p, compartments = c(
  "S", "E", "Pro",
  "I_chainid1X1", "I_chainid1X2", "R"
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
  c("S" = 999, "E" = 0, "I" = 1, "R" = 0, "H" = 0),
  sei2rhcompiled,
  c(beta = 2, tau = 1, taue = .5, probH = .2),
  25,
  10
)
plot_stoch_model(dyn5, compartments = c(
  "S", "E", "I_chainid1X1", "I_chainid1X2",
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
  c("S" = 999, "E" = 0, "I" = 1, "R" = 0),
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
    metapopulation = c("R0", "R0times2"),
    scaletransitions = c(1, 2)
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
  add_group(
    groupname = c("Social", "Antisocial"),
    scaleprocessbyname = list(infection = c(2, 1))
  )
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
  "S_dummygroup1Social",
  "I_dummygroup1Social",
  "R_dummygroup1Social",
  "S_dummygroup1Antisocial",
  "I_dummygroup1Antisocial",
  "R_dummygroup1Antisocial"
))

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

dyn9 <- wrap_gillespie(
  c("S" = 999, "I" = 1, "R" = 0),
  sir2group,
  c(beta = 2, tau = 1),
  25,
  10
)

plot_stoch_model(dyn9,
  compartments = c(
    "S_AgeYoung",
    "I_AgeYoung",
    "R_AgeYoung",
    "S_AgeOld",
    "I_AgeOld",
    "R_AgeOld",
    "S_HospitalHCW",
    "I_HospitalHCW",
    "R_HospitalHCW",
    "S_HospitalPatient",
    "I_HospitalPatient",
    "R_HospitalPatient"
  ),
  colors = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(12)
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

dyn10 <- wrap_gillespie(
  c("S" = 999, "I" = 1, "R" = 0),
  outlistfull,
  c(beta = 2, tau = 1, mu = 0.1), # guess mu for travel,
  25,
  10
)

plot_stoch_model(dyn10,
  compartments = c(
    "I_AgeOld_metapopulationUK_chainid1X1",
    "I_AgeOld_metapopulationUK_chainid1X2",
    "I_AgeOld_metapopulationUK_chainid1X3",
    "I_AgeOld_metapopulationUSA_chainid1X1",
    "I_AgeOld_metapopulationUSA_chainid1X2",
    "I_AgeOld_metapopulationUSA_chainid1X3",
    "I_AgeYoung_metapopulationUK_chainid1X1",
    "I_AgeYoung_metapopulationUK_chainid1X2",
    "I_AgeYoung_metapopulationUK_chainid1X3",
    "I_AgeYoung_metapopulationUSA_chainid1X1",
    "I_AgeYoung_metapopulationUSA_chainid1X2",
    "I_AgeYoung_metapopulationUSA_chainid1X3",
    "I_HospitalHCW_metapopulationUK_chainid1X1",
    "I_HospitalHCW_metapopulationUK_chainid1X2",
    "I_HospitalHCW_metapopulationUK_chainid1X3",
    "I_HospitalHCW_metapopulationUSA_chainid1X1",
    "I_HospitalHCW_metapopulationUSA_chainid1X2",
    "I_HospitalHCW_metapopulationUSA_chainid1X3",
    "I_HospitalPatient_metapopulationUK_chainid1X1",
    "I_HospitalPatient_metapopulationUK_chainid1X2",
    "I_HospitalPatient_metapopulationUK_chainid1X3",
    "I_HospitalPatient_metapopulationUSA_chainid1X1",
    "I_HospitalPatient_metapopulationUSA_chainid1X2",
    "I_HospitalPatient_metapopulationUSA_chainid1X3",
    "S_AgeYoung_metapopulationUK",
    "R_AgeYoung_metapopulationUK",
    "S_AgeOld_metapopulationUK",
    "R_AgeOld_metapopulationUK",
    "S_HospitalPatient_metapopulationUK",
    "R_HospitalPatient_metapopulationUK",
    "S_HospitalHCW_metapopulationUK",
    "R_HospitalHCW_metapopulationUK",
    "S_AgeYoung_metapopulationUSA",
    "R_AgeYoung_metapopulationUSA",
    "S_AgeOld_metapopulationUSA",
    "R_AgeOld_metapopulationUSA",
    "S_HospitalPatient_metapopulationUSA",
    "R_HospitalPatient_metapopulationUSA",
    "S_HospitalHCW_metapopulationUSA",
    "R_HospitalHCW_metapopulationUSA"
  ),
  colors = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(40)
)
