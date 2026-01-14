library(devtools)
load_all("compModels")

base_states <- c("S", "E", "I", "R")
infector_state <- "I"
infectee_state <- "S"

seirmodel <- define_states(base_states) |>
  add_infection(infector_state, infectee_state, "E", "beta") |>
  add_transition("E", "I", "alpha") |>
  add_transition("I", "R", "gamma")

seircompiled <- compilemodel(seirmodel)

initial_states_vec <- c("S" = 999, "I" = 1, "R" = 0, "E" = 0)
param_vector <- c("beta" = 0.3, "alpha" = 3, "gamma" = 5)
nsims <- 100
ntimesteps <- 100

odin_mod_results <- wrap_odin(
  init_vals = initial_states_vec,
  compiledmodel = seircompiled,
  param_vector = param_vector,
  n_timesteps = ntimesteps,
  n_sims = nsims
)
ggplot(odin_mod_results, aes(x = step, y = I, group = sim)) +
  geom_line() +
  labs(title = "SEIR odin simulation", x = "Time", y = "Infected") +
  theme_minimal()
