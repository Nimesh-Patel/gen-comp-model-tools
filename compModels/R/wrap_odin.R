#' Takes input from compiledmodel and additional model information,
#' creates model in format for odin (needs to be saved as an R file),
#' runs odin model, returns results in dataframe
#'
#' @param init_vals named vector of initial values for each state
#' Example: init_vals <- c("S" = 999, "I" = 1, "R" = 0, "E" = 0)
#' @param compiledmodel compiled model object
#' created by compModels function compilemodel()
#' @param param_vector vector of parameter names and values
#' Example: param_vector <-c("beta"= 2, "alpha" = 1, "gamma" = 1)
#' @param n_timesteps number of timesteps to run the model for
#' @param n_sims number of simulations to run
#' @param t0 initial time for model run, default is 0
#'
#' @return dataframe of results from odin model run
#' @family simulation
#' @export
wrap_odin <- function(init_vals,
                      compiledmodel,
                      param_vector,
                      n_timesteps,
                      n_sims,
                      t0 = 0) {
  model_states <- compiledmodel$modeloutstructions$updatedstates
  model_pop <- compiledmodel$modelinstructions$tblntotal

  odin_init <- odin_initial_states(model_states, init_vals)
  odin_tot <- odin_tot_pop(model_pop)
  core_eqns <- odin_core_eqns(compiledmodel, model_states)
  out_probs <- odin_probs(compiledmodel)
  out_draws <- odin_draws(compiledmodel, model_states)


  format_params <- paste0(names(param_vector), " <- ", param_vector)

  full_odin_text <- c(
    odin_init,
    odin_tot,
    core_eqns,
    out_probs,
    out_draws,
    format_params
  )

  model_generator <- odin::odin(full_odin_text)

  x <- model_generator$new()
  ncols <- length(model_states) + 2

  results_df <- as.data.frame(matrix(nrow = 0, ncol = ncols))
  for (each_sim in seq(1, n_sims)) {
    sim_results <- as.data.frame(x$run(t0:(n_timesteps + t0)))
    sim_results["sim"] <- each_sim
    results_df <- rbind(results_df, sim_results)
  }
  results_df
}
