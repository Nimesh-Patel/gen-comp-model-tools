#'  Run generalized stochastic compartmental model with intervention
#'
#' Runs basic stochastic compartmental model with an intervention using
#' [GillespieSSA::ssa()] given a user specified number and names of
#' compartments, transition matrix between compartments (change matrix),
#' disease-related parameters, number of time steps, initial starting
#' conditions for each compartment (state), and the start and end time for the
#' intervention. No seed is set in the function so results are expected to vary
#' each time the function is run. Parameter values must be a vector, as
#' required by GillespieSSA, and ssa() function returns given compartment names.
#'
#' @param parms_vec vector of parameter values
#' @param propensity_fns the transition rate equations
#' @param init_vals the starting values for populations in each compartment
#' @param n_timesteps number of time steps
#' @param change_matrix matrix that gives transition information between states
#' @param n_sims number of simulations desired
#' @param intervention_start_time time at which intervention will start
#' @param intervention_end_time time at which intervention will end
#' @param modifier_matrix value to be multiplied against transition matrix
#' @return list of dataframes object of the data part of GillespieSSA object
#' containing the time and states of the simulation. Number of elements should
#' match the number of simulations
#' @export
#' @examples
#' \dontrun{
#' initial_vals <- c(s = 990, e = 100, i1 = 10, i2 = 0, r = 0)
#' change_mtx <- matrix(
#'   c(
#'     -1, 0, 0, 0, # S -> E
#'     +1, -1, 0, 0, # E -> I1
#'     0, +1, -1, 0, # I1 -> I2
#'     0, 0, +1, -1, # I2 -> R
#'     0, 0, 0, +1
#'   ),
#'   nrow = length(initial_vals),
#'   byrow = TRUE,
#'   dimnames = list(names(initial_vals), NULL)
#' )
#' # Modifier matrix to reduce beta by half during the intervention period
#' modifier_matrix <- matrix(c(0.5),
#'   nrow = 1,
#'   dimnames = list(c("beta"), NULL)
#' )
#' results <- fn_run_gen_stochastic_int(
#'   parms_vec = c(
#'     beta = 0.001, sigma = 0.1,
#'     alpha = 0.05, gamma = 0.1
#'   ),
#'   propensity_fns = c(
#'     "beta * s * i1",
#'     "sigma * e",
#'     "alpha * i1",
#'     "gamma * i2"
#'   ),
#'   init_vals = initial_vals,
#'   n_timesteps = 100,
#'   change_matrix = change_mtx,
#'   n_sims = 10,
#'   intervention_start_time = 20,
#'   intervention_end_time = 40,
#'   modifier_matrix = modifier_matrix
#' )
#' }
fn_run_gen_stochastic_int <- function(parms_vec, propensity_fns, init_vals,
                                      n_timesteps, change_matrix, n_sims,
                                      intervention_start_time = NULL,
                                      intervention_end_time = NULL,
                                      modifier_matrix = NULL) {
  nu <- change_matrix
  t <- n_timesteps

  sims <- lapply(1:n_sims, function(i) {
    current_parms_vec <- parms_vec

    sim_data <- GillespieSSA::ssa(
      x0 = init_vals,
      a = propensity_fns,
      nu = nu,
      parms = current_parms_vec,
      tf = t,
      method = GillespieSSA::ssa.d(),
      simName = "General Stochastic Model With Intervention"
    )$data

    if (!is.null(intervention_start_time) && !is.null(intervention_end_time)) {
      sim_data_modified <- sim_data

      for (j in seq_len(nrow(sim_data))) {
        current_time <- sim_data[j, "t"]
        next_event_index <- j + 1

        if (is_intervention_period(
          current_time,
          intervention_start_time,
          intervention_end_time
        )) {
          updated_parms_vec <- update_parameters(
            current_parms_vec,
            modifier_matrix
          )
          next_event_time <- ifelse(next_event_index <= nrow(sim_data),
            sim_data[next_event_index, "t"],
            intervention_end_time
          )

          temp_sim_data <- GillespieSSA::ssa(
            x0 = stats::setNames(as.numeric(sim_data[j, -1]), names(init_vals)),
            a = propensity_fns,
            nu = nu,
            parms = updated_parms_vec,
            tf = min(next_event_time - current_time, t - current_time),
            method = GillespieSSA::ssa.d(),
            simName = "General Stochastic Model with Intervention"
          )$data

          temp_sim_data <- as.data.frame(temp_sim_data)
          temp_sim_data$t <- temp_sim_data$t + current_time

          if (next_event_index <= nrow(sim_data)) {
            sim_data_remaining <- sim_data[next_event_index:nrow(sim_data), ]
            sim_data_combined <- rbind(temp_sim_data, sim_data_remaining)
            sim_data_modified <- rbind(
              sim_data_modified[1:j, ],
              sim_data_combined
            )
          } else {
            sim_data_modified <- rbind(sim_data_modified[1:j, ], temp_sim_data)
          }
        } else {
          current_parms_vec <- parms_vec
        }

        if (current_time > intervention_end_time) {
          current_parms_vec <- parms_vec
        }

        if (next_event_index > nrow(sim_data)) {
          break
        }
      }
      return(data.frame(sim_data_modified))
    } else {
      return(data.frame(sim_data))
    }
  })

  lapply(sims, function(x) data.frame(x))
}
