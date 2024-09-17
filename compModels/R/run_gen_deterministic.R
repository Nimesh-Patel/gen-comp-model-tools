#' Run generalized deterministic compartmental model
#'
#' Runs basic compartmental model using [deSolve::ode()] given a user specified
#' number and names of compartments, transition matrix between compartments,
#' disease-related parameters, times, and starting conditions
#'
#' @param init_vals the starting values for populations in each compartment

#' @param times vector of time steps to run the model over
#' @param comp_names names of the compartments
#' @param trans_matrix matrix of transition rates between each of the
#' compartments
#' @param intervention_start_time time at which intervention will start
#' @param intervention_end_time time at which intervention will end
#' @param modifier_matrix value to be multiplied against transition matrix
#' @return data frame with columns for time and each compartment
#' @export
#' @examples
#' \dontrun{
#' initial_vals <- c(A = 1e5 - 1, B = 1, C = 0, D = 0, E = 0)
#' time_seq <- seq(0, 100, length.out = 101)
#' compartment_nms <- c("A", "B", "C", "D", "E")
#' transition_mtx <- matrix(
#'   c(
#'     -0.1, 0.1, 0, 0, 0,
#'     0, -0.1, 0.1, 0, 0,
#'     0, 0, -0.1, 0.1, 0,
#'     0, 0, 0, -0.1, 0.1,
#'     0.1, 0, 0, 0, -0.1
#'   ),
#'   nrow = length(compartment_nms),
#'   byrow = TRUE,
#'   dimnames = list(compartment_nms, compartment_nms)
#' )
#' start <- 25
#' end <- 30
#' modify <- 0.5
#' mod_mat <- matrix(
#'   c(
#'     1, 0.5, 1, 1, 1,
#'     1, 0.5, 0.25, 1, 1,
#'     1, 1, 0.25, .1, 1,
#'     1, 1, 1, .1, 1,
#'     1, 1, 1, 1, 1
#'   ),
#'   nrow = length(compartment_nms),
#'   byrow = TRUE,
#'   dimnames = list(compartment_nms, compartment_nms)
#' )
#' modelout <- run_gen_deterministic(
#'   initial_vals, time_seq, compartment_nms,
#'   transition_mtx,
#'   start, end,
#'   mod_mat
#' )
#' }
run_gen_deterministic <- function(init_vals, times, comp_names, trans_matrix,
                                  intervention_start_time = NULL,
                                  intervention_end_time = NULL,
                                  modifier_matrix = NULL) {
  validate_gen_determ_input(init_vals, times, comp_names, trans_matrix)

  model_ode <- function(time, state, parameters) {
    with(as.list(c(state)), {
      # Determine whether to apply intervention based on time
      if (is_intervention_period(
        time, intervention_start_time,
        intervention_end_time
      )) {
        current_trans_matrix <- modify_trans_mtx(
          trans_matrix,
          modifier_matrix
        )
      } else {
        current_trans_matrix <- trans_matrix
      }

      # Calculate changes based on current transition matrix
      # (original or modified)
      change_rates <- calculate_change_rates(
        state, comp_names,
        current_trans_matrix
      )

      list(change_rates)
    })
  }

  parms <- list()

  # Set initial states with names from comp_names argument
  init_state <- stats::setNames(init_vals, comp_names)

  # Solve ODEs using lsoda from deSolve package
  modelout <- deSolve::ode(
    y = init_state, times = times,
    func = model_ode, parms = parms
  )

  as.data.frame(modelout)
}
