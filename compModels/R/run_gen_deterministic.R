#' Run generalized deterministic compartmental model
#'
#' Runs basic compartmental model using [deSolve::ode()] given a user specified
#' number and names of compartments, transition matrix between compartments,
#' disease-related parameters, times, and starting conditions
#'
#' @param init_vals the starting values for populations in each compartment
#' @param times vector (type: double) containing time steps to run the model
#' over
#' @param comp_names names of the compartments
#' @param trans_matrix matrix of transition rates between each of the
#' compartments
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
#' modelout <- run_gen_deterministic(
#'   initial_vals, time_seq, compartment_nms,
#'   transition_mtx
#' )
#' }
run_gen_deterministic <- function(init_vals, times, comp_names, trans_matrix) {
  validate_gen_determ_input(init_vals, times, comp_names, trans_matrix)

  model_ode <- function(time, state, parms) {
    # Passing parms above so that this function works with ode(),
    # the parms are passed to ode() by their existence in the
    # enclosed environment of run_comp_model2()

    # Initialize change rates
    change_rates <- numeric(length(state))

    # Calculate changes based on transition matrix
    for (i in seq_along(comp_names)) {
      inflow <- sum(state * trans_matrix[, i], na.rm = TRUE)
      outflow <- sum(state[i] * trans_matrix[i, ], na.rm = TRUE)
      change_rates[i] <- inflow - outflow
    }

    list(change_rates)
  }

  # Set initial states with names from comp_names argument
  init_state <- stats::setNames(init_vals, comp_names)

  # Solve ODEs using lsoda from deSolve package
  modelout <- deSolve::ode(
    y = init_state, times = times,
    func = model_ode
  )

  as.data.frame(modelout)
}
