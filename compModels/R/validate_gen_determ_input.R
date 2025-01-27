#' Checks validity of input to generalized compartmental model
#'
#' Helper function that runs checks on generalized compartmental model input
#' parameters, times, and initial conditions
#'
#' @param init_vals the starting values for populations in each compartment
#' @param times vector of time steps to run the model over
#' @param comp_names names of the compartments
#' @param trans_matrix matrix of transition rates between each of the
#' @return stops with information or indicates Input Checks Passed
#' @examples
#' \dontrun{
#' validate_gen_determ_input(
#'   init_vals = c(1e05 - 1, 10, 5, 1, 0),
#'   times = seq(0.1, 100, by = 0.1),
#'   comp_names = c("A", "B", "C", "D", "E"),
#'   trans_matrix <- matrix(
#'     c(
#'       -0.1, 0.1, 0, 0, 0,
#'       0, -0.1, 0.1, 0, 0,
#'       0, 0, -0.1, 0.1, 0,
#'       0, 0, 0, -0.1, 0.1,
#'       0.1, 0, 0, 0, -0.1
#'     ),
#'     nrow = length(c("A", "B", "C", "D", "E")),
#'     byrow = TRUE,
#'     dimnames = list(c("A", "B", "C", "D", "E"), c("A", "B", "C", "D", "E"))
#'   )
#' )
#' }
validate_gen_determ_input <- function(init_vals, times,
                                      comp_names, trans_matrix, subgroups_list,
                                      intervention_start_time,
                                      intervention_end_time,
                                      modifier) {
  check_init_vals(init_vals, comp_names)
  check_times(times)
  check_trans_mtx(trans_matrix, comp_names)
  check_subgroups_list(subgroups_list)
  check_intervention_times(intervention_start_time, intervention_end_time)
  check_modifier2(modifier)

  print("All input checks passed")
}
