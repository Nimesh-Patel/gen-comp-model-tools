#' Checks validity of input to generalized compartmental model
#'
#' Helper function that runs checks on generalized compartmental model input
#' parameters, times, and initial conditions
#'
#' @param init_vals the starting values for populations in each compartment
#' @param times vector of time steps to run the model over
#' @param comp_names names of the compartments
#' @param trans_matrix matrix of transition rates between each of the
#' compartments
#' @param subgroups_list list of one or more vectors of subgroup names to be
#' used for stratification
#' @param intervention_start_time time at which intervention will start
#' @param intervention_end_time time at which intervention will end
#' @param modifier single value, matrix, or function to be applied to transition
#' matrix during an intervention
#' @return stops with information or indicates Input Checks Passed
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
