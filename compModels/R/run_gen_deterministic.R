#' Run generalized deterministic compartmental model
#'
#' Runs basic compartmental model using [deSolve::ode()] given a user specified
#' number and names of compartments, transition matrix between compartments,
#' disease-related parameters, times, and starting conditions. Uses helper
#' function [compModels::gen_rates()] to calculate modified transition rates
#' due to intervention, if intervention is present.
#'
#' @param init_vals_base the starting values for populations in each compartment
#' @param times vector of time steps to run the model over
#' @param comp_names_base names of the compartments
#' @param trans_matrix_base matrix of transition rates between each of the
#' compartments
#' @param subgroups_list list of one or more vectors of subgroup names to be
#' used for stratification
#' @param intervention_start_time time at which intervention will start
#' @param intervention_end_time time at which intervention will end
#' @param modifier single value, matrix, or function to be applied to transition
#' matrix during an intervention
#' @return data frame with columns for time and each compartment
#' @export
#' @examples
#' \dontrun{
#' # Example 1: Generalized model with intervention
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
#'   subgroups_list = NULL,
#'   transition_mtx,
#'   start, end,
#'   mod_mat
#' )
#'
#' # Example 2: With Stratifications
#' init_vals_base <- c(s = 99000, e = 0, i = 10000, r = 0)
#' comp_names_base <- c("s", "e", "i", "r")
#' beta <- .9 # transmission
#' sigma <- .1 # progression E->I
#' gamma <- 0.05 # recovery
#' trans_matrix_base <- matrix(
#'   c(
#'     -beta, 0, 0, 0,
#'     beta, -sigma, 0, 0,
#'     0, sigma, -gamma, 0,
#'     0, 0, gamma, 0
#'   ),
#'   nrow = length(comp_names_base), byrow = TRUE,
#'   dimnames = list(comp_names_base, comp_names_base)
#' )
#' times <- seq(0, 100, by = 1)
#' age_grps_vector <- c("child", "adult", "elderly")
#' vacc_status_vector <- c("vaccinated", "unvaccinated")
#' gender_cat_vector <- c("male", "female")
#' # Modification options
#' single_modifier_value <- 1
#' base_size_modifier_mtx <- matrix(single_modifier_value,
#'   nrow = length(comp_names_base),
#'   ncol = length(comp_names_base),
#'   dimnames = list(
#'     comp_names_base,
#'     comp_names_base
#'   )
#' )
#' get_time_dep_modifier <- function(time, subgroup_combinations) {
#'   num_subgroups <- nrow(subgroup_combinations)
#'   size <- num_subgroups * length(comp_names_base)
#'   if (time < 30) {
#'     modifier_value <- 0.8
#'   } else {
#'     modifier_value <- 0.5
#'   }
#'   matrix(modifier_value, nrow = size, ncol = size)
#' }
#' output1 <- run_gen_deterministic(
#'   init_vals_base = init_vals_base,
#'   times = times,
#'   comp_names_base = comp_names_base,
#'   trans_matrix_base = trans_matrix_base,
#'   subgroups_list =
#'     list(
#'       age_groups = age_grps_vector,
#'       vaccination_statuses = vacc_status_vector,
#'       gender_categories = gender_cat_vector
#'     ),
#'   intervention_start_time = 20,
#'   intervention_end_time = 40,
#'   modifier = get_time_dep_modifier
#' )
#' head(output1)
#' plot_determin_model(output1, stratify_by = age_grps_vector)
#'
#' output2 <- run_gen_deterministic(
#'   init_vals_base = init_vals_base,
#'   times = times,
#'   comp_names_base = comp_names_base,
#'   trans_matrix_base = trans_matrix_base,
#'   subgroups_list =
#'     list(
#'       age_groups = age_grps_vector,
#'       vaccination_statuses = vacc_status_vector,
#'       gender_categories = gender_cat_vector
#'     ),
#'   intervention_start_time = NULL,
#'   intervention_end_time = NULL,
#'   modifier = NULL
#' )
#' head(output2)
#' plot_determin_model(output2, stratify_by = c("child", "adult", "elderly"))
#' }
run_gen_deterministic <- function(init_vals_base, times, comp_names_base,
                                  trans_matrix_base,
                                  subgroups_list = list(),
                                  intervention_start_time = NULL,
                                  intervention_end_time = NULL,
                                  modifier = NULL) {
  validate_gen_determ_input(
    init_vals_base, times, comp_names_base,
    trans_matrix_base, subgroups_list,
    intervention_start_time,
    intervention_end_time,
    modifier
  )

  subgroups_list <- Filter(Negate(is.null), subgroups_list)
  subgroup_combinations <- expand.grid(subgroups_list)
  subgroup_combinations[] <- lapply(subgroup_combinations, as.character)
  comp_names <- generate_comp_names(
    comp_names_base,
    subgroup_combinations
  )
  init_vals <- expand_init_vals(
    init_vals_base,
    nrow(subgroup_combinations)
  )
  names(init_vals) <- comp_names
  trans_matrix <- expand_trans_mtx(
    trans_matrix_base,
    nrow(subgroup_combinations)
  )
  rownames(trans_matrix) <- colnames(trans_matrix) <- comp_names

  parms <- list(
    comp_names = comp_names,
    trans_matrix = trans_matrix,
    intervention_start_time = intervention_start_time,
    intervention_end_time = intervention_end_time,
    modifier = modifier,
    subgroup_combinations = subgroup_combinations
  )

  init_state <- stats::setNames(init_vals, comp_names)

  modelout <- deSolve::ode(
    y = init_state, times = times,
    func = gen_rates, parms = parms
  )

  as.data.frame(modelout)
}
