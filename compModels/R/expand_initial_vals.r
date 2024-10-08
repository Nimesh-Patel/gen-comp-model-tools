#' Helper function to expand initial values based on subgroup combinations
#'
#' Takes in a vector of initial values for the user defined names and number of
#' model states (compartments) and determines all the combinations of states and
#' user defined subgroups and distributes individuals evenly across all
#' combinations.
#'
#' @param init_vals_base the initial values supplied by the user for each state
#' (compartment) before stratification
#' @param num_combinations the total number of all subgroup combinations
#' @return an expanded vector of initial values for each state (compartment) and
#' subgroup combination
expand_init_vals <- function(init_vals_base, num_combinations) {
  if (num_combinations == 0) {
    return(init_vals_base)
  }
  base_allocations <- init_vals_base %/% num_combinations
  remainders <- init_vals_base %% num_combinations
  expanded_init_vals <- rep(base_allocations, times = num_combinations)
  for (i in seq_along(remainders)) {
    indices <-
      ((i - 1) * num_combinations + 1):
      ((i - 1) * num_combinations + remainders[i])
    expanded_init_vals[indices] <- expanded_init_vals[indices] + 1
  }
  expanded_init_vals
}
