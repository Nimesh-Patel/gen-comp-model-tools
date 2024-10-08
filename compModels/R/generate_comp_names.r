#' Helper function to generate compartment names for all subgroup combinations
#'
#' Takes in a vector of initial values for the user defined names of model
#' states (compartments) and determines all the name combinations of states and
#' user defined subgroups
#'
#' @param comp_names_base the initial values of compartment names supplied by
#' the user for each state (compartment) before stratification
#' @param subgroup_combinations a data frame containing all combinations of
#' subgroups for stratification
#' @return an expanded vector of compartment names for each state (compartment)
#' and subgroup combination.
generate_comp_names <- function(comp_names_base, subgroup_combinations) {
  if (nrow(subgroup_combinations) == 0) {
    return(comp_names_base)
  }

  full_comp_names <- c()

  for (subgroup_idx in seq_len(nrow(subgroup_combinations))) {
    current_subgroup <- subgroup_combinations[subgroup_idx, ]
    current_subgroup_str <- paste(current_subgroup, collapse = "_")

    names_for_current_subgroup <- paste(comp_names_base, current_subgroup_str,
      sep = "_"
    )
    full_comp_names <- c(full_comp_names, names_for_current_subgroup)
  }
  full_comp_names
}
