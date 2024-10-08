#' Helper function to apply intervention modifiers to an expanded transition
#' matrix
#'
#' Creates a modifier matrix of values of the same dimension as the current
#' transition matrix based on the user supplied modifier function, and
#' multiplies the modifier and transition matrix to produce a matrix of modified
#' transition rates between states (compartments).
#'
#' @param expanded_trans_matrix an expanded matrix of transition rates between
#' for each state (compartment) and subgroup combination
#' @param modifier user supplied function or matrix to apply modifications to
#' transition matrix
#' @param time vector of time steps to run the model over
#' @param subgroup_combinations a data frame containing all combinations of
#' subgroups for stratification
#' @return a matrix of modified transition rates between states (compartments)
apply_modifiers <- function(expanded_trans_matrix, modifier,
                            time, subgroup_combinations) {
  modified_trans_mtx <- expanded_trans_matrix
  if (!is.null(modifier)) {
    if (is.function(modifier)) {
      print("Running modifier function")
      current_modifier_matrix <- modifier(time, subgroup_combinations)
    } else if (is.matrix) {
      print("Using provided modifier matrix")
      current_modifier_matrix <- modifier
    } else if (is.numeric(modifier) && length(modifier) == 1) {
      print("Using provided modifier value")
      current_modifier_matrix <- matrix(modifier,
        nrow = nrow(expanded_trans_matrix),
        ncol = ncol(expanded_trans_matrix)
      )
    }

    if (all(dim(current_modifier_matrix) == dim(expanded_trans_matrix))) {
      modified_trans_mtx <- expanded_trans_matrix * current_modifier_matrix
    } else {
      stop("Modifier matrix dimensions do not match transition
            matrix dimensions.")
    }
  }
  modified_trans_mtx
}
