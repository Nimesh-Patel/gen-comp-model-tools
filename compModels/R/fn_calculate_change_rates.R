#' @title calculate_change_rates
#' @description
#' Helper function to calculate change rates for in and out flows in the
#' compartmental model
#' @param state holds the current state of the compartments
#' @param comp_names names of the compartments
#' @param trans_matrix matrix of transition rates between each of the
#' compartments
#' @return calculated change rates
#' @export
#' @examples
#' \dontrun{
#' change_rates <- calculate_change_rates(
#'   state, comp_names,
#'   current_trans_matrix
#' )
#' }
calculate_change_rates <- function(state, comp_names, trans_matrix) {
  change_rates <- numeric(length(state))
  for (i in seq_along(comp_names)) {
    inflow <- sum(state * trans_matrix[, i], na.rm = TRUE)
    outflow <- sum(state[i] * trans_matrix[i, ], na.rm = TRUE)
    change_rates[i] <- inflow - outflow
  }
  change_rates
}
