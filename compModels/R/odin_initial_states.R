#' Create initial states for odin model
#'
#' Format needs to be initial(state) <- value for each value
#' Need to provide a vector of states and a named vector of initial values.
#' If a state is provided with no value, it will be set to 0.
#'
#' @param states_vec vector of states in model as character strings
#' example: base_states <- c("S", "E", "I", "R")
#' @param initial_states_vec named vector of initial values for each state
#' example: initial_states_vec <- c("S" = 999, "E" = 0, "I" = 1, "R" = 0)
#' @return vector of strings formatted for odin initial states
#' @family initial_conditions
#' @export
odin_initial_states <- function(states_vec, initial_states_vec) {
  glue::glue(
    "initial({states_vec}) <-",
    " {ifelse(states_vec %in% names(initial_states_vec),",
    " initial_states_vec[states_vec], 0)}"
  )
}
