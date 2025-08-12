#' Create core equations for odin model
#'
#' Core equations describes the movements between states in the model
#'
#' Format needs to be: update(state) <- state - n_exit + n_enter
#'
#' @param compiled_model compiled model object
#' created by compModels function compilemodel()
#' @param states_vec vector of states in model as character strings
#' example: base_states <- c("S", "E", "I", "R")
#' @return vector of strings formatted to describe core equation
#' @family model_building
#' @export
odin_core_eqns <- function(compiled_model, states_vec) {
  get_ops <- function(state) {
    compiled_model$modelinstructions$tblprocesses |>
      dplyr::filter(states_up == state | states_down == state) |>
      glue::glue_data(
        "{ifelse(states_up == state, '+', '-')} n_{states_down}{states_up}"
      ) |>
      glue::glue_collapse(sep = " ")
  }
  purrr::map_chr(
    states_vec,
    ~ glue::glue("update({.}) <- {.} {get_ops(.)}")
  )
}
