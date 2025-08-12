#' Draws distributions for numbers changing between compartments/states
#' Refer to https://mrc-ide.github.io/odin.dust/articles/sir_models.html
#' for format examples
#'
#' In progress: need to add functionality for multinomial draws
#'
#' @param compiled_model compiled model object
#' created by compModels function compilemodel()
#' @param states_vec vector of states in model as character strings
#' example: base_states <- c("S", "E", "I", "R")
#' @return vector of strings formatted to describe draws of numbers
#' changing between states
#' @family model_building
#' @export
odin_draws <- function(compiled_model, states_vec) {
  get_draws <- function(state) {
    d <- compiled_model$modelinstructions$tblprocesses |>
      dplyr::filter(states_down == state)

    stopifnot("multinomial draws not yet implemented" = NROW(d) <= 1)

    d |>
      glue::glue_data(
        "n_{states_down}{states_up} <-",
        " rbinom({states_down}, p_{states_down}{states_up})"
      )
  }

  # unlist because map_chr gets mad if get_draws returns chr(0)
  # which it does in the zero case
  unlist(purrr::map(states_vec, get_draws))
}
