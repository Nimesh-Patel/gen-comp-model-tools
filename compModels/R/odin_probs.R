#' Calculate individual transition probabilities for odin model
#' Refer to https://mrc-ide.github.io/odin.dust/articles/sir_models.html
#' for format examples
#'
#' @param compiled_model compiled model object
#' created by compModels function compilemodel()
#' @return vector of strings formatted to describe transition probabilities
#' @family model_building
#' @export
odin_probs <- function(compiled_model) {
  compiled_model$modelinstructions$tblprocesses |>
    dplyr::mutate(
      odin_probs =
        glue::glue("p_{states_down}{states_up} <- 1 - exp(-{percapitarate})")
    ) |>
    dplyr::pull(odin_probs)
}
