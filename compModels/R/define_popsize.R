#' Set population states for compiled
#' model using all states
#'
#' Sets conditions globally based on current states.
#'
#' @param compiledmodel compiled model output
#' @param inputpops named list/vector vector specifying updatedstate population
#' sizes
#' Unspecified populations are set to 0
#' @return tibble with columns of updated state names
#' (updatedstates) and current conditions (popsize)
#' @export
#' @importFrom rlang .data
define_popsize <- function(compiledmodel, inputpops = c()) {
  # check input
  compiledstates <- compiledmodel$modelinstructions$tblupdatedstates
  namescompiled <- compiledmodel$modeloutstructions$updatedstates

  # preserve order of states
  tblupdatestate <-
    tibble::tibble(updatedstate = namescompiled) |>
    dplyr::left_join(compiledstates, by = "updatedstate")

  if (length(inputpops) > 0) {
    namesinput <- names(inputpops)
    if (length(setdiff(namescompiled, namesinput)) > 0) {
      message("Input pops does not assign values to every state. Setting
        others to 0")
    }
    if (length(setdiff(namesinput, namescompiled)) > 0) {
      message("Input pops contains named elements not found in compiled model.
        These values will be ignored.")
    }

    tblright <- tibble::enframe(inputpops) |>
      dplyr::rename(updatedstate = "name", popsize = "value")

    tblout <- tblupdatestate |>
      dplyr::left_join(tblright, by = "updatedstate") |>
      tidyr::replace_na(list(popsize = 0))
  } else {
    message("No population sizes input. Initializing all states to 0.")
    tblout <- tblupdatestate |> dplyr::mutate(popsize = 0)
  }
  # convert to integer if consistent with input
  tblout <- maintainpopsizetype(tblout, "popsize")

  tblout
}
