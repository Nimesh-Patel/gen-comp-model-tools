#' Set population states for compiled
#' model using all states
#'
#' Sets conditions globally based on current states.
#'
#' @param outlist compiled model output
#' @param namedvector numeric vector specifying initial populations
#' for named basestates
#' default is c() which specifies all 0 population
#' for each basestate
#' @return tibble with columns of updated state names
#' (updatedstates) and current conditions (X0)
#' @export
#' @importFrom rlang .data
define_currentstate <- function(outlist, namedvector = c()) {
  tblupdatestate <- tibble::tibble(
    updatedstate =
      outlist$modeloutstructions$updatedstates
  ) |>
    dplyr::left_join(outlist$modelinstructions$tblupdatedstates |>
      dplyr::select(
        -.data$interactionscale,
        -.data$transitionscale,
        -.data$environment_names
      ) |>
      dplyr::distinct(.data$updatedstate,
        .keep_all = TRUE
      ), by = "updatedstate")
  tblout <- tblupdatestate |> dplyr::mutate(X0 = 0)
  if (length(namedvector) > 0) {
    currstates <- tblupdatestate |>
      dplyr::distinct(.data$updatedstate) |>
      dplyr::pull()
    # check that the input names are uniquely defined states
    if (FALSE %in% (names(namedvector) %in% currstates)) {
      stop("condition vector must have names that match current states")
    }
    if (length(names(namedvector)) != length(unique(names(namedvector)))) {
      stop("current states must be uniquely")
    }
    jointhistbl <- tibble::enframe(namedvector) |>
      dplyr::rename(updatedstate = .data$name, X0 = .data$value)
    tblout <- tblupdatestate |>
      dplyr::left_join(jointhistbl, by = "updatedstate") |>
      tidyr::replace_na(list(X0 = 0))
  } else {
    tblout <- tblout |> dplyr::mutate(X0 = 0)
  }
  return(tblout)
}
