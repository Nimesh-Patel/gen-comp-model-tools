#' Set initial population states for compiled
#' model using base states
#'
#' Sets globally based on base states. Update for more flexibility.
#'
#' @param outlist compiled model output
#' @param namedvector numeric vector specifying inital populations
#' for named basestates
#' default is c() which specifies all 0 population
#' for each basestate
#' @return tibble with columns of updated state names
#' (updatedstates) and initial conditions (X0)
#' @export
#' @importFrom rlang .data
define_initialstate <- function(outlist, namedvector = c()) {
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
    currbasestates <- tblupdatestate |>
      dplyr::distinct(.data$basestates) |>
      dplyr::pull()
    # check that the input names are uniquely defined basestates
    if (FALSE %in% (names(namedvector) %in% currbasestates)) {
      stop("initial condition vector must have names that match basestates")
    }
    if (length(names(namedvector)) != length(unique(names(namedvector)))) {
      stop("basestates must be uniquely")
    }
    jointhistbl <- tibble::enframe(namedvector) |>
      dplyr::rename(basestates = .data$name, X0 = .data$value)
    tblout <- tblupdatestate |>
      dplyr::left_join(jointhistbl, by = "basestates") |>
      tidyr::replace_na(list(X0 = 0))
  } else {
    tblout <- tblout |> dplyr::mutate(X0 = 0)
  }
  return(tblout)
}
