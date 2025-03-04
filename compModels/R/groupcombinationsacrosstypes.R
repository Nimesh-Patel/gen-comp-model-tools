#' Generates all combinations groups with different grouptypes
#'
#' Useful when stratifying populations in multiple ways.
#'
#' @param combinethesetypes list of character vectors that
#' specify grouptypes to combine groups.
#' @param tblgroup tibble specifying groups in each group type
#' @return tibble where each row is a combination of groups
#' with grouptype columns and scaled process rates
#' @importFrom rlang .data
#' @importFrom rlang :=
groupcombinationsacrosstypes <- function(combinethesetypes, tblgroup) {
  if (length(combinethesetypes) > 1) {
    currtblgroup <- tblgroup |>
      dplyr::filter(.data$grouptype %in% combinethesetypes)
    currtblgrouptypes <- currtblgroup |>
      dplyr::distinct(.data$grouptype) |>
      dplyr::pull()
    listoftypes <- currtblgroup |> dplyr::group_split(.data$grouptype)
    names(listoftypes) <- currtblgrouptypes
    currbasestates <- currtblgroup$basestates[[1]]
    # types to get combinations
    list2cross <- lapply(listoftypes, function(x) {
      x |> dplyr::pull("groupname")
    })
    # multiply transition scales
    trans2cross <- lapply(listoftypes, function(x) {
      x |> dplyr::pull(.data$transitionscale)
    })
    # multiply interaction scales
    int2cross <- lapply(listoftypes, function(x) {
      x |> dplyr::pull(.data$interactionscale)
    })
    # eventually make this for any given name scale
    crossframe <- do.call(expand.grid, list2cross)
    trans2crossframe <- do.call(expand.grid, trans2cross) |>
      tidyr::unite("transitionscale", sep = "")
    int2crossframe <- do.call(expand.grid, int2cross) |>
      tidyr::unite("interactionscale", sep = "")

    tblout <- dplyr::bind_cols(crossframe, trans2crossframe, int2crossframe) |>
      dplyr::mutate(basestates = list(currbasestates))
  } else {
    tblout <- tblgroup |>
      dplyr::filter(.data$grouptype == combinethesetypes) |>
      dplyr::select(
        "groupname", "interactionscale",
        "transitionscale", "basestates"
      ) |>
      dplyr::rename(!!combinethesetypes := .data$groupname)
  }
  return(tblout)
}
