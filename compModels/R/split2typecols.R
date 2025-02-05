#' Expand table to include group name combinations
#'
#' Helper function simplifying combining group names
#'
#' @param currtbl tibble to expand
#' @param combinelist list with character vector values
#' specifying grouptypes to make group combinations
#' @param listgroupcombinations the combined
#' groupnames based on group types
#' @return tibble
#' @importFrom rlang .data
split2typecols <- function(currtbl, combinelist, listgroupcombinations) {
  currtypelist <- currtbl$grouptype[[1]]
  # look through combo list and grab named list
  namedlistidx <- NA
  for (combinelistidx in seq_along(combinelist)) {
    currcombinelist <- combinelist[[combinelistidx]]
    if (!(FALSE %in% (currtypelist %in% currcombinelist))) {
      namedlistidx <- combinelistidx
    }
  }
  if (is.na(namedlistidx)) {
    stop("no match in group combo list")
  }

  listgroupcombo_unnest <- listgroupcombinations[[namedlistidx]] |>
    tidyr::unnest(.data$basestates)
  # use to loop
  currtbl2loop <- currtbl |>
    dplyr::select(-.data$grouptype)
  # loop because they're easier, clean later
  bindlist <- list()
  for (rowidx in seq_len(nrow(currtbl2loop))) {
    currtblrow <- currtbl2loop[rowidx, ]
    currfromstate <- currtblrow$fromstate[[1]]
    currlistgroupcombo <- listgroupcombo_unnest |>
      dplyr::filter(.data$basestates == currfromstate)
    bindlist[[rowidx]] <- dplyr::bind_cols(currtblrow, currlistgroupcombo)
  }
  tblout <- dplyr::bind_rows(bindlist)
  return(tblout)
}
