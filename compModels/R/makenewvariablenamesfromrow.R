#' Generates new state names from basestates
#'
#' Pastes basestate names with column names and column
#' values, delimited by _
#'
#' @param tibblerow tibble row to generate state names from
#' Ignore columns with all NA.
#' @return string of updated name
#' @importFrom rlang .data
makenewvariablenamesfromrow <- function(tibblerow) {
  currbasestate <- tibblerow |> dplyr::pull(.data$basestates)
  newname <- currbasestate
  # remove columns with only NA
  tibblerow <- tibblerow |>
    dplyr::select(-.data$basestates) |>
    dplyr::select_if(~ !any(is.na(.)))
  if (ncol(tibblerow) >= 1) {
    for (currcol in colnames(tibblerow)) {
      newname <- paste0(
        newname, "_",
        currcol,
        as.character(tibblerow |> dplyr::pull(currcol))
      )
    }
  }
  return(newname)
}
