#' Convert groupnames input to named list
#'
#' groupnames input can be vector or named list, depending on model complexity.
#' This function coerces input so that it is a named list with names matching
#' table columns that contain input values. The resulting consistent format
#' simplifies manipulating user input across functions.
#'
#' @param groupnames vector or named list of groupnames. values are elements of
#' input table
#' @param tblwgroupcols table that has grouptype columns with groupnames as
#' values
#' @return groupnames coerced into named list with respective grouptype as name
#' @importFrom rlang .data
fixinput_groupnames <- function(groupnames, tblwgroupcols) {
  if (!identical(groupnames, NA)) {
    tblwgroupcols <- tblwgroupcols |>
      dplyr::select(-tidyselect::any_of(c("updatedstate", "popsize"))) |>
      suppressWarnings()
    if (is.null(names(groupnames))) {
      matchingcol_logic <-
        sapply(tblwgroupcols, function(col) any(col %in% groupnames))
      matchingcol <- colnames(tblwgroupcols)[matchingcol_logic]
      if (length(matchingcol) != 1) {
        stop("Either multiple group columns contain input group values. Or none
        contain. Check spelling and/or input explicitly named list")
      }
      groupnames <- as.list(groupnames)
      names(groupnames) <- matchingcol
    }
    groupnames <- as.list(groupnames)
  }
  groupnames
}
