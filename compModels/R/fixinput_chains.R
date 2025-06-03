#' Convert chains input to named list
#'
#' chains input can be vector or named list, depending on model complexity.
#' This function coerces input so that it is a named list with names matching
#' table columns that contain input values. The resulting consistent format
#' simplifies manipulating user input across functions.
#'
#' @param chains vector or named list of chains values are elements of
#' input table
#' @param tblwgroupcols table that has chainid columns with chains as
#' values
#' @return chains coerced into named list with respective chainid as name
#' @importFrom rlang .data
fixinput_chains <- function(chains, tblwgroupcols) {
  if (!identical(chains, NA)) {
    tblwgroupcols <- tblwgroupcols |>
      dplyr::select(-tidyselect::any_of(c("updatedstate", "popsize"))) |>
      suppressWarnings()
    tblwgroupcols_colnames <- names(tblwgroupcols)
    logic_tblchaincols <- grepl("chainid", tblwgroupcols_colnames)

    chain_cols <- tblwgroupcols_colnames[logic_tblchaincols]

    tblwgroupcols <- tblwgroupcols[, chain_cols]

    if (ncol(tblwgroupcols) == 0) {
      stop("No chain columns in input table")
    }

    tblwgroupcols_nona <- tblwgroupcols |>
      dplyr::select(tidyselect::where(~ !all(is.na(.))))
    if (ncol(tblwgroupcols_nona) == 0) {
      message("Chains exist, but no chain step values in this particular
    subsetted tblpopsize. Ignoring filtering by chain steps here.")
      chains <- NA
    } else {
      if (is.null(names(chains))) {
        matchingcol_logic <-
          sapply(
            tblwgroupcols,
            function(col) any(col %in% as.character(chains))
          )
        matchingcol <- chain_cols[matchingcol_logic]
        if (length(matchingcol) != 1) {
          stop("Either multiple or no chain columns detected
            based on values. Check
        input or use explicit named list to stratify by chain steps.")
        }
        chains <- as.list(chains)
        names(chains) <- matchingcol
      }
      chains <- lapply(chains, as.character)
    }
  }
  chains
}
