#' Set popsize for named states
#'
#' @param tblpopsize compiled updatedstate table with popsize value
#' @param namevec named list/vector vector specifying updatedstate popsize
#' values to set. Does not impact unnamed states.
#' @return tibble with columns of updated state names
#' (updatedstates) and current conditions (popsize)
#' @export
#' @importFrom rlang .data
setpopsize_byname <- function(tblpopsize, namevec) {
  check_inputpopsize(tblpopsize)
  intlogic_before <- intlogic(tblpopsize[["popsize"]])

  # check input
  namescompiled <- tblpopsize[["updatedstate"]]
  namesinput <- names(namevec)

  if (length(setdiff(namesinput, namescompiled)) > 0) {
    warning("Input pops contains named elements not found in compiled model.
      These values will be ignored.")
  }

  for (curridx in seq_along(namevec)) {
    currname <- namesinput[[curridx]]
    currval <- namevec[[curridx]]
    currlogic <- tblpopsize[["updatedstate"]] == currname
    tblpopsize[["popsize"]][currlogic] <- currval
  }

  tblpopsize <- maintainpopsizetype(tblpopsize, "popsize")

  intlogic_after <- intlogic(tblpopsize[["popsize"]])

  if (intlogic_before != intlogic_after) {
    warning("popsize changed between non-integer and integer. Ensure this is
    intentional.")
  }

  tblpopsize
}
