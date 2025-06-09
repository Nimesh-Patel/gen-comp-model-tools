#' Remove popsize for named states
#'
#' @param tblpopsize compiled updatedstate table with popsize value
#' @param namevec named list/vector vector specifying updatedstate popsize
#' values to set. Does not impact unnamed states.
#' @return tibble with columns of updated state names
#' (updatedstates) and current conditions (popsize)
#' @export
#' @importFrom rlang .data
removepopsize_byname <- function(tblpopsize, namevec) {
  tblout <- addpopsize_byname(tblpopsize, -namevec)

  if (TRUE %in% (tblout[["popsize"]] < 0)) {
    stop("Removing population led to negative population.")
  }
  tblout
}
