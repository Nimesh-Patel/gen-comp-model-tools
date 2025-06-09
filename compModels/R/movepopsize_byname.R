#' Move popsize for named states
#'
#' Move populations between some number of named updatedstates and to another
#' single named state.
#'
#' @param tblpopsize compiled updatedstate table with popsize column
#' @param popsizeval value to move between states
#' @param namebefore updatedstate value(s) in tblupdatedstate to remove popsize
#' @param nameafter updatedstate value in tblupdatedstate to add popsize
#' @return tibble with columns of updated state names
#' (updatedstates) and current conditions (popsize)
#' @export
#' @importFrom rlang .data
movepopsize_byname <- function(tblpopsize, popsizeval, namebefore, nameafter) {
  check_inputpopsize(tblpopsize)
  if (length(nameafter) != 1) {
    stop("nameafter argument must be length 1")
  }

  intlogic_before <- intlogic(tblpopsize[["popsize"]])

  namescompiled <- tblpopsize[["updatedstate"]]

  if (FALSE %in% (c(namebefore, nameafter) %in% namescompiled)) {
    stop("Input names not in updatedstates column. Check spelling.")
  }
  logicbefore <- tblpopsize[["updatedstate"]] %in% namebefore
  tblpopsize[["popsize"]][logicbefore] <-
    tblpopsize[["popsize"]][logicbefore] - popsizeval

  logicafter <- tblpopsize[["updatedstate"]] == nameafter
  tblpopsize[["popsize"]][logicafter] <-
    tblpopsize[["popsize"]][logicafter] + length(namebefore) * popsizeval

  # check for negative populations
  if (min(tblpopsize[["updatedstate"]]) < 0) {
    stop("Moving populations caused negative population. Consider using
    movefractionpopsize_byname.")
  }
  tblpopsize <- maintainpopsizetype(tblpopsize, "popsize")

  intlogic_after <- intlogic(tblpopsize[["popsize"]])

  if (intlogic_before != intlogic_after) {
    message("popsize changed between non-integer and integer. Ensure this is
    intentional.")
  }

  tblpopsize
}
