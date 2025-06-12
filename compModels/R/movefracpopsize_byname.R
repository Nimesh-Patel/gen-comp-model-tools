#' Move fraction of popsize for named states
#'
#' Move fraction of populations between some number of named updatedstates to
#' another single named state.
#'
#' @param tblpopsize compiled updatedstate table with popsize column
#' @param fracpopsize fraction of population to move between states
#' @param namebefore updatedstate value(s) in tblupdatedstate to remove popsize
#' @param nameafter updatedstate value in tblupdatedstate to add popsize
#' @param samplerandomly logical specifying whether population moved is
#' determined by binomial sampling (TRUE and integer populations) or by the
#' fractional amount (FALSE) with rounding when integer populations.
#' @return tibble with columns of updated state names
#' (updatedstates) and current conditions (popsize)
#' @export
#' @importFrom rlang .data
movefracpopsize_byname <- function(tblpopsize,
                                   fracpopsize,
                                   namebefore,
                                   nameafter,
                                   samplerandomly) {
  check_inputpopsize(tblpopsize)
  if (length(nameafter) != 1) {
    stop("nameafter argument must be length 1")
  }

  intlogic_before <- intlogic(tblpopsize)

  if (is.na(samplerandomly)) {
    samplerandomly <- intlogic_before
  }

  if (!intlogic_before && samplerandomly) {
    stop("Random sampling of moving populations specified, but populations
      are non-integer. Set samplerandomly to FALSE to move exact fraction.")
  }

  namescompiled <- tblpopsize[["updatedstate"]]

  if (FALSE %in% (c(namebefore, nameafter) %in% namescompiled)) {
    stop("Input names not in updatedstates column. Check spelling.")
  }

  logicbefore <- tblpopsize[["updatedstate"]] %in% namebefore
  popbefore <- tblpopsize[["popsize"]][logicbefore]
  if (samplerandomly) {
    movepopsize <-
      sapply(popbefore, function(x) {
        sum(stats::rbinom(x, 1, fracpopsize))
      })
  } else {
    movepopsize <- popbefore * fracpopsize
    if (intlogic_before) {
      movepopsize <- round(movepopsize)
    }
  }

  tblpopsize[["popsize"]][logicbefore] <-
    tblpopsize[["popsize"]][logicbefore] - movepopsize

  logicafter <- tblpopsize[["updatedstate"]] == nameafter
  tblpopsize[["popsize"]][logicafter] <-
    tblpopsize[["popsize"]][logicafter] + sum(movepopsize)

  # check for negative populations
  if (min(tblpopsize[["updatedstate"]]) < 0) {
    stop("Moving populations caused negative population. Consider using
    movefractionpopsize_byname.")
  }

  intlogic_after <- intlogic(tblpopsize)

  if (intlogic_before != intlogic_after) {
    message("popsize changed between non-integer and integer. Ensure this is
    intentional.")
  }

  tblpopsize
}
