#' Generate column for petersen matrix for a process
#'
#' Input up/down states can be repeated to
#' increase populations
#' @param statesdown State names to decrease in value
#' (by one)
#' @param statesup statesup State names to increase
#' in value (by one)
#' @param statelist character vector of state names
#' that specifies order in output matrix
#' @return Column matrix with length(statelist) rows and
#' values the total integer changes to population.
states2petercolumn <- function(statesdown, statesup, statelist) {
  outcol <- Matrix::Matrix(0, nrow = length(statelist), ncol = 1, sparse = TRUE)
  for (currstate in statesdown) {
    outcol[statelist == currstate, ] <- outcol[statelist == currstate, ] - 1
  }
  for (currstate in statesup) {
    outcol[statelist == currstate, ] <- outcol[statelist == currstate, ] + 1
  }
  return(outcol)
}
