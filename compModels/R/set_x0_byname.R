#' Generate initial conditions specifying value
#' for state names
#'

#' @param tblx0 tibble with "updatedstate" column
#' of state names
#' @param inputvec named list with state names as names and
#' values for population abundances
#' @return Same tibble with "X0" column added with values
set_x0_byname <- function(tblx0, inputvec) {
  if (length(inputvec) == 0) {
    stop("list specifying initial conditions can't be empty")
  }
  if (FALSE %in% (names(inputvec) %in% tblx0$updatedstate)) {
    stop("input names must match updatedstate")
  }
  inputnames <- names(inputvec)
  for (currx0idx in seq_along(inputvec)) {
    currnames <- inputnames[currx0idx]
    tblx0$X0[tblx0$updatedstate == currnames] <- inputvec[[currnames]]
  }
  return(tblx0)
}
