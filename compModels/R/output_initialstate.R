#' Converts table of initial conditions to named list
#'
#' Useful as input to deSolve
#'
#' @param x0tbl tibble of initial conditions
#' with columns
#' x0 and updatedstate
#' @return named list with names updated state and values x0
output_initialstate <- function(x0tbl) {
  x0 <- x0tbl$X0
  names(x0) <- x0tbl$updatedstate
  return(x0)
}
