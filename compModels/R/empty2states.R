#' Expands instructions to apply globally
#'
#' Replace "" with given list. Write unit
#' test and check where used for difference
#' between vectors and list.
#'
#' @param x character vector.
#' @param states list of values to replace
#' @return character vector with replaced values
empty2states <- function(x, states) {
  if ("" %in% x) {
    x <- c(x[x != ""], states)
  }
  return(x)
}
