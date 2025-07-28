#' Replace "" in a vector with a given vector
#'
#' Return the given x vector replacing all occurrences of "" with the given
#' states vector. Ordering of "" replacement is not maintained since it's a
#' many-to-one replacement.
#'
#' @param x character vector.
#' @param states list of values to replace
#' @return character vector with replaced values
#' @family internal
empty2states <- function(x, states) {
  if ("" %in% x) {
    x <- c(x[x != ""], states)
  }
  return(x)
}
