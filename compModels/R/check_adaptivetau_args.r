#' Validate adaptivetau argument values
#'
#' Helper function to validate adaptivetau argument values, used as part of a
#' validation function
#' @param transitions matrix object with transition information between states
#' @return stops with information or indicates checks passed
check_adaptivetau_args <- function(transitions) {
  if (!is.matrix(transitions)) stop("transitions must be a matrix")

  print("adaptivetau argument value checks passed")
}
