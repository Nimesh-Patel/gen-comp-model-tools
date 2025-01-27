#' Validate modifier 2
#'
#' Helper function to validate modifier values, used as part
#' of a validation function. A modifier should be one of: a single numeric
#' value, a matrix of numeric values, or a function.
#' @param modifier single value, matrix, or function to be applied to transition
#' matrix during an intervention
#' @return stops with information or indicates checks passed
#' }
check_modifier2 <- function(modifier) {
  if (!check_modifier1(modifier)) {
    stop("Modifier must be either a single numeric value, a matrix of numeric
          values, or a user-defined function.")
  }

  print("modifier checks 1 and 2 passed")
}
