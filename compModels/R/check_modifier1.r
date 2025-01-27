#' Validate modifier 1
#'
#' Helper function to validate modifier values, used as part
#' of a validation function. A modifier should be one of: a single numeric
#' value, a matrix of numeric values, or a function.
#' @param modifier single value, matrix, or function to be applied to transition
#' matrix during an intervention
#' @return stops with information or indicates checks passed
#' }
check_modifier1 <- function(modifier) {
  is_null <- is.null(modifier)
  is_single_numeric <- is.numeric(modifier) && length(modifier) == 1
  is_numeric_matrix <- is.matrix(modifier) && is.numeric(modifier)
  is_function <- is.function(modifier)

  return(is_null || is_single_numeric || is_numeric_matrix || is_function)
  print("modifier check 1 passed")
}
