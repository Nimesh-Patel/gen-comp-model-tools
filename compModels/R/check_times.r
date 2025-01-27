#' Validate time values
#'
#' Helper function to validate time values, used as part of a validation
#' function
#' @param param description
#' @return stops with information or indicates checks passed
#' }
check_times <- function(times) {
  if (length(times) < 1) {
    stop("The times vector must have at least one element.")
  }

  if (!is.numeric(times) || !all(diff(times) > 0)) {
    stop("The times vector must be numeric and strictly increasing.")
  }

  print("times values checks passed")
}
