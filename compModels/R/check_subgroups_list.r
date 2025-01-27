#' Validate subgroups list
#'
#' Helper function to validate subgroups, used as part of a validation
#' function
#' @param param description
#' @return stops with information or indicates checks passed
#' }

check_subgroups_list <- function(subgroups_list) {
  if (!is.null(subgroups_list)) {
    if (!is.list(subgroups_list)) {
      stop("subgroups_list must be a list.")
    }

    valid_subgroup <-
      any(sapply(subgroups_list, function(x) is.vector(x) && length(x) >= 2))

    if (!valid_subgroup) {
      stop("subgroups_list must contain at least one element with at least
             two values")
    }
  }

  print("subgroups list checks passed")
}
