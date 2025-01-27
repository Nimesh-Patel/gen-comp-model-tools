#' Validate initial values
#'
#' Helper function to validate initial values, used as part of a validation
#' function
#' @param init_vals the starting values for populations in each compartment
#' @param comp_names names of the compartments
#' @return stops with information or indicates checks passed
#' }


check_init_vals <- function(init_vals, comp_names) {
  if (length(init_vals) != length(comp_names)) {
    stop("The number of initial values must match the number of compartment
          names.")
  }

  if (!all(is.numeric(init_vals))) {
    stop("All initial values must be numeric.")
  }

  if (any(init_vals < 0)) {
    stop("All initial values must be non-negative.")
  }

  print("initial values checks passed")
}
