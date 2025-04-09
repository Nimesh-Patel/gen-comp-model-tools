#' Modify state values
#'
#' Enables user to make fine grained updates to compartment state values.
#'
#' @param vec a named numeric vector containing the current values of states
#' to be modified
#' @param modifications (optional) a user supplied list of named state
#' compartments and corressponding new values for the named compartments,
#' default is empty list, list()
#' @param adjust_names a character string or vector of character strings of
#' the names or prefixes of the state compartments for which the user would
#' like to adjust the corresponding values
#' @param adjust_value a numeric value or vector of numeric values of
#' the revised values of the state compartments for which the user would like to
#' adjust the corresponding values
#' @return named vector of compartment state values
modify_states <- function(vec, modifications = list(),
                          adjust_names = NULL, adjust_value = 0) {
  for (name in names(modifications)) {
    if (name %in% names(vec)) {
      vec[name] <- modifications[[name]]
    } else {
      warning(paste("Name", name, "not found in state names."))
    }
  }
  if (!is.null(adjust_names) && !is.null(adjust_value)) {
    if (length(adjust_names) != length(adjust_value)) {
      stop("Length of adjust_names and adjust_value must be the same.")
    }
    for (i in seq_along(adjust_names)) {
      if (grepl("^", adjust_names[i])) { # Check if it's a prefix
        matching_names <- grep(paste0("^", adjust_names[i]),
          names(vec),
          value = TRUE
        )
        vec[matching_names] <- vec[matching_names] + adjust_value[i]
      } else { # It's a specific name
        if (adjust_names[i] %in% names(vec)) {
          vec[adjust_names[i]] <- vec[adjust_names[i]] + adjust_value[i]
        } else {
          warning(paste("Name", adjust_names[i], "not found in state names."))
        }
      }
    }
  }
  return(vec)
}
