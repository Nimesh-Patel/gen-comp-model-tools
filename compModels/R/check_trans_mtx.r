#' Validate transition matrix
#'
#' Helper function to validate transition matrix, used as part of a validation
#' function
#' @param param description
#' @return stops with information or indicates checks passed
#' }

check_trans_mtx <- function(trans_matrix, comp_names) {
  if (!all(
    rownames(trans_matrix) %in% comp_names
  ) ||
    !all(colnames(trans_matrix) %in% comp_names)) {
    stop("Transition matrix row and column names must match the provided
          compartment names.")
  }

  if (!is.matrix(trans_matrix) || !is.numeric(trans_matrix)) {
    stop("The transition matrix must be a numeric matrix.")
  }

  print("transition matrix checks passed")
}
