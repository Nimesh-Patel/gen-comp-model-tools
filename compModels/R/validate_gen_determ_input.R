#' Checks validity of input to generalized compartmental model
#'
#' Helper function that runs checks on generalized compartmental model input
#' parameters, times, and initial conditions
#'
#' @param init_vals the starting values for populations in each compartment
#' @param times vector of time steps to run the model over
#' @param comp_names names of the compartments
#' @param trans_matrix matrix of transition rates between each of the
#' @return stops with information or indicates Input Checks Passed
#' @examples
#' \dontrun{
#' validate_gen_determ_input(
#'   init_vals = c(1e05 - 1, 10, 5, 1, 0),
#'   times = seq(0.1, 100, by = 0.1),
#'   comp_names = c("A", "B", "C", "D", "E"),
#'   trans_matrix <- matrix(
#'     c(
#'       -0.1, 0.1, 0, 0, 0,
#'       0, -0.1, 0.1, 0, 0,
#'       0, 0, -0.1, 0.1, 0,
#'       0, 0, 0, -0.1, 0.1,
#'       0.1, 0, 0, 0, -0.1
#'     ),
#'     nrow = length(c("A", "B", "C", "D", "E")),
#'     byrow = TRUE,
#'     dimnames = list(c("A", "B", "C", "D", "E"), c("A", "B", "C", "D", "E"))
#'   )
#' )
#' }
validate_gen_determ_input <- function(init_vals, times,
                                      comp_names, trans_matrix) {
  if (length(init_vals) != length(comp_names)) {
    stop("The number of initial values must match the number of
    compartment names.")
  }

  if (!all(
    rownames(trans_matrix) %in% comp_names
  ) ||
    !all(colnames(trans_matrix) %in% comp_names)) {
    stop("Transition matrix row and column names must match the provided
    compartment names.")
  }
  print("Input checks passed")
}
