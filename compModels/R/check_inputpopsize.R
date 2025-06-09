#' Ensure updatedstates input has popsize column
#'
#' @param input input to check
#' @return error if missing column

check_inputpopsize <- function(input) {
  if (!("popsize" %in% colnames(input))) {
    stop("Input/piped tblpopsize argument is missing a popsize column. Run
      define_popsize on a compiled model first to generate the table. ")
  }
}
