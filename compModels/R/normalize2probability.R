#' normalize character/numeric vector to sum to 1
#'
#' @param inputvector character/numeric vector
#' @return normalized vector same type as input.
normalize2probability <- function(inputvector) {
  currsum <- sumvector(inputvector)
  if (is.numeric(inputvector)) {
    outputprobvector <- inputvector / currsum
  } else { # character
    outputprobvector <- paste0("((", inputvector, ")/", currsum, ")")
  }
  return(outputprobvector)
}
