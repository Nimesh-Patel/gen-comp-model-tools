#' Convert list/vector to input length, keeps values as lists if
#'
#' Replicates length one vector and ensures n-length vectors
#' are same size as given length. respects if input
#' values are lists.
#'
#' @param inputvector character/numeric vector or list
#' @param desiredlength output vector length
#' @return lengthened vector/list same type as input
matchandpaste_basestate <- function(inputvector, desiredlength) {
  # coerces a vector to match a desired shape
  outputlist <- inputvector
  if (is.list(inputvector)) {
    if (length(inputvector) == 1) {
      outputlist <- rep(inputvector, desiredlength)
    }
  } else {
    outputlist <- rep(list(inputvector), desiredlength)
  }
  if (length(outputlist) != desiredlength) {
    stop("input base states must be vector or a list of
     length 1 or that of the vector to match")
  }
  return(outputlist)
}
