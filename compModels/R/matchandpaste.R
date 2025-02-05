#' Convert vector to input length
#'
#' Replicates length one vector and ensures n-length
#' vectors are same size as given length.
#' Same as matchlength.R
#'
#' @param inputvector character or numeric vector.
#' @param desiredlength output vector length
#' @return lengthened vector same type as input
matchandpaste <- function(inputvector, desiredlength) {
  # coerces a vector to match a desired shape
  outputvector <- inputvector
  if (length(inputvector) == 1) {
    outputvector <- rep(inputvector, desiredlength)
  }
  if (length(outputvector) != desiredlength) {
    stop("input vector must be length 1 or same length as vector to match")
  }
  return(outputvector)
}
