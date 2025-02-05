#' Convert vector to input length
#'
#' Replicates length one vector and ensures n-length
#' vectors are same size as given length.
#' Same as matchandpaste.R
#'
#' @param vec character or numeric vector.
#' @param desiredlength output vector length
#' @return lengthened vector same type as input
matchlength <- function(vec, desiredlength) {
  if (length(vec) == 1) {
    vec <- rep(vec, desiredlength)
  }
  if (length(vec) != desiredlength) {
    stop("Can not coerce input vector.
    It's length is not 1 or the desired length.")
  }
  return(vec)
}
