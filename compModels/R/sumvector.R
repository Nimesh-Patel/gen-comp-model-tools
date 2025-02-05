#' Sum vector according to type
#'
#' Useful to normalize rates by total
#' population size when
#' states are input as numeric or character
#'
#'
#' @param inputvector numeric or character vector
#' @return 1 element numeric or character vector respecting input type
sumvector <- function(inputvector) {
  # sums a vector respecting its class, output is same class
  if (is.numeric(inputvector)) {
    currsum <- sum(inputvector)
  } else if (is.character(inputvector)) {
    currsum <- paste0("(", paste(inputvector, collapse = "+"), ")")
  } else {
    stop("Input vector to normalize to a probability distribution
    is neither numeric or a character.")
  }
  return(currsum)
}
