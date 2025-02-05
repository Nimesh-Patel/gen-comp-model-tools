#' Normalize vector to have mean 1
#'
#' @param inputvector Numeric or character vector.
#' @return normalized vector of same type
#'
normalize2mean1 <- function(inputvector) {
  currsum <- sumvector(inputvector)
  if (is.numeric(inputvector)) {
    renormfactor <- length(inputvector) / currsum
    outputvector <- renormfactor * inputvector
  } else { # character
    renormfactor <- paste0(
      "(",
      as.character(length(inputvector)),
      "/",
      currsum, ")"
    )
    outputvector <- paste0(renormfactor, "*", inputvector)
  }
  return(outputvector)
}
