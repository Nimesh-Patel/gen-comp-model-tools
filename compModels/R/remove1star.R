#' Simplifies math equations
#'
#' Remove "1*" from string only if theres
#' no decimal before 1
#'
#' @param string Character vector to remove "1*" substring
#' @return string with "1*" removed
remove1star <- function(string) {
  outstring <- gsub("(?<!\\d)1\\*", "", string, perl = TRUE)
  return(outstring)
}
