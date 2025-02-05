#' Replace empty strings in a dataframe column with a given vector
#'
#'
#' @param currtbl A dataframe/tibble with columns
#' @param replacecol tibble name to replace "" values
#' @param replacevec character vector to replace with
#' @return tibble with specified column values ""
#' replaced with list of the vector to facilitate unnesting.
replaceglobal <- function(currtbl, replacecol, replacevec) {
  currtbl[[replacecol]][currtbl[[replacecol]] == ""] <- list(replacevec)
  return(currtbl)
}
