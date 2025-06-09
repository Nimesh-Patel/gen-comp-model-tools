#' Converts column to int if values have no remainder and total is larger than 1
#'
#' Useful for determining how to sample populations
#'
#' @param tbl2maintain table with column of values
#' @param col2maintain named column to convert values to int
#' Unspecified populations are set to 0
#' @return input tbl with column converted to int if values have no remainder
#' @importFrom rlang .data
maintainpopsizetype <- function(tbl2maintain, col2maintain) {
  outvec <- tbl2maintain[[col2maintain]]
  if (intlogic(outvec)) {
    tbl2maintain[[col2maintain]] <- as.integer(outvec)
  }
  tbl2maintain
}
