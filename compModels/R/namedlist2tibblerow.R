#' Create 1 row tibble from list where input
#' values are treated as lists
#'
#' @param namedlist named list
#' @return 1 row tibble with names as columns
#' and values as lists
namedlist2tibblerow <- function(namedlist) {
  tibblerow <- tibble::as_tibble(purrr::map(namedlist, function(x) {
    if (length(x) == 1) {
      x
    } else {
      list(x)
    }
  }))
  return(tibblerow)
}
