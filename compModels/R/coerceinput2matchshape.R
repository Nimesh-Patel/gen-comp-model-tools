#' Converts input list to nested lists of lists
#'
#' Useful when assigning multiple environments
#' and metapopulations
#'
#' @param queryinput character vector to be replicated
#' @param listofvectors list of lists whose shape is
#' use to coerce character vector
#' @return list of lists with inputbasestate as values
coerceinput2matchshape <- function(
    queryinput,
    listofvectors) {
  # user inputs different input classes.
  # repeat and flatten accordingly
  # get shapes of listofvectors which
  # the query input will be matched to
  numvectors <- length(listofvectors)
  listofvectors_vectorlength <- sapply(listofvectors, length)
  total_elements <- sum(listofvectors_vectorlength)
  # check that if list of same size as environment names
  if ((length(queryinput) != numvectors) && (length(queryinput) > 1)) {
    stop("Input list is not the same length and
         cannot be coerced into proper shape")
  }
  if (is.list(queryinput)) {
    if (length(queryinput) > 1) {
      queryinput_vectorlength <- sapply(queryinput, length)
      # check that element lengths are the same,
      # or coercible (i.e., replicate length 1 vectors)
      queryoutput <- unlist(
        mapply(matchandpaste, queryinput, listofvectors_vectorlength)
      )
    } else {
      # switch to vector
      queryinput <- as.vector(queryinput)
    }
  }
  if (!is.list(queryinput)) {
    queryoutput <- unlist(
      mapply(rep, queryinput,
        listofvectors_vectorlength,
        USE.NAMES = FALSE
      )
    )
  }
  return(queryoutput)
}
