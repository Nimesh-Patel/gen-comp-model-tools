#' Converts given base state list to nested lists of lists.
#'
#' Useful when assigning multiple environments and metapopulations
#'
#' @param inputbasestate character vector to be replicated
#' @param listofvectors list of lists whose shape is
#' use to coerce character vector
#' @return list of lists with inputbasestate as values
coerceinput2matchshape_basestate <- function(inputbasestate, listofvectors) { # nolint: object_length_linter
  # user inputs different input classes. repeat
  # and flatten accordingly
  # get shapes of listofvectors which
  # the query input will be matched to
  numvectors <- length(listofvectors)
  listofvectors_vectorlength <- sapply(listofvectors, length)
  total_elements <- sum(listofvectors_vectorlength)
  # check that if list of same size as environment names
  if (is.list(inputbasestate)) {
    basestate_output <- unlist(mapply(
      matchandpaste_basestate, inputbasestate, listofvectors_vectorlength
    ), recursive = FALSE)
  } else {
    basestate_output <- rep(list(inputbasestate), total_elements)
  }
  return(basestate_output)
}
