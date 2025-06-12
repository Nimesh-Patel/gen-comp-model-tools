#' Add popsize for states that match input features
#'
#' Add to popsize values for states that match input features. Match using OR
#' logic within each input list/vector of basestates, groupnames,
#' metapopulations, and chains. While match using AND logic between input
#' values.
#'
#' @param tblpopsize compiled updatedstate table with popsize value
#' @param addpopsize value to make all popsizes that match input features
#' @param basestates list/vector of basestates to match
#' default is NA which designates all basestates
#' @param groupnames named list/vector of groups to match. Names correspond to
#' grouptypes and the values correspond to groupnames. If names are not given
#' then grouptypes are automatically identified using tblpopsize.
#' default is NA which designates all groups
#' @param metapopulation list/vector of metapopulation names to match
#' default is NA which designates all metapopulations
#' @param chains list/vector of step numbers for chains to match. Names
#' correspond to chainid columns. If no chainid is given then it is
#' automatically identified from tblpopsize
#' default is NA which designates all chain steps
#' @return tibble with columns of updated state names
#' (updatedstates) and current conditions (popsize)
#' @export
#' @importFrom rlang .data
addpopsize_byfeature <- function(tblpopsize,
                                 addpopsize,
                                 basestates = NA,
                                 groupnames = NA,
                                 metapopulation = NA,
                                 chains = NA) {
  # change names to avoid issues when calling columns
  check_inputpopsize(tblpopsize)

  if (length(addpopsize) > 1) {
    stop("addpopsize must be length 1.")
  }

  intlogic_before <- intlogic(tblpopsize)

  tblpopsize_filter <-
    filterpopsize_byfeature(tblpopsize,
      basestates = basestates,
      groupnames = groupnames,
      metapopulation = metapopulation,
      chains = chains
    )

  for (rowidx in seq_len(nrow(tblpopsize_filter))) {
    currname <- tblpopsize_filter[["updatedstate"]][[rowidx]]
    currlogic <- tblpopsize[["updatedstate"]] == currname
    tblpopsize[["popsize"]][currlogic] <-
      tblpopsize[["popsize"]][currlogic] + addpopsize
  }

  intlogic_after <- intlogic(tblpopsize)

  if (intlogic_before != intlogic_after) {
    message("popsize changed between non-integer and integer. Ensure this is
    intentional.")
  }

  tblpopsize
}
