#' Remove popsize for states that match input features
#'
#' Match features based on input. Match using OR logic within each input
#' list/vector of basestates, groupnames, metapopulations, and chains. While
#' match using AND logic between input values.
#'
#' @param tblpopsize compiled updatedstate table with popsize value
#' @param removepopsize value to make all popsizes that match input features
#' @param basestates list/vector of basestates to match
#' default is "" which designates all basestates
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
removepopsize_byfeature <- function(tblpopsize,
                                    removepopsize,
                                    basestates = NA,
                                    groupnames = NA,
                                    metapopulation = NA,
                                    chains = NA) {
  tblout <- addpopsize_byfeature(tblpopsize,
    -removepopsize,
    basestates = basestates,
    groupnames = groupnames,
    metapopulation = metapopulation,
    chains = chains
  )

  # check for negative values
  if (TRUE %in% (tblout[["popsize"]] < 0)) {
    stop("Removing population led to negative population.")
  }
  tblout
}
