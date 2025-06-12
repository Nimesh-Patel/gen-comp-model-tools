#' Set popsize for states that match input features
#'
#' Match features based on input. Match using OR logic within each input
#' list/vector of basestates, groupnames, metapopulations, and chains. While
#' match using AND logic between input values.
#'
#' @param tblpopsize compiled updatedstate table with popsize value
#' @param newpopsize value to make all popsizes that match input features
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
setpopsize_byfeature <- function(tblpopsize,
                                 newpopsize,
                                 basestates = NA,
                                 groupnames = NA,
                                 metapopulation = NA,
                                 chains = NA) {
  # change names to avoid issues when calling columns
  check_inputpopsize(tblpopsize)

  if (length(newpopsize) > 1) {
    stop("newpopsize must be length 1.")
  }

  intlogic_before <- intlogic(tblpopsize)

  tblpopsize_filter <-
    filterpopsize_byfeature(tblpopsize,
      basestates = basestates,
      groupnames = groupnames,
      metapopulation = metapopulation,
      chains = chains
    )

  tblpopsize_filter <-
    tblpopsize_filter |>
    dplyr::mutate(popsizenew = newpopsize) |>
    dplyr::select("updatedstate", "popsizenew")

  tblout <- tblpopsize |>
    dplyr::left_join(tblpopsize_filter, by = "updatedstate") |>
    dplyr::mutate(popsize = dplyr::coalesce(.data$popsizenew, .data$popsize)) |>
    dplyr::select(-"popsizenew")

  intlogic_after <- intlogic(tblout)

  if (intlogic_before != intlogic_after) {
    message("popsize changed between non-integer and integer. Ensure this is
    intentional.")
  }

  tblout
}
