#' Filter popsize table to match input features
#'
#' Match features based on input. Match using OR logic within each input
#' list/vector of basestates, groupnames, metapopulations, and chains. While
#' match using AND logic between input values.
#'
#' @param tblpopsize compiled updatedstate table with popsize value
#' @param basestates list/vector of basestates to match
#' default is NA which designates all basestates
#' @param groupnames named list/vector of groups to match. grouptypes are names
#' and groupnames are values.
#' default is NA which designates all groups
#' @param metapopulation list/vector of metapopulation names to match
#' default is NA which designates all metapopulations
#' @param chains list/vector of step numbers for chains to match. Names
#' correspond to chainid columns
#' default is NA which designates all chain steps
#' @return tibble with columns of updated state names
#' (updatedstates) and current conditions (popsize)
#' @importFrom rlang .data
filterpopsize_byfeature <- function(tblpopsize, # nolint: cyclocomp_linter.
                                    basestates = NA,
                                    groupnames = NA,
                                    metapopulation = NA,
                                    chains = NA) {
  # change names to avoid issues when calling columns
  basestates2 <- basestates
  metapopulation2 <- metapopulation
  # check input
  colnamescompiled <- colnames(tblpopsize)

  tblpopsize_filter <- tblpopsize
  # filter by each input
  # basestates
  if (!identical(basestates, NA)) {
    # check that it's got elements
    if (FALSE %in% (basestates2 %in% tblpopsize[["basestates"]])) {
      warning("Some input basestates are not found in table.
        These will be ignored, but check spelling")
    }
    tblpopsize_filter <-
      tblpopsize_filter |>
      dplyr::filter(.data$basestates %in% basestates2)
  }

  # metapopulation
  if (!identical(metapopulation, NA)) {
    # check that it's got elements
    if (FALSE %in% (metapopulation2 %in% tblpopsize$metapopulation)) {
      warning("Some input metapopulations are not found in table.
        These will be ignored, but check spelling")
    }
    tblpopsize_filter <-
      tblpopsize_filter |>
      dplyr::filter(.data$metapopulation %in% metapopulation2)
  }

  # groups
  if (!identical(groupnames, NA)) {
    # convert non-named list input to named list, grab column name
    groupnames <- fixinput_groupnames(groupnames, tblpopsize_filter)
    typenames <- names(groupnames)

    # warn if names not in colnames
    if (FALSE %in% (typenames %in% colnamescompiled)) {
      warning("Some input group types are not found in table.
        These will be ignored, but check spelling.")
    }

    for (currtype in typenames) {
      currgroups <- groupnames[[currtype]]
      if (currtype %in% colnamescompiled) {
        if (FALSE %in% (currgroups %in% tblpopsize_filter[[currtype]])) {
          warning(paste0(
            "Some input groups do not exist in the ",
            currtype,
            "type. These will be ignored, but check spelling."
          ))
        }
        currfilterlogic <- tblpopsize_filter[[currtype]] %in% currgroups
        tblpopsize_filter <- tblpopsize_filter[currfilterlogic, ]
      }
    }
  }

  # chains
  if (!identical(chains, NA)) {
    # coerece input into consistent named list form
    chains <- fixinput_chains(chains, tblpopsize_filter)

    chainnames <- names(chains)
    for (currchainname in chainnames) {
      currsteps <- chains[[currchainname]]
      if (currchainname %in% colnamescompiled) {
        if (FALSE %in% (currsteps %in% tblpopsize_filter[[currchainname]])) {
          warning(paste0(
            "Some chain steps do not exist for ",
            currchainname,
            ". These will be ignored, but check spelling."
          ))
        }
        currfilterlogic <- tblpopsize_filter[[currchainname]] %in% currsteps
        tblpopsize_filter <- tblpopsize_filter[currfilterlogic, ]
      } else {
        warning(paste0(currchainname, " is not found in input table and will
          be ignored when filtering. Check spelling."))
      }
    }
  }

  tblpopsize_filter
}
