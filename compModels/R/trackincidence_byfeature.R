#' Update compiled model with incident populations aggregated by feature
#'
#' @param compiledmodel list of instructions for piping |>
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
#' @param trackname character vector to name tracked population
#' @return updated compiled model with added populations
#' @export
trackincidence_byfeature <- function(compiledmodel,
                                     basestates = NA,
                                     groupnames = NA,
                                     metapopulation = NA,
                                     chains = NA,
                                     trackname = NA) {
  tblnames <- compiledmodel[["modelinstructions"]][["tblupdatedstates"]]
  peterstates <- compiledmodel[["modeloutstructions"]][["updatedstates"]]
  petermat <- compiledmodel[["modeloutstructions"]][["petermatrix"]]

  tblnames_filter <- filterpopsize_byfeature(tblnames,
    basestates = basestates,
    groupnames = groupnames,
    metapopulation = metapopulation,
    chains = chains
  )

  filternames <- tblnames_filter[["updatedstate"]]

  compiledmodel <- trackincidence_byname(
    compiledmodel,
    filternames,
    trackname
  )

  compiledmodel
}
