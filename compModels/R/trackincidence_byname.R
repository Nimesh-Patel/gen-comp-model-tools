#' Update compiled model with incident populations aggregated by name
#'
#' @param compiledmodel list of instructions for piping |>
#' @param namevec vector of updated state names to aggregate
#' @param trackname character vector to name tracked population
#' @return updated compiled model with added populations
#' @export
trackincidence_byname <- function(compiledmodel, namevec, trackname = NA) {
  tblnames <- compiledmodel[["modelinstructions"]][["tblupdatedstates"]]
  peterstates <- compiledmodel[["modeloutstructions"]][["updatedstates"]]
  petermat <- compiledmodel[["modeloutstructions"]][["petermatrix"]]

  # check that names are in vector
  logic_checknames <- namevec %in% peterstates
  if (FALSE %in% logic_checknames) {
    stop("Named states to aggregate incidence not found in model. Check spelling
      and compare to names found in modelinstructions$tblupdatedstates.")
  }

  logic_aggregate <- peterstates %in% namevec
  # get named columns
  petermat_select <- petermat[logic_aggregate, , drop = FALSE]
  # focus on flow in (incidence)
  petermat_select[petermat_select <= 0] <- 0
  # aggregate all populations
  newpetermat <- rbind(petermat, Matrix::colSums(petermat_select))

  # get name
  if (is.na(trackname)) {
    logic_incnames <- grepl("^incidence_", peterstates)
    num_priorinc <- sum(logic_incnames)
    trackname <- paste0("incidence_", as.character(num_priorinc + 1))
    message(paste0(
      "No provided name for tracked incidence. Naming it ",
      trackname, "."
    ))
  }


  compiledmodel[["modeloutstructions"]][["updatedstates"]] <-
    c(peterstates, trackname)

  compiledmodel[["modeloutstructions"]][["petermatrix"]] <-
    newpetermat

  compiledmodel
}
